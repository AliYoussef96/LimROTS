LimROTS_microb <- function(
    assay.type = NULL,
    niter = 1000,
    K = NULL,
    meta.info,
    BPPARAM  = NULL,
    group,
    formula.str,
    robust = TRUE,
    trend = TRUE,
    pairwise_contrasts,
    model,
    strata,
    prev_filter = 0,
    is_winsor = TRUE,
    outlier_pct = 0.03,
    pseudo_count = 0.5){
    ### Perform sanity checks on the input parameters and data.
    SanityChecK.list <- SanityChecK(
        x,
        assay.type = assay.type,
        niter = niter,
        K = K,
        meta.info = meta.info,
        group = group,
        model = model,
        formula.str = formula.str,
        pairwise_contrasts = pairwise_contrasts,
        correlation_block = NULL,
        strata = strata,
        survival = FALSE,
        competing_risks = FALSE
    )

    meta.info <- SanityChecK.list[["meta.info"]]
    data <- SanityChecK.list[["data"]]

    ### Extra Perform sanity checks on the data

    if (anyNA(data) || any(!is.finite(data))) {
        stop("data contains missing or non-finite values.")
    }

    if (any(data < 0)) {
        stop("data cannot contain negative values.")
    }

    if (any(data != round(data))) {
        stop("data must contain original integer counts.")
    }

    tse <- TreeSummarizedExperiment(
        assays = list(counts = data),
        colData = DataFrame(meta.info)
    )

    pre_microb <- .preprocess_microb(tse = tse, prev_filter = prev_filter, 
    is_winsor = is_winsor, outlier_pct = outlier_pct, pseudo_count = pseudo_count)

    
    SanityChecK.list <- SanityChecK(
        SummarizedExperiment(assays = list(clr = pre_microb[["clr_features"]]),
        colData = DataFrame(pre_microb[["metadata_use"]])),
        assay.type = assay.type,
        niter = niter,
        K = K,
        meta.info = colnames(pre_microb[["metadata_use"]]),
        group = group,
        model = model,
        formula.str = formula.str,
        pairwise_contrasts = pairwise_contrasts,
        correlation_block = NULL,
        strata = strata,
        survival = FALSE,
        competing_risks = FALSE
    )

    meta.info <- SanityChecK.list[["meta.info"]]
    data <- SanityChecK.list[["data"]]
    group <- SanityChecK.list[["group"]]
    K <- SanityChecK.list[["K"]]

    ### bootstrap and permutation samples (from the LimROTS.R script)
    list_boot_perm <- .generate_bootstrap_permutations(
        niter_2 = 2 * niter,
        meta.info = meta.info,
        group = group,
        correlation_block = NULL,
        strata = strata
    )

    ### parallelized bootstrapping and permutation fit
    results_list <- parallel_microb(
        BPPARAM  = BPPARAM ,
        samples = list_boot_perm[["boot_Samples"]], 
        data = data,
        formula.str = formula.str,
        group = group,
        meta.info = meta.info,
        pSamples = list_boot_perm[["perm_Samples"]],
        model = model,
        pairwise_contrasts = pairwise_contrasts
    )


}


.preprocess_microb <- function(tse, prev_filter, is_winsor, outlier_pct, pseudo_count) {

    tse <- subsetByPrevalent(
            tse,
            assay.type = "counts",
            detection = 0,
            prevalence = prev_filter,
            include.lowest = TRUE,
            update.tree = FALSE
        )

        if (nrow(tse) == 0L) {
            stop("No features remain after prevalence filtering.")
        }

        filtered_counts <- as.matrix(
            SummarizedExperiment::assay(tse, "counts")
        )


        library_sizes <- colSums(filtered_counts)

        tse <- mia::transformAssay(
            tse,
            assay.type = "counts",
            method = "relabundance",
            name = "relabundance",
            pseudocount = FALSE
        )


        if (is_winsor) {
            if (!is.numeric(outlier_pct) || length(outlier_pct) != 1L ||
                is.na(outlier_pct) || outlier_pct < 0 || outlier_pct >= 1) {
                stop("outlier_pct must be a single number in [0, 1).")
            }

            feature_names <- rownames(relative_abundance)
            sample_names <- colnames(relative_abundance)

            relative_abundance <- t(
                apply(relative_abundance, 1, function(x) {
                    upper_limit <- quantile(
                        x,
                        probs = 1 - outlier_pct,
                        names = FALSE,
                        na.rm = TRUE,
                        type = 7
                    )

                    Winsorize(
                        x,
                        val = c(min(x, na.rm = TRUE), upper_limit)
                    )
                })
            )

            dimnames(relative_abundance) <- list(
                feature_names,
                sample_names
            )

            processed_counts <- sweep(
                relative_abundance,
                2,
                library_sizes,
                "*"
            )

            processed_counts <- round(processed_counts)
        } else {
            processed_counts <- filtered_counts
        }

        assay(
            tse,
            "winsorized_counts"
        ) <- processed_counts

        tse <- transformAssay(
            tse,
            assay.type = "winsorized_counts",
            method = "pseudocount",
            name = "processed_counts",
            pseudocount = pseudo_count
        )

        processed_counts <- as.matrix(
            assay(tse, "processed_counts")
        )

        if (any(processed_counts <= 0)) {
            stop("All processed counts must be positive before log2 transformation.")
        }

        log_counts <- log2(processed_counts)
        clr_samples_by_features <- t(log_counts) - colMeans(log_counts)
        clr_samples_by_features <- data.frame(t(clr_samples_by_features), 
        check.names = FALSE,
        check.rows = FALSE)

        metadata_use <- as.data.frame(
            colData(tse)
        )

        metadata_use <- metadata_use[
            colnames(clr_samples_by_features),
            ,
            drop = FALSE
        ]

        if (!identical(colnames(clr_samples_by_features), rownames(metadata_use))) {
            stop("Processed counts and metadata are not aligned.")
        }


        return(list(
            clr_features = clr_samples_by_features,
            metadata_use = metadata_use
        ))

}