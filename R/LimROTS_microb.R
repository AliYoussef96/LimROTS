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

    row_length <- nrow(as.matrix(data))
    col_length <- nrow(list_boot_perm[["boot_Samples"]])

    D <- matrix(nrow = row_length, ncol = col_length)
    S <- matrix(nrow = row_length, ncol = col_length)
    pD <- matrix(nrow = row_length, ncol = col_length)
    pS <- matrix(nrow = row_length, ncol = col_length)

    for (i in seq_along(results_list)) {
        D[, i] <- results_list[[i]][["ds"]][["d_result"]]
        S[, i] <- results_list[[i]][["ds"]][["s_result"]]
        pD[, i] <- results_list[[i]][["pdps"]][["pd_result"]]
        pS[, i] <- results_list[[i]][["pdps"]][["ps_result"]]
    }

    if (any(vapply(list(D, S, pD, pS), anyNA, logical(1)))) {
        stop("D, S, pD, and pS must not contain NA values.")
    }
    
    if (!all(vapply(list(S, pD, pS), function(x) identical(dim(x), dim(D)), logical(1)))) {
        stop("D, S, pD, and pS must all have identical dimensions.")
    }

    if (any(vapply(list(D, S, pD, pS), function(x) any(!is.finite(x)), logical(1)))) {
    stop("D, S, pD, and pS must not contain NA, NaN, Inf, or -Inf.")
    }
    
    if (any(vapply(list(D, S, pD, pS), function(x) any(x < 0), logical(1)))) {
    stop("D, S, pD, and pS must contain only non-negative values.")
    }

    ### optimize parameters
    optimized.parameters <-
        Optimizing(niter = niter, K = K, D = D, S = S, pD = pD, pS = pS)
    a1 <- optimized.parameters[["a1"]]
    a2 <- optimized.parameters[["a2"]]
    k <- optimized.parameters[["k"]]
    R <- optimized.parameters[["R"]]
    Z <- optimized.parameters[["Z"]]
    ztable <- optimized.parameters[["ztable"]]

    fit_microb_final <- .final_fit_microb(
        data = data,
        meta.info = meta.info,
        group = group,
        formula.str = formula.str,
        trend = trend,
        robust = robust,
        a1 = a1, a2 = a2,
        pairwise_contrasts = pairwise_contrasts,
        model = model,
        pD = pD, pS = pS
    )

    d = fit_microb_final[["d"]]
    #pD = fit_microb_final[["pD"]]
    p_value = fit_microb_final[["p_value"]]
    FDR = fit_microb_final[["FDR"]]
    corrected.logfc = fit_microb_final[["corrected.logfc"]]
    q_values = fit_microb_final[["q_values"]]
    BH.pvalue = fit_microb_final[["BH.pvalue"]]

    ### generate output
    ### this from LimROTS main function file
    LimROTS.output <- .generate_output(
        x = x,
        data = data,
        d = d,
        p_value = p_value,
        FDR = FDR,
        corrected.logfc = corrected.logfc,
        q_values = q_values,
        BH.pvalue = BH.pvalue,
        a1 = a1, a2 = a2, k = k, Z = Z, R = R, ztable = ztable, pD = pD
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


.final_fit <- function(data, meta.info, group, formula.str, trend, robust,
 a1, a2, pairwise_contrasts, model,
    pD , pS) {
        # Fit the limma model on the original data and compute statistics, p-values, and FDR.
        fit <- fit_original_microb(
            x = data,
            group = group,
            meta.info = meta.info,
            formula.str = formula.str,
            trend = trend,
            robust = robust,
            pairwise_contrasts = pairwise_contrasts,
            model = model
        )

        d <- fit$d / (a1 + a2 * fit$s)
        pD <- pD / (a1 + a2 * pS)

        p_value <- empPvals(
            stat = d,
            stat0 = pD,
            pool = TRUE
        )
        FDR <- calculateFalseDiscoveryRate(d, pD)
        corrected.logfc <- fit$corrected.logfc
        q_values <- tryCatch(
            {
                    qvalue(
                        p_value,
                        pi0.method = "bootstrap",
                        lambda = seq(0.01, 0.95, 0.01))
            },
            error = function(e) {
                message("qvalue() failed (return NULL): ", e$message)
                NULL
            }
        )
        BH.pvalue <- p.adjust(p_value, method = "BH")
    return(list(
        d = d,
        pD = pD,
        p_value = p_value,
        FDR = FDR,
        corrected.logfc = corrected.logfc,
        q_values = q_values,
        BH.pvalue = BH.pvalue
    ))

    }
