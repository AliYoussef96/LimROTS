#' LimROTS Microbiome Differential-Abundance Analysis
#'
#' Analyze microbiome count data using prevalence filtering, compositional
#' preprocessing, and LimROTS bootstrap and permutation statistics.
#'
#' @param x A `SummarizedExperiment` containing feature counts in rows and
#' samples in columns.
#' @param assay.type Character string or numeric index selecting the assay in
#' `x` to analyze. If `NULL`, the first assay is used.
#' @param niter Integer. Number of bootstrap or permutation replicates to
#' generate.
#' @param K Integer. Number of top-ranked features to use. If `NULL`, defaults
#' to one quarter of the number of features.
#' @param meta.info Character vector naming the sample-level metadata columns
#' in `colData(x)` required by the model.
#' @param BPPARAM A `BiocParallelParam` object specifying the parallel backend.
#' If `NULL`, the parallel helper selects a two-worker backend appropriate to
#' the operating system.
#' @param group Character string naming the column in `meta.info` that
#' identifies each sample's group.
#' @param formula.str Character string containing the model formula. It must
#' include the group and any covariates to model.
#' @param robust Logical indicating whether to use robust empirical Bayes
#' moderation. See \link[limma]{eBayes}.
#' @param trend Logical indicating whether to fit an intensity-dependent trend
#' in the prior variance. See \link[limma]{eBayes}.
#' @param pairwise_contrasts Character vector of contrast expressions, such as
#' `"group1 - group2"`, to fit.
#' @param model Character string specifying the fitting method, such as
#' `"limma"` or `"dream"`.
#' @param strata Character vector or `NULL`. Names of character or factor
#' metadata columns used to define strata for resampling and permutation.
#' @param prev_filter Numeric scalar specifying the minimum prevalence required
#' for a feature to be retained.
#' @param is_winsor Logical indicating whether to apply winsorization during
#' microbiome preprocessing.
#' @param outlier_pct Numeric scalar specifying the upper proportion of values
#' to winsorize within each feature.
#' @param pseudo_count Positive numeric scalar added to processed counts before
#' the log transformation.
#'
#' @importFrom mia subsetByPrevalent transformAssay
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @importFrom SummarizedExperiment assay colData
#' @importFrom S4Vectors DataFrame
#' @importFrom DescTools Winsorize
#' @importFrom BiocParallel SnowParam MulticoreParam bplapply bpoptions
#' @importFrom limma makeContrasts lmFit contrasts.fit eBayes
#' @importFrom modeest mlv
#' @importFrom stats as.formula complete.cases formula median model.matrix
#' setNames quantile
#' @importFrom stringr fixed str_split_fixed
#' @importFrom utils combn
#'
#' @return The input `SummarizedExperiment` with feature-level statistics,
#' empirical p-values, q-values, FDR estimates, corrected log-fold changes,
#' and BH-adjusted p-values added to `rowData()`. Optimization parameters and
#' null statistics are added to `metadata()`.
#'
#' @details The microbiome compositional bias correction is adapted from the
#' LinDA method implemented in \pkg{MicrobiomeStat}. LimROTS combines this
#' correction with its bootstrap and permutation-based inference.
#'
#' @seealso \code{\link[MicrobiomeStat]{linda}}
#'
#' @references
#' Huijuan Zhou, Kejun He, Jun Chen, and Xianyang Zhang. (2021). LinDA: Linear
#' Models for Differential Abundance Analysis of Microbiome Compositional Data.
#'
#' @export
LimROTS_microb <- function(
    x,
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
        model = "limma",
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
        model = "limma",
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


#' Preprocess Microbiome Counts for LimROTS
#'
#' Filter features by prevalence and prepare transformed microbiome data and
#' aligned sample metadata for LimROTS analysis.
#'
#' @param tse A `TreeSummarizedExperiment` containing non-negative feature
#' counts in the `"counts"` assay.
#' @param prev_filter Numeric scalar specifying the minimum prevalence required
#' for a feature to be retained.
#' @param is_winsor Logical indicating whether to apply winsorization during
#' microbiome preprocessing.
#' @param outlier_pct Numeric scalar specifying the upper proportion of values
#' to winsorize within each feature.
#' @param pseudo_count Positive numeric scalar added to processed counts before
#' the log transformation.
#'
#' @return A list containing the transformed feature data `clr_features` and
#' the aligned sample metadata `metadata_use`.
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

        relative_abundance <- as.matrix(
            SummarizedExperiment::assay(tse, "relabundance")
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


#' Fit Final Microbiome Models and Calculate Inference
#'
#' Fit the model on the original microbiome data, scale observed and null
#' statistics, and calculate empirical p-values and false-discovery measures.
#'
#' @param data Feature-by-sample matrix with sample names as column names.
#' @param meta.info Data frame of sample-level metadata, with sample names as
#' row names.
#' @param group Character string naming the column in `meta.info` that
#' identifies each sample's group.
#' @param formula.str Character string containing the model formula.
#' @param trend Logical indicating whether to fit an intensity-dependent trend
#' in the prior variance. See \link[limma]{eBayes}.
#' @param robust Logical indicating whether to use robust empirical Bayes
#' moderation. See \link[limma]{eBayes}.
#' @param a1 Numeric first statistic-scaling parameter.
#' @param a2 Numeric second statistic-scaling parameter.
#' @param pairwise_contrasts Character vector of contrast expressions, such as
#' `"group1 - group2"`, to fit.
#' @param model Character string specifying the fitting method, such as
#' `"limma"` or `"dream"`.
#' @param pD Matrix of null-model statistics from the permutation fits.
#' @param pS Matrix of null-model standard errors from the permutation fits.
#'
#' @return A list containing scaled observed statistics `d`, scaled null
#' statistics `pD`, empirical p-values `p_value`, FDR estimates `FDR`,
#' corrected log-fold changes `corrected.logfc`, q-values `q_values`, and
#' BH-adjusted p-values `BH.pvalue`.
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
