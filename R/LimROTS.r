# Copyright (C) 2025 Ali Youssef.
# This program is licensed under GPL (>= 2).
#
#'LimROTS: A Hybrid Method Integrating Empirical Bayes and 
#'Reproducibility-Optimized Statistics for Robust Differential Expression Analysis
#' 
#' Analyze feature-by-sample data using linear models, empirical Bayes
#' moderation, and reproducibility-optimized bootstrap statistics.
#'
#' @param x A `SummarizedExperiment` containing features in rows and samples in
#' columns.
#' @param assay.type Character string  selecting the assay in
#' `x` to analyze. If `NULL`, the first assay is used.
#' @param niter Integer. Number of bootstrap or permutation replicates to
#' generate.
#' @param K Integer. Number of top-ranked features to use. If `NULL`, defaults
#' to one quarter of the number of features.
#' @param meta.info Character vector naming the sample-level metadata columns
#' in `colData(x)` required by the model.
#' @param group Character string naming the column in `meta.info` that
#' identifies each sample's group.
#' @param BPPARAM A `BiocParallelParam` object specifying the parallel backend.
#' If `NULL`, the parallel helper selects a two-worker backend appropriate to
#' the operating system.
#' @param formula.str Character string containing the model formula. It must
#' include the group and any covariates to model.
#' @param robust Logical indicating whether to use robust empirical Bayes
#' moderation. See \link[limma]{eBayes}.
#' @param trend Logical indicating whether to fit an intensity-dependent trend
#' in the prior variance. See \link[limma]{eBayes}.
#' @param correlation_block Character string or `NULL`. Name of the column in
#' `meta.info` identifying correlated samples or blocks. When supplied, blocks
#' are preserved during resampling and within-block correlation is modeled.
#' @param pairwise_contrasts Character vector of contrast expressions, such as
#' `"group1 - group2"`, to fit.
#' @param model Character string specifying the fitting method: `"limma"` for
#' fixed effects or `"dream"` for random effects.
#' @param strata Character vector or `NULL`. Names of character or factor
#' metadata columns used to define strata for resampling and permutation.
#'
#' @return The input `SummarizedExperiment` with feature-level statistics,
#' empirical p-values, q-values, FDR estimates, corrected log-fold changes,
#' and BH-adjusted p-values added to `rowData()`. Optimization parameters and
#' null statistics are added to `metadata()`.
#'
#' @examples
#' # Example usage:
#'
#' data <- data.frame(matrix(rnorm(500), nrow = 100, ncol = 10))
#' # Simulated data
#' meta.info <- data.frame(
#'     group = factor(rep(1:2, each = 5)),
#'     row.names = colnames(data)
#' )
#' formula.str <- "~ 0 + group"
#' result <- LimROTS(data,
#'     meta.info = meta.info, group = "group",
#'     formula.str = formula.str, niter = 10
#' )
#'
#' @importFrom stats model.matrix formula p.adjust
#' @importFrom dplyr bind_cols
#' @importFrom qvalue empPvals qvalue
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom BiocParallel SnowParam MulticoreParam bplapply bpoptions
#' bpRNGseed
#' @importFrom S4Vectors DataFrame metadata
#'
#' @details LimROTS fits linear models with \pkg{limma} or dream and applies empirical
#' Bayes moderation to estimate feature coefficients and standard errors.
#' Reproducibility-optimized parameters \eqn{\alpha_1} and \eqn{\alpha_2} are
#' selected by maximizing the overlap of top-ranked features across
#' group-preserving bootstrap datasets, following the \link[ROTS]{ROTS}
#' approach. The final statistic for each feature is:
#'
#' \deqn{t_{\alpha_{(p)}} = \frac{\beta_{(p)}}
#'          {\alpha_1 + \alpha_2 \times SE_{post(p)}}}
#'
#' where \eqn{\beta_{(p)}} is the model coefficient and
#' \eqn{SE_{post(p)}} is its empirical Bayes-adjusted standard error.
#' Empirical p-values are calculated from permutation-based null statistics.
#' The function also estimates FDR, q-values, and Benjamini-Hochberg adjusted
#' p-values.
#'
#' This function processes a dataset using parallel computation. It
#' leverages the \pkg{BiocParallel} framework to distribute tasks
#' across multiple workers, which  can significantly reduce runtime for
#' large datasets.
#'
#' @references
#'   Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W., and Smyth,
#'   G.K. (2015). limma powers differential expression analyses for
#'   RNA-sequencing and microarray studies. Nucleic Acids Research 43(7), e47.
#'
#'   Suomi T, Seyednasrollah F, Jaakkola M, Faux T, Elo L (2017). “ROTS: An
#'   R package for reproducibility-optimized statistical testing.
#'   ” _PLoS computational biology_, *13*(5), e1005562.
#'   \url{doi:10.1371/journal.pcbi.1005562}
#'   \url{https://doi.org/10.1371/journal.pcbi.1005562},
#'   \url{http://www.ncbi.nlm.nih.gov/pubmed/28542205}
#'
#'   Elo LL, Filen S, Lahesmaa R, Aittokallio T. Reproducibility-optimized test
#'   statistic for ranking genes in microarray studies.
#'   IEEE/ACM Trans Comput Biol Bioinform. 2008;5(3):423-431.
#'   \url{doi:10.1109/tcbb.2007.1078}
#'
#' @export


LimROTS <- function(x,
    assay.type = NULL,
    niter = 1000,
    K = NULL,
    meta.info,
    BPPARAM  = NULL,
    group,
    formula.str,
    robust = TRUE,
    trend = TRUE,
    correlation_block = NULL,
    pairwise_contrasts,
    model,
    strata
    ) {
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
        correlation_block = correlation_block,
        strata = strata,
        survival = FALSE,
        competing_risks = FALSE
    )


    meta.info <- SanityChecK.list[["meta.info"]]
    data <- SanityChecK.list[["data"]]
    group <- SanityChecK.list[["group"]]
    K <- SanityChecK.list[["K"]]

    ### bootstrap and permutation samples
    list_boot_perm <- .generate_bootstrap_permutations(
        niter_2 = 2 * niter,
        meta.info = meta.info,
        group = group,
        correlation_block = correlation_block,
        strata = strata
    )

    ### parallelized bootstrapping and permutation fit
    results_list <- Boot_parallel(
        BPPARAM  = BPPARAM ,
        samples = list_boot_perm[["boot_Samples"]], 
        data = data,
        formula.str = formula.str,
        group = group,
        meta.info = meta.info,
        pSamples = list_boot_perm[["perm_Samples"]],
        correlation_block = correlation_block,
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

    ### optimize parameters
    optimized.parameters <-
        Optimizing(niter =  niter, K = K,
        N = N, D = D, S = S, pD = pD, pS = pS)
    a1 <- optimized.parameters[["a1"]]
    a2 <- optimized.parameters[["a2"]]
    k <- optimized.parameters[["k"]]
    R <- optimized.parameters[["R"]]
    Z <- optimized.parameters[["Z"]]
    ztable <- optimized.parameters[["ztable"]]

    ### fit limma on the original data

    fit_limma <- .final_fit(
        data = data,
        meta.info = meta.info,
        group = group,
        formula.str = formula.str,
        trend = trend,
        robust = robust,
        correlation_block = correlation_block,
        a1 = a1, a2 = a2,
        pairwise_contrasts = pairwise_contrasts,
        model = model,
        pD = pD, pS = pS
    )
    d = fit_limma[["d"]]
    #pD = fit_limma[["pD"]]
    p_value = fit_limma[["p_value"]]
    FDR = fit_limma[["FDR"]]
    corrected.logfc = fit_limma[["corrected.logfc"]]
    q_values = fit_limma[["q_values"]]
    BH.pvalue = fit_limma[["BH.pvalue"]]

    ### generate output

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
    
    return(LimROTS.output)
}



#' Generate Bootstrap and Permutation Samples
#'
#' Generate group-preserving bootstrap samples and corresponding permuted
#' metadata, using optional correlation blocks or strata.
#'
#' @param niter_2 Integer. Number of bootstrap and permutation replicates to
#' generate.
#' @param meta.info Data frame of sample-level metadata, with sample names as
#' row names.
#' @param group Character string naming the column in `meta.info` that
#' identifies each sample's group.
#' @param correlation_block Character string or `NULL`. Name of the column in
#' `meta.info` identifying correlated samples or blocks.
#' @param strata Character vector or `NULL`. Names of character or factor
#' metadata columns used to define strata for resampling and permutation.
#'
#' @return A list containing `boot_Samples`, a matrix of bootstrap sample
#' names, and `perm_Samples`, a list of permuted metadata data frames.
.generate_bootstrap_permutations <- function(niter_2, meta.info, group, 
                                             correlation_block, strata) {
    

    if(!is.null(correlation_block)){
        
        message("..........")
        
        samples <- bootstrapSamples_limRots_block(
            niter = niter_2,
            meta.info = meta.info,
            group = group,
            correlation_block = correlation_block
        )
        
        pSamples <- permutationSamples_limRots(niter = niter_2, meta.info = meta.info, 
                                   group = group, strata = correlation_block)
    
        return(list(boot_Samples = samples, perm_Samples = pSamples))
        
    }else if(is.null(correlation_block)){
        message("..........")
        
        if(!is.null(strata)){
            message("..........")
            
            samples <- bootstrapSamples_limRots(
                niter = niter_2,
                meta.info = meta.info,
                group = group, strata = strata
            )
            
            pSamples <- permutationSamples_limRots(niter = niter_2, meta.info = meta.info, 
                                                   group = group, strata = strata)
            
            return(list(boot_Samples = samples, perm_Samples = pSamples))
            
        }else if(is.null(strata)){
            message("..........")
            
            samples <- bootstrapS(niter = niter_2, meta.info = meta.info, group = group)

            pSamples <- vector("list", niter_2)
            for (i in seq_len(niter_2)) {
                shuffle_df <- meta.info
                ## Shuffle only the target group labels across all samples
                shuffle_df[[group]] <- sample(
                    meta.info[[group]],
                    size = nrow(meta.info),
                    replace = FALSE
                )
                pSamples[[i]] <- shuffle_df
            }
            
            return(list(boot_Samples = samples, perm_Samples = pSamples))
            
        }
        
    }
    
}

#' Fit Original Data and Calculate Inference
#'
#' Fit the selected linear model to the original data, scale observed and null
#' statistics, and calculate empirical p-values and false-discovery measures.
#'
#' @param data Feature-by-sample matrix with sample names as column names.
#' @param meta.info Data frame of sample-level metadata, with sample names as
#' row names.
#' @param group Character string naming the column in `meta.info` that
#' identifies each sample's group.
#' @param formula.str Character string containing the model formula.
#' @param trend Logical indicating whether to fit an intensity-dependent trend
#' in the prior variance.
#' @param robust Logical indicating whether to use robust empirical Bayes
#' moderation.
#' @param correlation_block Character string or `NULL`. Name of the column in
#' `meta.info` identifying correlated samples or blocks.
#' @param a1 Numeric first statistic-scaling parameter.
#' @param a2 Numeric second statistic-scaling parameter.
#' @param pairwise_contrasts Character vector of contrast expressions, such as
#' `"group1 - group2"`, to fit.
#' @param model Character string specifying the fitting method: `"limma"` for
#' fixed effects or `"dream"` for random effects.
#' @param pD Matrix of null-model statistics from the permutation fits.
#' @param pS Matrix of null-model standard errors from the permutation fits.
#'
#' @return A list containing scaled observed statistics `d`, scaled null
#' statistics `pD`, empirical p-values `p_value`, FDR estimates `FDR`,
#' corrected log-fold changes `corrected.logfc`, q-values `q_values`, and
#' BH-adjusted p-values `BH.pvalue`.
.final_fit <- function(data, meta.info, group, formula.str, trend, robust,
    correlation_block, a1, a2, pairwise_contrasts, model,
    pD , pS) {
        # Fit the limma model on the original data and compute statistics, p-values, and FDR.
        fit <- fit_original(
            x = data,
            group = group,
            meta.info = meta.info,
            formula.str = formula.str,
            trend = trend,
            robust = robust,
            correlation_block = correlation_block,
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

#' Add LimROTS Results to a SummarizedExperiment
#'
#' Add feature-level statistics to `rowData()` and optimization parameters and
#' null statistics to `metadata()`.
#'
#' @param x A `SummarizedExperiment` whose rows correspond to analyzed features.
#' @param data Feature-by-sample matrix used for fitting, with feature names as
#' row names.
#' @param d Numeric vector of scaled feature statistics.
#' @param p_value Numeric vector of empirical p-values.
#' @param FDR Numeric vector of estimated false discovery rates.
#' @param corrected.logfc Numeric vector of signed model coefficients.
#' @param q_values Q-value estimation result from `qvalue()`, or `NULL` if
#' estimation failed.
#' @param BH.pvalue Numeric vector of Benjamini-Hochberg adjusted p-values.
#' @param a1 Numeric first statistic-scaling parameter.
#' @param a2 Numeric second statistic-scaling parameter.
#' @param k Integer. Number of top-ranked features selected during optimization.
#' @param Z Numeric optimization z-score.
#' @param R Numeric reproducibility score from optimization.
#' @param ztable Matrix of z-scores across the optimization parameter grid.
#' @param pD Matrix of scaled null statistics from the permutation fits.
#'
#' @return The `SummarizedExperiment` `x` with differential-expression results
#' added to its row data and metadata.
.generate_output <- function(x, data, d, p_value, FDR, corrected.logfc,
    q_values, BH.pvalue, a1, a2, k, Z, R, ztable, pD) {
        # Generate the output SummarizedExperiment object with the LimROTS results.
    new_rowData <- DataFrame(
        statistics = d,
        pvalue = p_value,
        qvalue = if (is.null(q_values)) rep(NA_real_, length(d))
                    else q_values$qvalues,
        FDR = FDR,
        corrected.logfc = corrected.logfc,
        BH.pvalue = BH.pvalue,
        row.names = row.names(data)
    )

    new_rowData <- new_rowData[match(rownames(x), 
                                            rownames(new_rowData)), ]
    if (!identical(rownames(new_rowData), rownames(x))) {
        stop("Can not add the LimROTS results to the 
                                                SummarizedExperiment")
    }
    correct.order <- rownames(x)
    rowData(x) <- cbind(rowData(x), new_rowData)
    if (!identical(correct.order, rownames(x))) {
        stop("Can not add the LimROTS results to the 
                                                SummarizedExperiment")
    }
    metadata(x) <- c(metadata(x), list(
        a1 = a1,
        a2 = a2,
        k = k,
        Z = Z,
        R = R,
        ztable = ztable,
        q_values = q_values,
        null.statistics = pD
    ))

    LimROTS.output <- x
    return(LimROTS.output)
    }
