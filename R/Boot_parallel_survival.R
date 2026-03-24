# Copyright (C) 2025 Ali Youssef.
# This program is licensed under GPL (>= 2).
#
#' Parallel processing handling function for LimROTS survival
#'
#' @param data A \code{SummarizedExperiment} object or a matrix where rows
#' represent features (e.g., genes, proteins) and columns represent samples.
#' The values should be log-transformed.
#' @param a1 Optional numeric value used in the optimization process.
#' If defined by the user, no optimization occurs.
#' @param a2 Optional numeric value used in the optimization process.
#' If defined by the user, no optimization occurs.
#' @param meta.info A data frame containing sample-level metadata, where each
#' row corresponds to a sample. It should include the grouping variable
#' specified in \code{group.name}. If \code{x} is a \code{SummarizedExperiment}
#' object, \code{meta.info} must be a vector of the metadata needed for the
#' model to run and can be retrieved using \code{colData()}.
#' @param BPPARAM A parallel BPPARAM object for distributed computation.
#' @param formula.str A formula string used when covariates are present in meta.
#' info for modeling. It should include "~ 0 + ..." to exclude the
#' intercept from the model.
#' @param samples bootstrapped samples matrix
#' @param groups groups information from `meta.info`
#' @param pSamples a permutated list of samples
#' @param competing_risks Logical. If \code{TRUE}, the Fine\u2013Gray competing
#' risks model via \code{crr} from \code{cmprsk} is used instead of the
#' standard Cox proportional hazards model.
#'
#'
#' @return A list containing: \code{D, S, pD, pS} for bootstrapped data and
#'  for permuted data.
#'
#' @importFrom BiocParallel SnowParam MulticoreParam bplapply
#' @importFrom stats model.matrix formula p.adjust
#'


Boot_parallel_survival <- function(BPPARAM = NULL,
                          samples,
                          data,
                          formula.str,
                          meta.info,
                          a1,
                          a2,
                          pSamples,
                          competing_risks) {
    if (is.null(BPPARAM)) {
        if (.Platform$OS.type == "windows") {
            BPPARAM <- SnowParam(workers = 2)
            message("Using SnowParam (Windows) with two workers.")
        } else {
            BPPARAM <- MulticoreParam(workers = 2)
            message("Using MulticoreParam (Unix-like OS) with two workers.")
        }
    } else {
        message("Using provided parallel backend BPPARAM.")
    }
    if (inherits(BPPARAM, "SnowParam")) {
        BPPARAM$exportglobals <- FALSE
    }
    export_vars <- list(
        samples = samples,
        data = data,
        formula.str = formula.str,
        meta.info = meta.info,
        a1 = a1,
        a2 = a2,
        pSamples = pSamples,
        competing_risks = competing_risks
    )
    export_funcs <- list(bootstrap_survival = bootstrap_survival,
                         permutating_survival = permutating_survival)
    for (name in names(export_vars)) {
        assign(name, export_vars[[name]])
    }
    for (name in names(export_funcs)) {
        assign(name, export_funcs[[name]])
    }
    
    results_list <- bplapply(seq_len(nrow(samples)), function(i) {
        samples.R <- split(samples[i, ], rep(1, length(samples[i, ])))
        pSamples_i <- pSamples[[i]]
        d_result <- s_result <- pd_result <- ps_result <- NULL
        if (is.null(a1) | is.null(a2)) {
            fit <- bootstrap_survival(
                x = lapply(samples.R, function(x)
                    data[, x]
                ),
                meta.info = meta.info,
                formula.str = formula.str,
                competing_risks = competing_risks
            )
        }
        d_result <- fit$d
        s_result <- fit$s
        df1 <- data.frame(d_result = d_result, s_result = s_result)
        pFit <- permutating_survival(
            x = data,
            meta.info = pSamples_i,
            formula.str = formula.str,
            competing_risks = competing_risks
        )
        pd_result <- pFit$d
        ps_result <- pFit$s
        df2 <- data.frame(pd_result = pd_result, ps_result = ps_result)
        list(ds = df1, pdps = df2)
    }, BPPARAM = BPPARAM, 
    BPOPTIONS = bpoptions(packages = c("utils", "stringr", "stats", "survival")))
    return(results_list)
}