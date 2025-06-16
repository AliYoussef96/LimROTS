#' Parallel processing handling function
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
#' @param group.name A string specifying the column in \code{meta.info} that
#' represents the groups or conditions for comparison.
#' @param BPPARAM A parallel BPPARAM object for distributed computation.
#' @param formula.str A formula string used when covariates are present in meta.
#' info for modeling. It should include "~ 0 + ..." to exclude the
#' intercept from the model.
#' @param samples bootstrapped samples matrix
#' @param groups groups information from `meta.info`
#' @param pSamples a permutated list of samples
#'
#'
#' @return A list containing: \code{D, S, pD, pS} for bootstrapped data and
#'  for permuted data.
#'
#' @importFrom BiocParallel SnowParam MulticoreParam bplapply
#' @importFrom stats model.matrix formula p.adjust
#'


Boot_parallel <- function(BPPARAM = NULL,
    samples,
    data,
    formula.str,
    group.name,
    groups,
    meta.info,
    a1,
    a2,
    pSamples) {
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
        group.name = group.name,
        groups = groups,
        meta.info = meta.info,
        a1 = a1,
        a2 = a2,
        pSamples = pSamples
    )
    export_funcs <- list(Limma_bootstrap = Limma_bootstrap,
                                    Limma_permutating = Limma_permutating)
    for (name in names(export_vars)) {
        assign(name, export_vars[[name]])
    }
    for (name in names(export_funcs)) {
        assign(name, export_funcs[[name]])
    }
    
    results_list <- bplapply(seq_len(nrow(samples)), function(i) {
        samples.R <- split(samples[i, ], groups)
        pSamples_i <- pSamples[[i]]
        d_result <- s_result <- pd_result <- ps_result <- NULL
        if (is.null(a1) | is.null(a2)) {
            fit <- Limma_bootstrap(
                x = lapply(samples.R, function(x)
                    data[, x]
                ),
                group.name = group.name,
                meta.info = meta.info,
                formula.str = formula.str
            )
        }
        d_result <- fit$d
        s_result <- fit$s
        df1 <- data.frame(d_result = d_result, s_result = s_result)
        pFit <- Limma_permutating(
            x = data,
            group.name = group.name,
            meta.info = pSamples_i,
            formula.str = formula.str
            )
        pd_result <- pFit$d
        ps_result <- pFit$s
        df2 <- data.frame(pd_result = pd_result, ps_result = ps_result)
        list(ds = df1, pdps = df2)
    }, BPPARAM = BPPARAM, 
                BPOPTIONS = bpoptions(packages = c("utils", 
                                                            "stringr", "stats", 
                                                                    "limma")))
    return(results_list)
}
