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
#' @param seed.cl An integer specifying the seed for randomization;
#' if not provided, the default is 1234.
#' @param cluster A parallel cluster object for distributed computation,
#' e.g., created by \code{makeCluster()}. Default is 2.
#' @param formula.str A formula string used when covariates are present in meta.
#' info for modeling. It should include "~ 0 + ..." to exclude the
#' intercept from the model.
#' @param robust indicating whether robust fitting should be used.
#' Default is TRUE, see \link{eBayes}.
#' @param trend indicating whether to include trend fitting in the
#' differential expression analysis. Default is TRUE. see \link{eBayes}.
#' @param permutating.group Logical, If \code{TRUE}, the permutation for
#' calculating the null distribution is performed by permuting the target
#' group only specified in \code{group.name}. If FALSE, the entire
#' \code{meta.info} will be permuted (recommended to be set to FALSE).
#' @param samples bootstrapped samples matrix
#'
#'
#' @return A list containing: \code{D, S, pD, pS} for bootstrapped data and
#'  for permuted data.
#'
#' @importFrom parallel makeCluster clusterSetRNGStream clusterExport
#' stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @import doRNG
#'
#'
#'

Boot_parallel <- function(cluster, seed.cl , samples, data,
                                        formula.str, group.name, groups,
                                        meta.info, a1, a2, trend, robust,
                                        permutating.group) {
    if (is.null(cluster)) {
        cluster <- makeCluster(2)
        registerDoParallel(cluster)
        message("No cluster found; only two cores will be used!")
    } else {
        registerDoParallel(cluster)
    }
    clusterSetRNGStream(cluster, iseed = seed.cl)
    clusterExport( cluster,
            varlist = c("samples", "data",
                            "formula.str", "group.name", "groups", "meta.info",
                                                "a1", "a2", "trend", "robust" ,
                                                "permutating.group"),
                        envir = environment())
    results_list <- foreach(
        i = seq_len(nrow(samples)),
        .combine = "c", .packages = c("utils", "stringr", "stats", "limma"),
        .export = c("Limma_bootstrap",
                        "Limma_permutating")
    ) %dorng% {
        samples.R <- split(samples[i, ], groups)
        d_result <- s_result <- pd_result <- ps_result <- NULL
        if (is.null(a1) | is.null(a2)) {
            fit <- Limma_bootstrap(
                x = lapply(samples.R, function(x) data[, x]),
                    group.name = group.name, meta.info = meta.info,
                formula.str = formula.str, trend = trend, robust = robust
            )
        }
        d_result <- fit$d
        s_result <- fit$s
        df1 <- data.frame(d_result = d_result, s_result = s_result)
        pFit <- Limma_permutating(
            x = lapply(split(seq_len( length(groups) ), groups), function(x)
                    data[, x]), group.name = group.name, meta.info = meta.info,
                        formula.str = formula.str, trend = trend,
            robust = robust, permutating.group = permutating.group )
        pd_result <- pFit$d
        ps_result <- pFit$s
        df2 <- data.frame(pd_result = pd_result, ps_result = ps_result)
        # Return results for this iteration as a data frame
        list(ds = df1, pdps = df2)
    }
    stopCluster(cluster)
    return(results_list)
}
