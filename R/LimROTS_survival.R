# Copyright (C) 2025 Ali Youssef.
# This program is licensed under GPL (>= 2).
#
#' `LimROTS_survival`: A Hybrid Method Integrating Empirical Bayes and 
#' Reproducibility-Optimized Statistics for Robust 
#' survival analysis in Omics Data
#'
#' @param x A \code{SummarizedExperiment} object, where rows
#' represent features (e.g., proteins, metabolites) and columns
#' represent samples.
#' The values should be log-transformed.
#' @param niter An integer representing the amount of bootstrap iterations.
#' Default is 1000.
#' @param K An optional integer representing the top list size for ranking.
#' If not specified, it is set to one-fourth of the number of features.
#' @param a1 Optional numeric value used in the optimization process.
#' If defined by the user, no optimization occurs.
#' @param a2 Optional numeric value used in the optimization process.
#' If defined by the user, no optimization occurs.
#' @param verbose Logical, indicating whether to display messages during the
#' function's execution. Default is \code{TRUE}.
#' @param meta.info a character vector of the metadata needed for the
#' model to run and can be retrieved using \code{colData()}.
#' @param BPPARAM   A \code{BiocParallelParam} object specifying the
#' parallelization backend (e.g., \code{MulticoreParam}, \code{SnowParam}).
#' The default depends on the operating system: if the user is on Windows,
#' \code{SnowParam(workers = 2)} is used; otherwise,
#' \code{MulticoreParam(workers = 2)}.
#' @param formula.str A formula string for modeling.
#' It should include "~ 0 + ..." to exclude the intercept from the model.
#' All the model parameters must be present in \code{meta.info}.
#' @param competing_risks Logical. If \code{TRUE}, the function will fit a
#' competing risks model using the \code{crr} function from the
#' \code{cmprsk} package instead of the standard Cox proportional hazards
#' model. Default is \code{FALSE}.
#' @param correlation_block Character or NULL. The name of a column in
#' `meta.info` that defines correlation blocks. Samples sharing the same
#' value in this column are always resampled together as a unit. If NULL,
#' the function behaves identically to \code{bootstrapSamples_limRots}.
#'
#' @return A \code{SummarizedExperiment} when \code{x} is a
#' \code{SummarizedExperiment} (results added to \code{rowData} and
#' \code{metadata}), or a list with the following elements:
#' \item{data}{The original data matrix.}
#' \item{niter}{The number of bootstrap samples used.}
#' \item{statistics}{The optimized statistics for each feature.}
#' \item{pvalue}{P-values computed based on the permutation samples.}
#' \item{FDR}{False discovery rate estimates.}
#' \item{a1}{Optimized parameter used in survival ranking.}
#' \item{a2}{Optimized parameter used in survival ranking.}
#' \item{k}{Top list size used for ranking.}
#' \item{R}{Reproducibility score for the optimized parameters.}
#' \item{Z}{Z-score corresponding to the optimized parameters.}
#' \item{ztable}{Table of z-scores across the parameter grid.}
#' \item{exp_coef}{Exponentiated coefficient (hazard ratio estimate) from
#' the Cox proportional hazards model for each feature.}
#' \item{q_values}{Estimated q-values using the \code{qvalue} package.}
#' \item{BH.pvalue}{Benjamini-Hochberg adjusted p-values.}
#' \item{null.statistics}{The optimized null statistics for each feature.}
#'
#' @importFrom stats model.matrix formula p.adjust
#' @importFrom dplyr bind_cols
#' @importFrom qvalue empPvals qvalue
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom BiocParallel SnowParam MulticoreParam bplapply bpoptions
#' bpRNGseed
#' @importFrom S4Vectors DataFrame metadata
#'
#' @details **LimROTS_survival** applies the reproducibility-optimized
#' statistic framework to survival analysis. For each bootstrap resample,
#' a Cox proportional hazards model (or a Fine\u2013Gray competing risks model
#' via \code{crr} from \pkg{cmprsk} when \code{competing_risks = TRUE})
#' is fitted per feature using the supplied \code{formula.str}. The
#' coefficient \eqn{\beta_{(p)}} and its standard error \eqn{SE_{(p)}} are
#' extracted for each feature across bootstrap datasets.
#'
#' Reproducibility-optimized parameters \eqn{\alpha_1} and \eqn{\alpha_2}
#' are then selected by maximizing the overlap of top-ranked features across
#' group-preserving bootstrap datasets, following the
#' \link[ROTS]{ROTS} approach. The final statistic for each feature is:
#'
#' \deqn{t_{\alpha_{(p)}} = \frac{\beta_{(p)}}
#'          {\alpha_1 + \alpha_2 \times SE_{(p)}}}
#'
#' where \eqn{t_{\alpha_{(p)}}} is the final statistic, \eqn{\beta_{(p)}}
#' is the Cox model coefficient, and \eqn{SE_{(p)}} is its standard error.
#'
#' P-values are computed from permutation-based null distributions using
#' \link[qvalue]{empPvals} from the \pkg{qvalue} package, alongside an
#' internal FDR implementation adapted from the ROTS package. Q-values are
#' estimated via \link[qvalue]{qvalue} with the proportion of true null
#' p-values set using the bootstrap method \link[qvalue]{pi0est}.
#'
#' This function processes a dataset using parallel computation. It
#' leverages the \pkg{BiocParallel} framework to distribute tasks
#' across multiple workers, which can significantly reduce runtime for
#' large datasets.
#'
#' @references
#'   Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W., and Smyth,
#'   G.K. (2015). limma powers differential expression analyses for
#'   RNA-sequencing and microarray studies. Nucleic Acids Research 43(7), e47
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
#' @examples
#' # Simulate a small SummarizedExperiment with survival metadata
#' library(SummarizedExperiment)
#' set.seed(123)
#' nsamples <- 20
#' nfeatures <- 50
#' sim_data <- matrix(rnorm(nfeatures * nsamples), nrow = nfeatures)
#' colnames(sim_data) <- paste0("sample", seq_len(nsamples))
#' rownames(sim_data) <- paste0("gene", seq_len(nfeatures))
#'
#' col_data <- DataFrame(
#'     time = abs(rnorm(nsamples, mean = 5, sd = 2)),
#'     event = sample(0:1, nsamples, replace = TRUE),
#'     group = factor(rep(seq_len(2), each = nsamples / 2))
#' )
#' rownames(col_data) <- colnames(sim_data)
#'
#' se <- SummarizedExperiment(
#'     assays = list(counts = sim_data),
#'     colData = col_data
#' )
#'
#' formula.str <- "~ Surv(time, event) + group"
#' result <- LimROTS_survival(
#'     x = se,
#'     meta.info = c("time", "event", "group"),
#'     formula.str = formula.str,
#'     niter = 10,
#'     verbose = FALSE,
#'     competing_risks = FALSE
#' )
#'
#'
#' @export


LimROTS_survival <- function(x,
    niter = 1000,
    K = NULL,
    a1 = NULL,
    a2 = NULL,
    verbose = TRUE,
    meta.info,
    BPPARAM  = NULL,
    formula.str,
    competing_risks = FALSE,
    correlation_block = NULL) {

    SanityChecK.list <- SanityChecK(
        x,
        niter = niter,
        K = K,
        meta.info = meta.info,
        verbose = verbose,
        log = TRUE,
        survival = TRUE,
        formula.str = formula.str
    )
    meta.info <- SanityChecK.list$meta.info
    data <- SanityChecK.list$data
    K <- SanityChecK.list$K
    formula.str <- SanityChecK.list$formula.str

    samples <- bootstrapSamples_limRots_cox(
        niter = 2 * niter,
        meta.info = meta.info,
        correlation_block = correlation_block)
    

    pSamples <- list()
    for (i in seq_len(nrow(samples)) ) {
        shuffle_df <- meta.info
        shuffle_df <- data.frame(meta.info[sample(nrow(meta.info)), ])
        colnames(shuffle_df) <- colnames(meta.info)
        pSamples[[i]] <-  shuffle_df
    }  
     
    
    D <- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))
    S <- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))
    pD <- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))
    pS <- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))
    
    results_list <- Boot_parallel_survival(
        BPPARAM  = BPPARAM ,
        samples = samples, data = data,
        formula.str = formula.str,
        meta.info = meta.info,
        a1 = a1, a2 = a2,
        pSamples = pSamples,
        competing_risks = competing_risks
    )

    for (i in seq_along(results_list)) {
        D[, i] <- results_list[[i]][["ds"]][["d_result"]]
        S[, i] <- results_list[[i]][["ds"]][["s_result"]]
        pD[, i] <- results_list[[i]][["pdps"]][["pd_result"]]
        pS[, i] <- results_list[[i]][["pdps"]][["ps_result"]]
    }
    
    D <- abs(D)
    pD <- abs(pD)
    
    if (is.null(a1) | is.null(a2)) {
        ssq <- c(seq(0, 20) / 100, seq(11, 50) / 50, seq(6, 25) / 5)
        N <- c(
            seq(1, 20) * 5,
            seq(11, 50) * 10,
            seq(21, 40) * 25,
            seq(11, 1000) * 100
        )
        K <- min(K, nrow(data))
        N <- N[N < K]
        optimized.parameters <-
            Optimizing(niter, ssq, N, D, S, pD, pS, verbose)
        a1 <- optimized.parameters$a1
        a2 <- optimized.parameters$a2
        k <- optimized.parameters$k
        R <- optimized.parameters$R
        Z <- optimized.parameters$Z
        ztable <- optimized.parameters$ztable
        fit <- fit_survival(
            x = data,
            meta.info = meta.info,
            formula.str = formula.str,
            competing_risks = competing_risks
        )
        d <- fit$d / (a1 + a2 * fit$s)
        pD <- pD / (a1 + a2 * pS)

        if (verbose) {
            message("Computing p-values and FDR")
        }
        p <- empPvals(
            stat = d,
            stat0 = pD,
            pool = TRUE
        )
        FDR <- calculateFalseDiscoveryRate(d, pD)
        exp_coef <- fit$exp_coef
        q_values <- tryCatch(
            {
                    qvalue(
                        p,
                        pi0.method = "bootstrap",
                        lambda = seq(0.01, 0.95, 0.01))
            },
            error = function(e) {
                message("qvalue() failed (return NULL): ", e$message)
                NULL
            }
        )
        BH.pvalue <- p.adjust(p, method = "BH")

        if (inherits(x, "SummarizedExperiment")) {
            new_rowData <- DataFrame(
                statistics = d,
                pvalue = p,
                qvalue = q_values$qvalues,
                FDR = FDR,
                exp_coef = exp_coef,
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
            remove(x)
            gc()
        } else {
            LimROTS.output <- list(
                data = data,
                niter = niter,
                statistics = d,
                pvalue = p,
                FDR = FDR,
                a1 = a1,
                a2 = a2,
                k = k,
                R = R,
                Z = Z,
                ztable = ztable,
                groups = groups,
                exp_coef = exp_coef,
                q_values = q_values,
                BH.pvalue = BH.pvalue,
                null.statistics = pD
            )
        }
    } else {
        fit <- Limma_fit(
            x = lapply(split(seq_len(length(
                groups
            )), groups), function(x) {
                data[, x]
            }),
            meta.info = meta.info,
            formula.str = formula.str
        )
        d <- fit$d / (a1 + a2 * fit$s)
        pD <- pD / (a1 + a2 * pS)
        if (verbose) {
            message("Calculating p-values and FDR")
        }
        p <- empPvals(
            stat = d,
            stat0 = pD,
            pool = TRUE
        )
        FDR <- calculateFalseDiscoveryRate(d, pD)
        exp_coef <- fit$exp_coef
        q_values <- tryCatch(
            {
                qvalue(
                    p,
                    pi0.method = "bootstrap",
                    lambda = seq(0.01, 0.95, 0.01)
                )
            },
            error = function(e) {
                message("qvalue() failed (return NULL): ", e$message)
                NULL
            }
        )
        BH.pvalue <- p.adjust(p, method = "BH")

        if (inherits(x, "SummarizedExperiment")) {
            new_rowData <- DataFrame(
                statistics = d,
                pvalue = p,
                qvalue = q_values$qvalues,
                FDR = FDR,
                exp_coef = exp_coef,
                BH.pvalue = BH.pvalue,
                row.names = row.names(data)
            )

            new_rowData <- new_rowData[match(rownames(x), 
                                                        rownames(new_rowData))
                                                            , ]
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
                k = NULL,
                Z = NULL,
                R = NULL,
                ztable = ztable,
                q_values = q_values,
                null.statistics = pD
            ))
            LimROTS.output <- x
            remove(x)
            gc()
        } else {
            LimROTS.output <- list(
                data = data,
                niter = niter,
                statistics = d,
                pvalue = p,
                FDR = FDR,
                a1 = a1,
                a2 = a2,
                k = NULL,
                R = NULL,
                Z = NULL,
                groups = groups,
                exp_coef = exp_coef,
                q_values = q_values,
                BH.pvalue = BH.pvalue,
                null.statistics = pD
            )
        }
    }
    return(LimROTS.output)
}