#' `LimROTS`: A Hybrid Method Integrating Empirical Bayes and 
#' Reproducibility-Optimized Statistics for Robust Differential 
#' Expression Analysis
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
#' @param log Logical, indicating whether the data is already log-transformed.
#' Default is \code{TRUE}.
#' @param verbose Logical, indicating whether to display messages during the
#' function's execution. Default is \code{TRUE}.
#' @param meta.info a character vector of the metadata needed for the
#' model to run and can be retrieved using \code{colData()}.
#' @param group.name A string specifying the column in \code{meta.info} that
#' represents the groups or conditions for comparison. group.name should be 
#' retrieved using \code{colData()} as factor.
#' @param BPPARAM   A \code{BiocParallelParam} object specifying the
#' parallelization backend (e.g., \code{MulticoreParam}, \code{SnowParam}).
#' The default depends on the operating system: if the user is on Windows,
#' \code{SnowParam(workers = 2)} is used; otherwise,
#' \code{MulticoreParam(workers = 2)}.
#' @param formula.str A formula string for modeling.
#' It should include "~ 0 + ..." to exclude the intercept from the model.
#' All the model parameters must be present in \code{meta.info}.
#' @param robust indicating whether robust fitting should be used.
#' Default is TRUE, see \link[limma]{eBayes}.
#' @param trend indicating whether to include trend fitting in the
#' differential expression analysis. Default is TRUE. see \link[limma]{eBayes}.
#' @param permutating.group Logical, If \code{TRUE}, the permutation for
#' calculating the null distribution is performed by permuting the target group
#' only specified in \code{group.name} Preserving all the other sample
#' information. If `FALSE`, the entire sample information retrieved from
#' \code{meta.info} will be permuted (recommended to be set to FALSE).
#'
#'
#' @return An object of class `"SummarizedExperiment"` with the 
#' following elements:
#' \item{data}{The original data matrix.}
#' \item{niter}{The number of bootstrap samples used.}
#' \item{statistics}{The optimized statistics for each feature.}
#' \item{logfc}{Log-fold change values between groups.}
#' \item{pvalue}{P-values computed based on the permutation samples.}
#' \item{FDR}{False discovery rate estimates.}
#' \item{a1}{Optimized parameter used in differential expression ranking.}
#' \item{a2}{Optimized parameter used in differential expression ranking.}
#' \item{k}{Top list size used for ranking.}
#' \item{corrected.logfc}{estimate of the log2-fold-change
#' corresponding to the effect corrected by the s model
#' see \link[limma]{topTable}.}
#' \item{q_values}{Estimated q-values using the `qvalue` package.}
#' \item{BH.pvalue}{Benjamini-Hochberg adjusted p-values.}
#' \item{null.statistics}{The optimized null statistics for each feature.}
#'
#'
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
#'     meta.info = meta.info, group.name = "group",
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
#' @details The **LimROTS** approach initially uses
#' \pkg{limma} package functionality to simulate the intensity data of
#' proteins and
#' metabolites. A linear model is subsequently fitted using the design matrix.
#' Empirical Bayes variance shrinking is then implemented. To obtain the
#' moderated t-statistics, the adjusted standard error
#' \eqn{SE_{post} = \sqrt{s^2_{\text{post}} } \times unscaled SD}
#' for each feature is computed, along with the regression
#' coefficient for each feature (indicating the impact of variations in the
#' experimental settings). Then, by adapting a reproducibility-optimized
#' technique known as \link[ROTS]{ROTS} to establish an optimality
#' based on the largest overlap of top-ranked features within group-preserving
#' bootstrap datasets, Finally based on the optimized parameters
#' \eqn{\alpha1} and
#' \eqn{\alpha2} this equation used to calculates the final statistics:
#'
#' \deqn{t_{\alpha_{(p)}} = \frac{\beta_{(p)}}
#'          {\alpha1 + \alpha2 \times SEpost_{(p)}}}where
#'          \eqn{t_{\alpha_{(p)}}} is the final statistics for each feature,
#'          \eqn{\beta_{(p)}} is the coefficient, and \eqn{SEpost_{(p)}}
#'          is the the adjusted
#'          standard error. LimROTS generates p-values from permutation samples
#'          using the implementation available in
#'          \link[qvalue]{qvalue} package, along with internal implementation of FDR
#'          adapted from ROTS package. Additionally, the qvalue package is used
#'          to calculate q-values, were the proportion of true null p-values is
#'          set to the bootstrap method \link[qvalue]{pi0est}. We recommend using
#'          permutation-derived p-values and qvalues.
#'
#' This function processes a dataset using parallel computation. It
#' leverages the \pkg{BiocParallel} framework to distribute tasks
#' across multiple workers, which  can significantly reduce runtime for
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
#'
#'
#' @export


LimROTS <- function(x,
    niter = 1000,
    K = NULL,
    a1 = NULL,
    a2 = NULL,
    log = TRUE,
    verbose = TRUE,
    meta.info,
    BPPARAM  = NULL,
    group.name,
    formula.str,
    robust = TRUE,
    trend = TRUE,
    permutating.group = FALSE) {
    SanityChecK.list <- SanityChecK(
        x,
        niter = niter,
        K = K,
        meta.info = meta.info,
        group.name = group.name,
        verbose = verbose,
        log = log
    )
    meta.info <- SanityChecK.list$meta.info
    data <- SanityChecK.list$data
    groups <- SanityChecK.list$groups
    event <- SanityChecK.list$event
    K <- SanityChecK.list$K
    
    
    if (length(unique(groups)) == 2) {
        group1_data <- data[, groups == 1]
        group2_data <- data[, groups == 2]
        if (log) {
            logfc <-
                rowMeans(group1_data, na.rm = TRUE) -
                rowMeans(group2_data, na.rm = TRUE)
        }
    } else {
        logfc <- rep(NA, nrow(data))
    }
    if (verbose) {
        message("Initiating limma on bootstrapped samples")
    }

    if (ncol(meta.info) > 1) {
        samples <- bootstrapSamples_limRots(
            niter = 2 * niter,
            meta.info = meta.info,
            group.name = group.name
            )
        if(permutating.group == TRUE){
            pSamples <- list()
            for (i in seq_len(nrow(samples)) ) {
                shuffle_df <- meta.info
                shuffle_df[, group.name] <- sample(shuffle_df[, group.name])
                colnames(shuffle_df) <- colnames(meta.info)
                pSamples[[i]] <-  shuffle_df
            }
        }else{
            pSamples <- list()
            for (i in seq_len(nrow(samples)) ) {
                shuffle_df <- meta.info
                shuffle_df <- data.frame(meta.info[sample(nrow(meta.info)), ])
                colnames(shuffle_df) <- colnames(meta.info)
                pSamples[[i]] <-  shuffle_df
            }   
        }
    } else {
        samples <- bootstrapS(2 * niter, meta.info, group.name)
        pSamples <- list()
        for (i in seq_len(nrow(samples)) ) {
            shuffle_df <- meta.info
            shuffle_df <- data.frame(meta.info[sample(nrow(meta.info)), ])
            colnames(shuffle_df) <- colnames(meta.info)
            pSamples[[i]] <-  shuffle_df
        }
        
    }
    D <- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))
    S <- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))
    pD <- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))
    pS <- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))
    results_list <- Boot_parallel(
        BPPARAM  = BPPARAM ,
        samples = samples, data = data,
        formula.str = formula.str,
        group.name = group.name,
        groups = groups,
        meta.info = meta.info,
        a1 = a1, a2 = a2,
        pSamples = pSamples
    )

    for (i in seq_along(results_list)) {
        D[, i] <- results_list[[i]][["ds"]][["d_result"]]
        S[, i] <- results_list[[i]][["ds"]][["s_result"]]
        pD[, i] <- results_list[[i]][["pdps"]][["pd_result"]]
        pS[, i] <- results_list[[i]][["pdps"]][["ps_result"]]
    }
    rm(samples)
    gc()
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
        fit <- Limma_fit(
            x = lapply(split(seq_len(length(
                groups
            )), groups), function(x) {
                data[, x]
            }),
            group.name = group.name,
            meta.info = meta.info,
            formula.str = formula.str,
            trend = trend,
            robust = robust
        )
        d <- fit$d / (a1 + a2 * fit$s)
        pD <- pD / (a1 + a2 * pS)
        rm(pS)
        gc()
        if (verbose) {
            message("Computing p-values and FDR")
        }
        p <- empPvals(
            stat = d,
            stat0 = pD,
            pool = TRUE
        )
        FDR <- calculateFalseDiscoveryRate(d, pD)
        corrected.logfc <- fit$corrected.logfc
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
                logfc = logfc,
                pvalue = p,
                qvalue = q_values$qvalues,
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
            remove(x)
            gc()
        } else {
            LimROTS.output <- list(
                data = data,
                niter = niter,
                statistics = d,
                logfc = logfc,
                pvalue = p,
                FDR = FDR,
                a1 = a1,
                a2 = a2,
                k = k,
                R = R,
                Z = Z,
                ztable = ztable,
                groups = groups,
                corrected.logfc = corrected.logfc,
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
            group.name = group.name,
            meta.info = meta.info,
            formula.str = formula.str,
            trend = trend,
            robust = robust
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
        corrected.logfc <- fit$corrected.logfc
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
                logfc = logfc,
                pvalue = p,
                qvalue = q_values$qvalues,
                FDR = FDR,
                corrected.logfc = corrected.logfc,
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
                logfc = logfc,
                pvalue = p,
                FDR = FDR,
                a1 = a1,
                a2 = a2,
                k = NULL,
                R = NULL,
                Z = NULL,
                groups = groups,
                corrected.logfc = corrected.logfc,
                q_values = q_values,
                BH.pvalue = BH.pvalue,
                null.statistics = pD
            )
        }
    }
    return(LimROTS.output)
}
