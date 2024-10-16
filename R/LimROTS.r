#' `LimROTS`: A Hybrid Method Integrating Empirical Bayes and
#' Reproducibility-Optimized Statistics for Robust Analysis of Proteomics and
#' Metabolomics Data
#'
#' @param x A \code{SummarizedExperiment} object or a matrix where rows
#' represent features (e.g., proteins, metabolites) and columns
#' represent samples.
#' The values should be log-transformed.
#' @param B An integer representing the amount of bootstrap iterations.
#' Default is 1000.
#' @param K An optional integer representing the top list size for ranking.
#' If not specified, it is set to one-fourth of the number of features.
#' @param a1 Optional numeric value used in the optimization process.
#' If defined by the user, no optimization occurs.
#' @param a2 Optional numeric value used in the optimization process.
#' If defined by the user, no optimization occurs.
#' @param log Logical, indicating whether the data is already log-transformed.
#' Default is \code{TRUE}.
#' @param progress Logical, indicating whether to display a progress bar the
#' function's execution. Default is \code{FALSE}.
#' @param verbose Logical, indicating whether to display messages during the
#' function's execution. Default is \code{TRUE}.
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
#' @param formula.str A formula string for modeling.
#' It should include "~ 0 + ..." to exclude the intercept from the model.
#' All the model parameters must be present in \code{meta.info}.
#' @param robust indicating whether robust fitting should be used.
#' Default is TRUE, see \link{eBayes}.
#' @param trend indicating whether to include trend fitting in the
#' differential expression analysis. Default is TRUE. see \link{eBayes}.
#' @param permutating.group Logical, If \code{TRUE}, the permutation for
#' calculating the null distribution is performed by permuting the target group
#' only specified in \code{group.name} Preserving all the other sample
#' information. If `FALSE`, the entire sample information retrieved from
#' \code{meta.info} will be permuted (recommended to be set to TRUE).
#'
#'
#' @return An object of class `"list"` with the following elements:
#' \item{data}{The original data matrix.}
#' \item{B}{The number of bootstrap samples used.}
#' \item{d}{The optimized statistics for each feature.}
#' \item{logfc}{Log-fold change values between groups.}
#' \item{pvalue}{P-values computed based on the permutation samples.}
#' \item{FDR}{False discovery rate estimates.}
#' \item{a1}{Optimized parameter used in differential expression ranking.}
#' \item{a2}{Optimized parameter used in differential expression ranking.}
#' \item{k}{Top list size used for ranking.}
#' \item{corrected.logfc}{estimate of the log2-fold-change
#' corresponding to the effect corrected by the s model
#' see \link{topTable}.}
#' \item{q_values}{Estimated q-values using the `qvalue` package.}
#' \item{BH.pvalue}{Benjamini-Hochberg adjusted p-values.}
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
#'     formula.str = formula.str, B = 10
#' )
#'
#' @importFrom stats model.matrix formula p.adjust
#' @importFrom dplyr bind_cols
#' @importFrom qvalue empPvals qvalue
#' @import utils
#'
#' @details The **LimROTS** approach initially uses
#' \link{limma} package functionality to simulate the intensity data of
#' proteins and
#' metabolites. A linear model is subsequently fitted using the design matrix.
#' Empirical Bayes variance shrinking is then implemented. To obtain the
#' moderated t-statistics, the adjusted standard error
#' \eqn{SEpost = √(s2.post)
#' \times unscaled SD} for each feature is computed, along with the regression
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
#'          \link{qvalue} package, along with internal implementation of FDR
#'          adapted from ROTS package. Additionally, the qvalue package is used
#'          to calculate q-values, were the proportion of true null p-values is
#'          set to the bootstrap method \link{pi0est}. We recommend using
#'          permutation-derived p-values and qvalues.
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
                    B = 1000,
                    K = NULL,
                    a1 = NULL,
                    a2 = NULL,
                    log = TRUE,
                    progress = FALSE,
                    verbose = TRUE,
                    meta.info = NULL,
                    cluster = NULL,
                    group.name = NULL,
                    formula.str = NULL,
                    seed.cl = 1234,
                    robust = TRUE,
                    trend = TRUE,
                    permutating.group = FALSE) {

    SanityChecK.list <- SanityChecK(
        x,
        B = B,
        K = K,
        a1 = a1,
        a2 = a2,
        meta.info = meta.info,
        group.name = group.name,
        formula.str = formula.str,
        verbose = verbose,
        log = log
    )
    meta.info <- SanityChecK.list$meta.info
    data <- SanityChecK.list$data
    groups <- SanityChecK.list$groups
    event <- SanityChecK.list$event
    K <- SanityChecK.list$K

    #### FC
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
    if (verbose)
        message("Bootstrapping samples")

    if (ncol(meta.info) > 1) {
        samples <- bootstrapSamples.limRots(B = 2 * B,
                                            meta.info = meta.info,
                                            group.name = group.name)
        pSamples <- NULL
    } else {
        samples <- bootstrapS(2 * B, meta.info, group.name)
    }
    D <- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))
    S <- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))
    pD <- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))
    pS <- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))
    pb <- txtProgressBar(
        min = 0,
        max = 100,
        style = 3
    )
    if (progress) {
        setTxtProgressBar(pb, 50)
    }
    results_list <- Boot_parallel(cluster = cluster, seed.cl = seed.cl,
                                        samples = samples, data = data,
                                        formula.str = formula.str,
                                        group.name = group.name,
                                        groups = groups,
                                        meta.info = meta.info,
                                        a1 = a1, a2 = a2,
                                        trend = trend, robust = robust ,
                                        permutating.group = permutating.group)

    if (progress) {setTxtProgressBar(pb, 80)}

    names(results_list) <- paste0(
        names(results_list),
        seq(1, length(names(results_list)))
    )
    j <- 0
    q <- 0
    for (i in seq_along(results_list)) {
        if (grepl("ds", names(results_list)[i], fixed = TRUE)) {
            j <- j + 1
            D[, j] <- results_list[[names(results_list)[i]]]$d_result
            S[, j] <- results_list[[names(results_list)[i]]]$s_result
        } else {
            q <- q + 1
            pD[, q] <- results_list[[names(results_list)[i]]]$pd_result
            pS[, q] <- results_list[[names(results_list)[i]]]$ps_result
        }
    }
    if (progress) {setTxtProgressBar(pb, 100)}
    if (progress){close(pb)}
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
            Optimizing(B, ssq, N, D, S, pD, pS, verbose, progress)
        a1 <- optimized.parameters$a1
        a2 <- optimized.parameters$a2
        k <- optimized.parameters$k
        R <- optimized.parameters$R
        Z <- optimized.parameters$Z
        ztable <- optimized.parameters$ztable
        fit <- Limma_fit(
            x = lapply(split(seq_len(length(
                groups
            )), groups), function(x)
                data[, x]),
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
        if (verbose)
            message("Calculating p-values")
        p <- empPvals(
            stat = d,
            stat0 = pD,
            pool = TRUE
        )
        if (verbose)
            message("Calculating FDR")
        FDR <- calculateFalseDiscoveryRate(d, pD, progress)
        corrected.logfc <- fit$corrected.logfc
        rm(pD)
        gc()
        q_values <- qvalue(p,
            pi0.method = "bootstrap",
            lambda = seq(0.01, 0.95, 0.01)
        )
        BH.pvalue <- p.adjust(p, method = "BH")
        LimROTS.output <- list(
            data = data,
            B = B,
            d = d,
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
            BH.pvalue = BH.pvalue
        )
    } else {
        fit <- Limma_fit(
            x = lapply(split(seq_len(length(
                groups
            )), groups), function(x)
                data[, x]),
            group.name = group.name,
            meta.info = meta.info,
            formula.str = formula.str,
            trend = trend,
            robust = robust
        )
        d <- fit$d / (a1 + a2 * fit$s)
        pD <- pD / (a1 + a2 * pS)
        if (verbose)
            message("Calculating p-values")
        p <- empPvals(
            stat = d,
            stat0 = pD,
            pool = TRUE
        )
        if (verbose)
            message("Calculating FDR")
        FDR <- calculateFalseDiscoveryRate(d, pD, progress)
        corrected.logfc <- fit$corrected.logfc
        q_values <- qvalue(p,
            pi0.method = "bootstrap",
            lambda = seq(0.01, 0.95, 0.01)
        )
        BH.pvalue <- p.adjust(p, method = "BH")
        LimROTS.output <- list(
            data = data,
            B = B,
            d = d,
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
            BH.pvalue = BH.pvalue
        )
    }
    return(LimROTS.output)
}
