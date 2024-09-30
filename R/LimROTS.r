#' LimROTS: An Extension of the ROTS Method with Limma Integration
#'
#' The `LimROTS` function employs a reproducibility-optimized test statistic utilising the limma and ROTS methodology to simulate complex experimental designs.
#'
#' @param x A \code{SummarizedExperiment} object or a matrix where rows represent features (e.g., genes, proteins) and columns represent samples. The values should be log-transformed.
#' @param B An integer specifying the number of bootstrap iterations. Default is 1000.
#' @param K An optional integer representing the top list size for ranking. If not specified, it is set to one-fourth of the number of features.
#' @param a1 Optional numeric value used in the optimization process. If defined by the user, no optimization occurs.
#' @param a2 Optional numeric value used in the optimization process. If defined by the user, no optimization occurs.
#' @param log Logical, indicating whether the data is already log-transformed. Default is \code{TRUE}.
#' @param progress Logical, indicating whether to display a progress bar during bootstrap sampling. Default is \code{FALSE}.
#' @param verbose Logical, indicating whether to display messages during the function's execution. Default is \code{TRUE}.
#' @param meta.info A data frame containing sample-level metadata, where each row corresponds to a sample. It should include the grouping variable specified in \code{group.name}. If \code{x} is a \code{SummarizedExperiment} object, \code{meta.info} must be a vector of the metadata needed for the model to run and can be retrieved using \code{colData()}.
#' @param group.name A string specifying the column in \code{meta.info} that represents the groups or conditions for comparison.
#' @param seed.cl An integer specifying the seed for randomization; if not provided, the default is 1234.
#' @param cluster A parallel cluster object for distributed computation, e.g., created by \code{makeCluster()}. Default is 2.
#' @param survival Logical, indicating whether to enable survival analysis. If \code{TRUE}, then \code{meta.info} should contain \code{time} and \code{event} columns.
#' @param paired Logical, indicating whether the data represent paired samples. Default is \code{FALSE}.
#' @param n.ROTS Logical. If \code{TRUE}, all parameters related to \code{LimROTS} will be ignored, and the original \code{ROTS} analysis will run. This must be \code{TRUE} when \code{survival} or \code{paired} is set to \code{TRUE}.
#' @param formula.str A formula string used when covariates are present in meta.info for modeling. It should include "~ 0 + ..." to exclude the intercept from the model.
#' @param robust indicating whether robust fitting should be used. Default is TRUE, see \link{eBayes}.
#' @param trend indicating whether to include trend fitting in the differential expression analysis. Default is TRUE. see \link{eBayes}.
#'
#'
#' @return A list of class `"list"` with the following elements:
#' \item{data}{The original data matrix.}
#' \item{B}{The number of bootstrap samples used.}
#' \item{d}{Differential expression statistics for each feature.}
#' \item{logfc}{Log-fold change values between groups.}
#' \item{pvalue}{P-values computed based on the permutation samples.}
#' \item{FDR}{False discovery rate estimates.}
#' \item{a1}{Optimized parameter used in differential expression ranking.}
#' \item{a2}{Optimized parameter used in differential expression ranking.}
#' \item{k}{Top list size used for ranking.}
#' \item{corrected.logfc}{Corrected log-fold change values, if applicable.}
#' \item{q_values}{Estimated q-values using the `qvalue` package.}
#' \item{BH.pvalue}{Benjamini-Hochberg adjusted p-values.}
#'
#'
#'
#' @examples
#' # Example usage:
#'
#' data <- data.frame(matrix(rnorm(500), nrow = 100, ncol = 10)) # Simulated data
#' meta.info <- data.frame(group = factor(rep(1:2, each = 5)), row.names = colnames(data))
#' formula.str <- "~ 0 + group"
#' result <- LimROTS(data, meta.info = meta.info, group.name = "group",
#'                                 formula.str = formula.str, B = 10)
#'
#' @importFrom limma voom lmFit eBayes
#' @importFrom stats model.matrix formula p.adjust
#' @importFrom dplyr bind_cols
#' @importFrom parallel makeCluster clusterSetRNGStream clusterExport stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @import doRNG
#' @importFrom qvalue empPvals qvalue
#' @import utils
#' @import SummarizedExperiment
#'
#' @details Differential expression analysis is a prevalent method utilised in the examination of diverse biological data.
#'     The reproducibility-optimized test statistic (ROTS) modifies a t-statistic based on the data's intrinsic characteristics and ranks features according to their statistical significance for differential expression between two or more groups, as shown by the f-statistic.
#'     Focussing on proteomics and metabolomics, the current ROTS implementation cannot account for technical or biological covariates such as MS batches or gender differences among the samples.
#'     Consequently, we developed LimROTS, which employs a reproducibility-optimized test statistic utilising the limma methodology to simulate complex experimental designs.
#'
#'
#' @references
#'   Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W., and Smyth, G.K. (2015). limma powers differential
#'   expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research 43(7), e47
#'
#'   Suomi T, Seyednasrollah F, Jaakkola M, Faux T, Elo L (2017). “ROTS: An R package for reproducibility-optimized
#'   statistical testing.” _PLoS computational biology_, *13*(5), e1005562. \url{doi:10.1371/journal.pcbi.1005562}
#'   \url{https://doi.org/10.1371/journal.pcbi.1005562}, \url{http://www.ncbi.nlm.nih.gov/pubmed/28542205}
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
                    survival = FALSE,
                    paired = FALSE,
                    n.ROTS = FALSE,
                    seed.cl = 1234,
                    robust = TRUE,
                    trend = TRUE)
{
    SanityChecK.list <- SanityChecK(
        x,
        B = B,
        K = K,
        a1 = a1,
        a2 = a2,
        meta.info = meta.info,
        group.name = group.name,
        formula.str = formula.str,
        survival = survival,
        paired = paired,
        n.ROTS = n.ROTS,
        verbose = verbose,
        log = log
    )

    meta.info <- SanityChecK.list$meta.info
    data <- SanityChecK.list$data
    groups <- SanityChecK.list$groups
    event <- SanityChecK.list$event
    K <- SanityChecK.list$K

    #### FC

    if (length(unique(groups)) == 2 && survival == FALSE) {
        group1_data <- data[, groups == 1]
        group2_data <- data[, groups == 2]

        if (log) {
            logfc <- rowMeans(group1_data, na.rm = TRUE) - rowMeans(group2_data, na.rm = TRUE)
        } else {
            logfc <- rowMeans(log2(group1_data + 1), na.rm = TRUE) - rowMeans(log2(group2_data + 1), na.rm = TRUE)
        }
    } else {
        logfc <- rep(NA, nrow(data))
    }


    if (verbose)
        message("Bootstrapping samples")


    if (n.ROTS == FALSE) {
        if (ncol(meta.info) > 1) {
            samples <- bootstrapSamples.limRots(
                B = 2 * B,
                meta.info = meta.info,
                group.name = group.name
            )
            pSamples <- NULL
        } else {
            paired <- FALSE
            samples <- bootstrapS(2 * B, meta.info, group.name, paired)
            pSamples <- NULL
        }
    } else {
        samples <- bootstrapS(2 * B, meta.info, group.name, paired)
        pSamples <- permutatedS(meta.info, 2 * B)
    }

    D <- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))
    S <- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))
    pD <- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))
    pS <- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))

    pb <- txtProgressBar(min = 0,
        max = 100,
        style = 3)

    if (is.null(cluster)) {
        cluster <- makeCluster(2)
        registerDoParallel(cluster)
        message("No cluster found; only two cores will be used!")
    } else {
        registerDoParallel(cluster)
    }

    clusterSetRNGStream(cluster, iseed = seed.cl)
    clusterExport(
        cluster,
        varlist = c(
            "pb",
            "samples",
            "pSamples",
            "D",
            "data",
            "S",
            "pD",
            "pS",
            "time",
            "formula.str",
            "group.name",
            "groups",
            "event",
            "meta.info",
            "a1",
            "a2",
            "trend",
            "robust",
            "n.ROTS",
            "survival"
        ),
        envir = environment()
    )



    if (progress) {
        setTxtProgressBar(pb, 50)
    }
    results_list <- foreach(
        i = seq_len(nrow(samples)),
        .combine = "c",
        .packages = c("utils", "stringr", "stats" , "limma"),
        .export = c("testStatSurvivalOptimized" , "testStatistic_with_covariates" , "testStatOptimized", "testStatistic_with_covariates_permutating")
    ) %dorng% {
        samples.R <- split(samples[i, ], groups)

        # Initialize placeholders for results
        d_result <- s_result <- pd_result <- ps_result <- NULL

        # Compute D and S if conditions are met
        if (is.null(a1) | is.null(a2)) {
            if (survival == TRUE) {
                fit <- testStatSurvivalOptimized(lapply(samples.R, function(x)
                    data[, x]),
                groups,
                event)

            } else if (n.ROTS == FALSE) {
                fit <- testStatistic_with_covariates(
                    x = lapply(samples.R, function(x)
                        data[, x]),
                    group.name = group.name,
                    meta.info = meta.info,
                    formula.str = formula.str,
                    trend =
                        trend,
                    robust = robust
                )

            } else {
                fit <- testStatOptimized(paired, lapply(samples.R, function(x)
                    data[, x]))
            }
            d_result <- fit$d
            s_result <- fit$s

            df1 <- data.frame(d_result = d_result, s_result = s_result)
        }

        # Compute pD and pS
        if (survival == TRUE) {
            pSamples.R <- split(pSamples[i, ], groups)

            pFit <- testStatSurvivalOptimized(lapply(pSamples.R, function(x)
                data[, x]), groups, event)
        } else if (n.ROTS == FALSE) {
            pFit <- testStatistic_with_covariates_permutating(
                x = lapply(split(seq_len(
                    length(groups)
                ), groups), function(x)
                    data[, x]),
                group.name = group.name,
                meta.info = meta.info,
                formula.str = formula.str,
                trend =
                    trend,
                robust = robust
            )
        } else {
            pSamples.R <- split(pSamples[i, ], groups)

            pFit <- testStatOptimized(paired, lapply(pSamples.R, function(x)
                data[, x]))
        }
        pd_result <- pFit$d
        ps_result <- pFit$s

        df2 <- data.frame(pd_result = pd_result, ps_result = ps_result)


        # Return results for this iteration as a data frame
        list(ds = df1, pdps = df2)
    }

    stopCluster(cluster)

    if (progress) {
        setTxtProgressBar(pb, 80)
    }

    names(results_list) <- paste0(names(results_list), seq(1, length(names(results_list))))

    j <- 0
    q <- 0
    # Populate matrices D, S, pD, pS from results
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


    if (progress) {
        setTxtProgressBar(pb, 100)
    }


    if (progress)
        close(pb)
    rm(samples, pSamples)
    gc()

    if (is.null(a1) | is.null(a2)) {
        ssq <- c(seq(0, 20) / 100, seq(11, 50) / 50, seq(6, 25) / 5)
        N <- c(seq(1, 20) * 5,
            seq(11, 50) * 10,
            seq(21, 40) * 25,
            seq(11, 1000) *
                100)
        K <- min(K, nrow(data))
        N <- N[N < K]

        optimized.parameters <- Optimizing(B, ssq, N, D, S, pD, pS, verbose, progress)


        a1 <- optimized.parameters$a1
        a2 <- optimized.parameters$a2
        k <- optimized.parameters$k
        R <- optimized.parameters$R
        Z <- optimized.parameters$Z
        ztable <- optimized.parameters$ztable

        if (survival == TRUE) {
            fit <- testStatSurvivalOptimized(lapply(split(seq_len(
                length(groups)
            ), groups), function(x)
                data[, x]), groups, event)
        } else if (n.ROTS == FALSE) {
            fit <- testStatistic_with_covariates_Fit(
                x = lapply(split(seq_len(
                    length(groups)
                ), groups), function(x)
                    data[, x]),
                group.name = group.name,
                meta.info = meta.info,
                formula.str = formula.str,
                trend = trend,
                robust = robust
            )
        } else {
            fit <- testStatOptimized(paired, lapply(split(seq_len(
                length(groups)
            ), groups), function(x)
                data[, x]))
        }
        d <- fit$d / (a1 + a2 * fit$s)
        pD <- pD / (a1 + a2 * pS)
        rm(pS)
        gc()
        if (verbose)
            message("Calculating p-values")
        p <- empPvals(stat = d,
            stat0 = pD,
            pool = TRUE)
        if (verbose)
            message("Calculating FDR")

        FDR <- calculateFalseDiscoveryRate(d, pD, progress)
        corrected.logfc <- fit$corrected.logfc
        rm(pD)
        gc()

        q_values <- qvalue(p,
            pi0.method = "bootstrap",
            lambda = seq(0.01, 0.95, 0.01))
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
    }
    else {
        if (survival == TRUE) {
            fit <- testStatSurvivalOptimized(lapply(split(seq_len(
                length(groups)
            ), groups), function(x)
                data[, x]), groups, event)
        } else if (n.ROTS == FALSE) {
            fit <- testStatistic_with_covariates_Fit(
                x = lapply(split(seq_len(
                    length(groups)
                ), groups), function(x)
                    data[, x]),
                group.name = group.name,
                meta.info = meta.info,
                formula.str = formula.str,
                trend = trend,
                robust = robust
            )
        } else {
            fit <- testStatOptimized(paired, lapply(split(seq_len(
                length(groups)
            ), groups), function(x)
                data[, x]))
        }

        d <- fit$d / (a1 + a2 * fit$s)
        pD <- pD / (a1 + a2 * pS)
        if (verbose)
            message("Calculating p-values")
        p <- empPvals(stat = d,
            stat0 = pD,
            pool = TRUE)
        if (verbose)
            message("Calculating FDR")
        FDR <- calculateFalseDiscoveryRate(d, pD, progress)
        corrected.logfc <- fit$corrected.logfc
        q_values <- qvalue(p,
            pi0.method = "bootstrap",
            lambda = seq(0.01, 0.95, 0.01))
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
