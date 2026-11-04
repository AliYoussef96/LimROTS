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
#' @param assay.type A character string or numeric index specifying the assay
#' to use if \code{x} is a \code{SummarizedExperiment}. Default is \code{NULL}
#' (uses the first assay).
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
#' a Cox proportional hazards model (or a competing risks model
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
#' set.seed(123)
#' nsamples <- 20
#' nfeatures <- 50
#' sim_data <- matrix(rnorm(nfeatures * nsamples), nrow = nfeatures)
#' colnames(sim_data) <- paste0("sample", seq_len(nsamples))
#' rownames(sim_data) <- paste0("gene", seq_len(nfeatures))
#'
#' meta_data <- data.frame(
#'     time = abs(rnorm(nsamples, mean = 5, sd = 2)),
#'     event = sample(0:1, nsamples, replace = TRUE),
#'     group = factor(rep(seq_len(2), each = nsamples / 2)),
#'     row.names = colnames(sim_data)
#' )
#'
#' formula.str <- "~ Surv(time, event) + group"
#' result <- LimROTS_survival(
#'     x = sim_data,
#'     meta.info = meta_data,
#'     formula.str = formula.str,
#'     niter = 10,
#'     verbose = FALSE,
#'     competing_risks = FALSE
#' )
#'
#'
#' @export


LimROTS_survival <- function(x,
    assay.type = NULL,
    niter = 1000,
    K = NULL,
    a1 = NULL,
    a2 = NULL,
    verbose = TRUE,
    meta.info,
    BPPARAM  = NULL,
    formula.str,
    competing_risks = FALSE) {
        
    ### Perform sanity checks on the input parameters and data.
    SanityChecK.list <- SanityChecK(
        x = x,
        assay.type = assay.type,
        niter = niter,
        K = K,
        meta.info = meta.info,
        group = "event",
        model = NULL,
        pairwise_contrasts = NULL,
        correlation_block = NULL,
        strata = NULL,
        survival = TRUE,
        formula.str = formula.str,
        competing_risks = competing_risks
    )
    meta.info <- SanityChecK.list[["meta.info"]]
    data <- SanityChecK.list[["data"]]
    K <- SanityChecK.list[["K"]]
    formula.str <- SanityChecK.list[["formula.str"]]
    

    ### Generate bootstrap and permutation samples.
    list_boot_perm <- .generate_bootstrap_permutations_survival(
        niter_2 = 2 * niter,
        meta.info = meta.info
    )

    ### Run the bootstrap and permutation fits in parallel.
    results_list <- Boot_parallel_survival(
        BPPARAM = BPPARAM,
        samples = list_boot_perm[["boot_Samples"]],
        data = data,
        formula.str = formula.str,
        meta.info = meta.info,
        a1 = a1,
        a2 = a2,
        pSamples = list_boot_perm[["perm_Samples"]],
        competing_risks = competing_risks
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

    D <- abs(D)
    pD <- abs(pD)

    ### Optimize parameters unless both values were supplied.
    if (is.null(a1) || is.null(a2)) {
        optimized.parameters <- Optimizing(
            niter = niter,
            K = K,
            N = NULL,
            D = D,
            S = S,
            pD = pD,
            pS = pS
        )
        a1 <- optimized.parameters[["a1"]]
        a2 <- optimized.parameters[["a2"]]
        k <- optimized.parameters[["k"]]
        R <- optimized.parameters[["R"]]
        Z <- optimized.parameters[["Z"]]
        ztable <- optimized.parameters[["ztable"]]
    } else {
        k <- R <- Z <- ztable <- NULL
    }

    ### Fit the survival model on the original data.
    fit_survival_result <- .final_fit_survival(
        data = data,
        meta.info = meta.info,
        formula.str = formula.str,
        competing_risks = competing_risks,
        a1 = a1,
        a2 = a2,
        pD = pD,
        pS = pS,
        verbose = verbose
    )

    ### Generate output.
    .generate_output_survival(
        x = x,
        data = data,
        fit = fit_survival_result,
        niter = niter,
        a1 = a1,
        a2 = a2,
        k = k,
        Z = Z,
        R = R,
        ztable = ztable
    )
}

.generate_bootstrap_permutations_survival <- function(niter_2, meta.info) {

    samples <- bootstrapS(niter = niter_2, meta.info = meta.info, group = "event")

    pSamples <- vector("list", niter_2)

    for (i in seq_len(niter_2)) {
        shuffle_df <- meta.info
        idx <- sample(seq_len(nrow(meta.info)))
        shuffle_df[["time"]] <-  meta.info[["time"]][idx]
        shuffle_df[["event"]] <- meta.info[["event"]][idx]
        pSamples[[i]] <- shuffle_df
    }

    list(boot_Samples = samples, perm_Samples = pSamples)
}

.final_fit_survival <- function(data, meta.info, formula.str,
    competing_risks, a1, a2, pD, pS, verbose) {
    fit <- fit_survival(
        x = data,
        meta.info = meta.info,
        formula.str = formula.str,
        competing_risks = competing_risks
    )
    d <- fit[["d"]] / (a1 + a2 * fit[["s"]])
    pD <- pD / (a1 + a2 * pS)

    if (verbose) {
        message("Computing p-values and FDR")
    }
    p_value <- empPvals(stat = d, stat0 = pD, pool = TRUE)
    FDR <- calculateFalseDiscoveryRate(d, pD)
    q_values <- tryCatch(
        qvalue(
            p_value,
            pi0.method = "bootstrap",
            lambda = seq(0.01, 0.95, 0.01)
        ),
        error = function(e) {
            message("qvalue() failed (return NULL): ", e$message)
            NULL
        }
    )
    BH.pvalue <- p.adjust(p_value, method = "BH")

    list(
        d = d,
        pD = pD,
        p_value = p_value,
        FDR = FDR,
        exp_coef = fit[["exp_coef"]],
        q_values = q_values,
        BH.pvalue = BH.pvalue
    )
}

.generate_output_survival <- function(x, data, fit, niter, a1, a2, k, Z,
    R, ztable) {
    if (inherits(x, "SummarizedExperiment")) {
        feat_names <- rownames(x)
        new_rowData <- DataFrame(
            statistics = setNames(fit[["d"]], rownames(data))[feat_names],
            pvalue = setNames(fit[["p_value"]], rownames(data))[feat_names],
            qvalue = if (is.null(fit[["q_values"]])) {
                rep(NA_real_, length(feat_names))
            } else {
                setNames(fit[["q_values"]]$qvalues, rownames(data))[feat_names]
            },
            FDR = setNames(fit[["FDR"]], rownames(data))[feat_names],
            exp_coef = setNames(
                as.numeric(fit[["exp_coef"]]),
                rownames(data)
            )[feat_names],
            BH.pvalue = setNames(fit[["BH.pvalue"]], rownames(data))[feat_names],
            row.names = feat_names
        )
        if (!identical(rownames(new_rowData), feat_names)) {
            stop("Can not add the LimROTS results to the SummarizedExperiment")
        }

        correct.order <- rownames(x)
        if (ncol(rowData(x)) == 0L) {
            rowData(x) <- new_rowData
        } else {
            rowData(x) <- cbind(rowData(x), new_rowData)
        }
        if (!identical(correct.order, rownames(x))) {
            stop("Can not add the LimROTS results to the SummarizedExperiment")
        }
        metadata(x) <- c(metadata(x), list(
            a1 = a1,
            a2 = a2,
            k = k,
            Z = Z,
            R = R,
            ztable = ztable,
            q_values = fit[["q_values"]],
            null.statistics = fit[["pD"]]
        ))
        return(x)
    }

    list(
        data = data,
        niter = niter,
        statistics = fit[["d"]],
        pvalue = fit[["p_value"]],
        FDR = fit[["FDR"]],
        a1 = a1,
        a2 = a2,
        k = k,
        R = R,
        Z = Z,
        ztable = ztable,
        exp_coef = fit[["exp_coef"]],
        q_values = fit[["q_values"]],
        BH.pvalue = fit[["BH.pvalue"]],
        null.statistics = fit[["pD"]]
    )
}