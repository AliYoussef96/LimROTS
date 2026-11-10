# Copyright (C) 2025 Ali Youssef.
# This program is licensed under GPL (>= 2).
#
#' Reproducibility-Optimized Survival Analysis
#'
#' Analyze feature-by-sample data using a Cox proportional hazards or
#' competing-risks model and reproducibility-optimized bootstrap statistics.
#'
#' @param x A `SummarizedExperiment` containing features in rows and samples in
#' columns.
#' @param assay.type Character string or numeric index selecting the assay in
#' `x` to analyze. If `NULL`, the first assay is used.
#' @param niter Integer. Number of bootstrap or permutation replicates to
#' generate.
#' @param K Integer. Number of top-ranked features to use. If `NULL`, defaults
#' to one quarter of the number of features.
#' @param meta.info Character vector naming the sample-level metadata columns
#' in `colData(x)` required for the survival model. Must include `time`,
#' `event`, and the formula covariates.
#' @param BPPARAM A `BiocParallelParam` object specifying the parallel backend.
#' If `NULL`, the parallel helper selects a two-worker backend appropriate to
#' the operating system.
#' @param formula.str Character string containing the survival model formula.
#' It must include `Surv(time,event)~y`; additional covariates must be columns
#' in `meta.info`.
#' @param competing_risks Logical. If `TRUE`, fit a Fine-Gray competing-risks
#' model with `cmprsk::crr()`; otherwise fit a Cox proportional hazards model
#' with `survival::coxph()`.
#'
#' @return The input `SummarizedExperiment` with feature-level statistics,
#' empirical p-values, q-values, FDR estimates, exponentiated coefficients,
#' and BH-adjusted p-values added to `rowData()`. Optimization results and
#' null statistics are added to `metadata()`.
#'
#' @importFrom stats model.matrix formula p.adjust setNames
#' @importFrom dplyr bind_cols
#' @importFrom qvalue empPvals qvalue
#' @importFrom utils txtProgressBar setTxtProgressBar combn
#' @importFrom BiocParallel SnowParam MulticoreParam bplapply bpoptions
#' bpRNGseed
#' @import survival
#' @importFrom cmprsk crr
#' @importFrom stringr str_split_fixed fixed str_detect
#' @importFrom S4Vectors DataFrame metadata
#'
#' @details The function applies the reproducibility-optimized statistic
#' framework to survival analysis. For each bootstrap resample, a Cox
#' proportional hazards model or Fine-Gray competing-risks model is fitted
#' per feature using `formula.str`. The feature coefficient
#' \eqn{\beta_{(p)}} and its standard error \eqn{SE_{(p)}} are extracted across
#' bootstrap datasets.
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
#' is the feature coefficient, and \eqn{SE_{(p)}} is its standard error.
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
#'     competing_risks = FALSE
#' )
#'
#' @export 
LimROTS_survival <- function(x,
    assay.type = NULL,
    niter = 1000,
    K = NULL,
    meta.info,
    BPPARAM  = NULL,
    formula.str,
    competing_risks = FALSE) {
    
    cli_alert_info("Validating survival-analysis inputs.")
    formula.str <- gsub("\\s+", "",formula.str)
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
    cli_alert_info("Generating survival bootstrap and permutation samples.")
    list_boot_perm <- .generate_bootstrap_permutations_survival(
        niter_2 = 2 * niter,
        meta.info = meta.info
    )

    ### Run the bootstrap and permutation fits in parallel.
    cli_alert_info("Fitting survival bootstrap and permutation replicates.")
    results_list <- Boot_parallel_survival(
        BPPARAM = BPPARAM,
        samples = list_boot_perm[["boot_Samples"]],
        data = data,
        formula.str = formula.str,
        meta.info = meta.info,
        pSamples = list_boot_perm[["perm_Samples"]],
        competing_risks = competing_risks
    )
    cli_alert_success("Survival bootstrap and permutation fits complete.")

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
        cli_abort("D, S, pD, and pS must not contain NA values.")
    }
    
    if (!all(vapply(list(S, pD, pS), function(x) identical(dim(x), dim(D)), logical(1)))) {
        cli_abort("D, S, pD, and pS must all have identical dimensions.")
    }

    if (any(vapply(list(D, S, pD, pS), function(x) any(!is.finite(x)), logical(1)))) {
    cli_abort("D, S, pD, and pS must not contain NA, NaN, Inf, or -Inf.")
    }

    if (any(vapply(list(D, S, pD, pS), function(x) any(x < 0), logical(1)))) {
    cli_abort("D, S, pD, and pS must contain only non-negative values.")
    }

    ### Optimize parameters unless both values were supplied.
    cli_alert_info("Optimizing survival reproducibility parameters.")
    optimized.parameters <-
        Optimizing(niter = niter, K = K, D = D, S = S, pD = pD, pS = pS)
    a1 <- optimized.parameters[["a1"]]
    a2 <- optimized.parameters[["a2"]]
    k <- optimized.parameters[["k"]]
    R <- optimized.parameters[["R"]]
    Z <- optimized.parameters[["Z"]]
    ztable <- optimized.parameters[["ztable"]]
    cli_alert_success("Survival parameter optimization complete.")


    ### Fit the survival model on the original data.
    cli_alert_info("Fitting the final survival model and calculating inference.")
    fit_survival_result <- .final_fit_survival(
        data = data,
        meta.info = meta.info,
        formula.str = formula.str,
        competing_risks = competing_risks,
        a1 = a1,
        a2 = a2,
        pD = pD,
        pS = pS
    )
    cli_alert_success("Final survival inference complete.")

    d = fit_survival_result[["d"]]
    #pD = fit_limma[["pD"]]
    p_value = fit_survival_result[["p_value"]]
    FDR = fit_survival_result[["FDR"]]
    q_values = fit_survival_result[["q_values"]]
    BH.pvalue = fit_survival_result[["BH.pvalue"]]

    ### generate output

    LimROTS.output <- .generate_output(
        x = x,
        data = data,
        d = d,
        p_value = p_value,
        FDR = FDR,
        corrected.logfc = rep(NA, length(d)),
        q_values = q_values,
        BH.pvalue = BH.pvalue,
        a1 = a1, a2 = a2, k = k, Z = Z, R = R, ztable = ztable, pD = pD
    )
    cli_alert_success("LimROTS survival analysis complete.")
}

#' Generate Survival Bootstrap and Permutation Samples
#'
#' Generate bootstrap sample names grouped by event status and permutation
#' replicates by jointly shuffling survival times and event labels.
#'
#' @param niter_2 Integer. Number of bootstrap and permutation replicates to
#' generate.
#' @param meta.info Data frame of sample-level metadata, with sample names as
#' row names and an `event` column.
#'
#' @return A list containing `boot_Samples`, a matrix of bootstrap sample names,
#' and `perm_Samples`, a list of permuted metadata data frames.
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

#' Fit Survival Models and Calculate Inference
#'
#' Fit per-feature survival models on the original data, scale observed and
#' null statistics with the supplied parameters, and calculate empirical
#' p-values and false-discovery measures.
#'
#' @param data Feature-by-sample matrix with sample names as column names.
#' @param meta.info Data frame of sample-level metadata, with sample names as
#' row names.
#' @param formula.str Character string containing the survival model formula.
#' @param competing_risks Logical. If `TRUE`, use a Fine-Gray competing-risks
#' model; otherwise use a Cox proportional hazards model.
#' @param a1 Numeric first statistic-scaling parameter.
#' @param a2 Numeric second statistic-scaling parameter.
#' @param pD Matrix of null-model coefficients from the permutation fits.
#' @param pS Matrix of null-model standard errors from the permutation fits.
#'
#' @return A list containing observed statistics `d`, scaled null statistics
#' `pD`, empirical p-values `p_value`, FDR estimates `FDR`, exponentiated
#' coefficients `exp_coef`, q-values `q_values`, and BH-adjusted p-values
#' `BH.pvalue`.
.final_fit_survival <- function(data, meta.info, formula.str,
    competing_risks, a1, a2, pD, pS) {
    fit <- fit_survival(
        x = data,
        meta.info = meta.info,
        formula.str = formula.str,
        competing_risks = competing_risks
    )
    d <- fit[["d"]] / (a1 + a2 * fit[["s"]])
    pD <- pD / (a1 + a2 * pS)

    p_value <- empPvals(
        stat = d, 
        stat0 = pD, 
        pool = TRUE
    )
    FDR <- calculateFalseDiscoveryRate(d, pD)
    q_values <- tryCatch(
        qvalue(
            p_value,
            pi0.method = "bootstrap",
            lambda = seq(0.01, 0.95, 0.01)
        ),
        error = function(e) {
            cli_alert_warning("qvalue() failed (return NULL): {e$message}")
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