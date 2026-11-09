#' Reproducibility-Optimized Peptide-Level Differential Expression Analysis
#'
#' Analyze peptide intensity data using protein-level mixed models,
#' reproducibility-optimized bootstrap statistics, and permutation-based
#' inference.
#'
#' @param x A `SummarizedExperiment` containing peptide features in rows and
#' samples in columns.
#' @param assay.type Character string or numeric index selecting the assay in
#' `x` to analyze. If `NULL`, the first assay is used.
#' @param niter Integer. Number of bootstrap or permutation replicates to
#' generate.
#' @param K Integer. Number of top-ranked features to use. If `NULL`, defaults
#' to one quarter of the number of features.
#' @param meta.info Character vector naming the sample-level metadata columns
#' in `colData(x)` required by the model.
#' @param peptide_metadata Data frame mapping each `peptide_id` to a
#' `protein_id`. Peptide IDs must match the row names of the intensity matrix.
#' @param BPPARAM A `BiocParallelParam` object specifying the parallel backend.
#' If `NULL`, the parallel helper selects a two-worker backend appropriate to
#' the operating system.
#' @param group Character string naming the column in `meta.info` that
#' identifies each sample's group.
#' @param formula.str Character string containing the fixed-effects portion of
#' the model formula. Its first variable must be `group`; random intercepts for
#' peptide and sample are added internally.
#' @param pairwise_contrasts Character vector of contrast expressions, such as
#' `"group1 - group2"`, to fit.
#' @param strata Character vector or `NULL`. Names of character or factor
#' metadata columns used to define strata for resampling and permutation.
#'
#' @return The input `SummarizedExperiment` with protein-level statistics,
#' empirical p-values, q-values, FDR estimates, corrected log-fold changes,
#' and BH-adjusted p-values added to `rowData()`. Optimization parameters and
#' null statistics are added to `metadata()`.
#'
#' @rawNamespace import(data.table, except = shift)
#' @importFrom BiocParallel SnowParam MulticoreParam bplapply bpoptions
#' @importFrom emmeans emmeans contrast
#' @importFrom lme4 lmer
#' @importFrom qvalue empPvals qvalue
#' @importFrom stats as.formula p.adjust setNames
#' @importFrom stringr fixed str_remove_all str_split_fixed
#' 
#' @export 
#' 
LimROTS_peptide <- function(x,  
                            assay.type = NULL,
                            niter = 1000,
                            K = NULL,
                            meta.info,
                            peptide_metadata,
                            BPPARAM  = NULL,
                            group,
                            formula.str,
                            pairwise_contrasts,
                            strata){
  
  ### Perform sanity checks on the input parameters and data.
  cli_alert_info("Validating peptide analysis inputs.")
  SanityChecK.list <- SanityChecK(
    x,
    assay.type = assay.type,
    niter = niter,
    K = K,
    meta.info = meta.info,
    group = group,
    model = "lmer", 
    formula.str = formula.str,
    pairwise_contrasts = pairwise_contrasts,
    correlation_block = NULL,
    strata = strata,
    survival = FALSE,
    competing_risks = FALSE
  )
  
  
  meta.info <- SanityChecK.list[["meta.info"]]
  data <- SanityChecK.list[["data"]]
  group <- SanityChecK.list[["group"]]
  K <- SanityChecK.list[["K"]]
  cli_alert_success("Peptide analysis inputs validated.")
  
  ### Extra Perform sanity checks on the peptide_metadata
  
  formula_variables <- as.formula(formula.str)
  formula_variables <- all.vars(formula_variables)
  
  if (formula_variables[1] != group) {
    cli_abort(
      paste0(
        "The first variable after '~' must be the group variable: ",
        group
      )
    )
  }
  
  formula.str <- paste0(
    "y" ,formula.str,
    " + (1 | peptide_id) + (1 | sample_id)"
  )
  
  model_formula <- as.formula(formula.str)
  
  cli_alert_info("Using the following formula for the model: {formula.str}")
  
  #============================================================
  # Check the intensity matrix
  #============================================================
  
  if (!is.matrix(data)) {
    cli_abort("data must be a matrix.")
  }
  
  if (!is.numeric(data)) {
    cli_abort("data must contain numeric intensities.")
  }
  
  if (is.null(rownames(data))) {
    cli_abort("data must have peptide IDs as row names.")
  }
  
  if (is.null(colnames(data))) {
    cli_abort("data must have sample IDs as column names.")
  }
  
  if (anyDuplicated(rownames(data))) {
    cli_abort("Peptide IDs in data must be unique.")
  }
  
  if (anyDuplicated(colnames(data))) {
    cli_abort("Sample IDs in data must be unique.")
  }
  #============================================================
  # Check peptide metadata
  #============================================================
  
  if (
    !all(
      c("peptide_id", "protein_id") %in%
      colnames(peptide_metadata)
    )
  ) {
    cli_abort("peptide_metadata must contain peptide_id and protein_id.")
  }
  
  if (anyNA(peptide_metadata[["peptide_id"]])) {
    cli_abort("peptide_metadata contains missing peptide IDs.")
  }
  
  if (anyNA(peptide_metadata[["protein_id"]])) {
    cli_abort("peptide_metadata contains missing protein IDs.")
  }
  
  if (anyDuplicated(peptide_metadata[["peptide_id"]])) {
    cli_abort("peptide IDs in peptide_metadata must be unique.")
  }
  
  if (
    !setequal(
      rownames(data),
      peptide_metadata[["peptide_id"]]
    )
  ) {
    cli_abort("Peptide IDs do not align between data and peptide_metadata.")
  }
  
  
  #============================================================
  # Check sample metadata
  #============================================================
  if (anyNA(meta.info)) {
    cli_abort("meta.info contains missing values.")
  }
  
  if (anyDuplicated(meta.info[["sample_id"]])) {
    cli_abort("sample IDs in meta.info must be unique.")
  }
  
  if("sample_id" %in% strata) {
    cli_abort("sample_id should not be a column in meta.info.")
  }

    if("peptide_id" %in% strata) {
    cli_abort("peptide_id should not be a column in meta.info.")
  }

  #============================================================
  # Reorder metadata to exactly match the matrix
  #============================================================
  
  peptide_metadata <- peptide_metadata[
    match(
      rownames(data),
      peptide_metadata[["peptide_id"]]
    ),
    ,
    drop = FALSE
  ]
  
  if(!identical(rownames(data), peptide_metadata[["peptide_id"]])) {
    cli_abort("Peptide IDs do not align between data and peptide_metadata after
         reordering.")
  }
  
  
  ### bootstrap and permutation samples
  cli_alert_info("Generating peptide bootstrap and permutation samples.")
  list_boot_perm <- .generate_bootstrap_permutations_peptide(
    niter_2 = 2 * niter,
    meta.info = meta.info,
    group = group,
    strata = strata
  )
  
  
  ### parallelized bootstrapping and permutation fit
  cli_alert_info("Fitting peptide bootstrap and permutation replicates.")
  results_list <- parallel_peptide(
    BPPARAM  = BPPARAM,
    samples = list_boot_perm[["boot_Samples"]], 
    data = data,
    model_formula = model_formula,
    group = group,
    meta.info = meta.info,
    pSamples = list_boot_perm[["perm_Samples"]],
    pairwise_contrasts = pairwise_contrasts,
    peptide_metadata = peptide_metadata
  )
  cli_alert_success("Peptide bootstrap and permutation fits complete.")

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

  ### optimize parameters
  cli_alert_info("Optimizing peptide reproducibility parameters.")
  optimized.parameters <-
        Optimizing(niter = niter, K = K, D = D, S = S, pD = pD, pS = pS)
    a1 <- optimized.parameters[["a1"]]
    a2 <- optimized.parameters[["a2"]]
    k <- optimized.parameters[["k"]]
    R <- optimized.parameters[["R"]]
    Z <- optimized.parameters[["Z"]]
    ztable <- optimized.parameters[["ztable"]]
    cli_alert_success("Peptide parameter optimization complete.")

    cli_alert_info("Fitting final peptide models and calculating inference.")
    final_peptide <- .fit_final_peptide(
        data = data,
        meta.info = meta.info,
        group = group,
        model_formula = model_formula,
        a1 = a1, a2 = a2,
        pairwise_contrasts = pairwise_contrasts,
        pD = pD, pS = pS, peptide_metadata = peptide_metadata
    )
    cli_alert_success("Final peptide inference complete.")

    d = final_peptide[["d"]]
    #pD = final_peptide[["pD"]]
    p_value = final_peptide[["p_value"]]
    FDR = final_peptide[["FDR"]]
    corrected.logfc = final_peptide[["corrected.logfc"]]
    q_values = final_peptide[["q_values"]]
    BH.pvalue = final_peptide[["BH.pvalue"]]
  
    ### generate output
    LimROTS.output <- .generate_output_peptide(
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
    cli_alert_success("LimROTS peptide analysis complete.")
    
    return(LimROTS.output)

}


#' Generate Peptide Bootstrap and Permutation Samples
#'
#' Generate group-preserving bootstrap samples and corresponding permuted
#' metadata, optionally stratifying samples by the requested metadata columns.
#'
#' @param niter_2 Integer. Number of bootstrap and permutation replicates to
#' generate.
#' @param meta.info Data frame of sample-level metadata, with sample names as
#' row names.
#' @param group Character string naming the column in `meta.info` that
#' identifies each sample's group.
#' @param strata Character vector or `NULL`. Names of character or factor
#' metadata columns used to define strata for resampling and permutation.
#'
#' @return A list containing `boot_Samples`, a matrix of bootstrap sample
#' names, and `perm_Samples`, a list of permuted metadata data frames.
#'
.generate_bootstrap_permutations_peptide <- function(niter_2, meta.info, group, strata) {
  
  
  if(!is.null(strata)){
    cli_alert_info("Generating bootstrap and permutation samples within strata.")
    
    samples <- bootstrapSamples_limRots(
      niter = niter_2,
      meta.info = meta.info,
      group = group, strata = strata
    )
    
    pSamples <- permutationSamples_limRots(niter = niter_2, meta.info = meta.info, 
                                           group = group, strata = strata)
    
    return(list(boot_Samples = samples, perm_Samples = pSamples))
    
  }else if(is.null(strata)){
    cli_alert_info("Generating bootstrap and permutation samples.")
    
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




#' Fit Final Peptide Models and Calculate Inference
#'
#' Fit protein-level mixed models on the original data, scale observed and null
#' statistics, and calculate empirical p-values and false-discovery measures.
#'
#' @param data Feature-by-sample matrix with sample names as column names and
#' peptide IDs as row names.
#' @param meta.info Data frame of sample-level metadata, with sample names as
#' row names.
#' @param group Character string naming the column in `meta.info` that
#' identifies each sample's group.
#' @param formula.str Formula for the peptide-level mixed model.
#' @param a1 Numeric first statistic-scaling parameter.
#' @param a2 Numeric second statistic-scaling parameter.
#' @param pairwise_contrasts Character vector of contrast expressions, such as
#' `"group1 - group2"`, to fit.
#' @param pD Matrix of null-model statistics from the permutation fits.
#' @param pS Matrix of null-model standard errors from the permutation fits.
#' @param peptide_metadata Data frame mapping each `peptide_id` to a
#' `protein_id`.
#'
#' @return A list containing scaled observed statistics `d`, scaled null
#' statistics `pD`, empirical p-values `p_value`, FDR estimates `FDR`,
#' corrected log-fold changes `corrected.logfc`, q-values `q_values`, and
#' BH-adjusted p-values `BH.pvalue`.
#'
.fit_final_peptide <- function(data, meta.info, group, model_formula,
        a1, a2, pairwise_contrasts, pD , pS, peptide_metadata) {
        # Fit the limma model on the original data and compute statistics, p-values, and FDR.
        fit <- fit_original_peptide(
            x = data,
            group = group,
            meta.info = meta.info,
            model_formula = model_formula,
            pairwise_contrasts = pairwise_contrasts,
            peptide_metadata = peptide_metadata
        )
        d <- abs(fit$d) / (a1 + a2 * fit$s)
        pD <- pD / (a1 + a2 * pS)

        p_value <- empPvals(
            stat = d,
            stat0 = pD,
            pool = TRUE
        )
        FDR <- calculateFalseDiscoveryRate(d, pD)
        corrected.logfc <- fit$d
        q_values <- tryCatch(
            {
                    qvalue(
                        p_value,
                        pi0.method = "bootstrap",
                        lambda = seq(0.01, 0.95, 0.01))
            },
            error = function(e) {
                cli_alert_warning("qvalue() failed (return NULL): {e$message}")
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


#' Add LimROTS Peptide Results to a SummarizedExperiment
#'
#' Add peptide-level statistics and associated results to `metadata()` of the
#' input `SummarizedExperiment`.
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
#' @return The `SummarizedExperiment` `x` with peptide-level results and
#' optimization metadata added to `metadata()`.
.generate_output_peptide <- function(x, data, d, p_value, FDR, corrected.logfc,
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
    row.names = names(d)
  )
  

  metadata(x) <- c(metadata(x), list(
    LimROTS_peptide = new_rowData,
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
