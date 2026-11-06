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
  SanityChecK.list <- LimROTS:::SanityChecK(
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
  
  ### Extra Perform sanity checks on the peptide_metadata
  
  formula_variables <- as.formula(formula.str)
  formula_variables <- all.vars(formula_variables)
  
  if (formula_variables[1] != group) {
    stop(
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
  
  message("Using the following formula for the model: ", formula.str)
  
  #============================================================
  # Check the intensity matrix
  #============================================================
  
  if (!is.matrix(data)) {
    stop("data must be a matrix.")
  }
  
  if (!is.numeric(data)) {
    stop("data must contain numeric intensities.")
  }
  
  if (is.null(rownames(data))) {
    stop("data must have peptide IDs as row names.")
  }
  
  if (is.null(colnames(data))) {
    stop("data must have sample IDs as column names.")
  }
  
  if (anyDuplicated(rownames(data))) {
    stop("Peptide IDs in data must be unique.")
  }
  
  if (anyDuplicated(colnames(data))) {
    stop("Sample IDs in data must be unique.")
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
    stop("peptide_metadata must contain peptide_id and protein_id.")
  }
  
  if (anyNA(peptide_metadata[["peptide_id"]])) {
    stop("peptide_metadata contains missing peptide IDs.")
  }
  
  if (anyNA(peptide_metadata[["protein_id"]])) {
    stop("peptide_metadata contains missing protein IDs.")
  }
  
  if (anyDuplicated(peptide_metadata[["peptide_id"]])) {
    stop("peptide IDs in peptide_metadata must be unique.")
  }
  
  if (
    !setequal(
      rownames(data),
      peptide_metadata[["peptide_id"]]
    )
  ) {
    stop("Peptide IDs do not align between data and peptide_metadata.")
  }
  
  
  #============================================================
  # Check sample metadata
  #============================================================
  if (anyNA(meta.info)) {
    stop("meta.info contains missing values.")
  }
  
  if (anyDuplicated(meta.info[["sample_id"]])) {
    stop("sample IDs in meta.info must be unique.")
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
    stop("Peptide IDs do not align between data and peptide_metadata after
         reordering.")
  }
  
  
  ### bootstrap and permutation samples
  list_boot_perm <- .generate_bootstrap_permutations_peptide(
    niter_2 = 2 * niter,
    meta.info = meta.info,
    group = group,
    strata = strata
  )
  
  
  ### parallelized bootstrapping and permutation fit
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

  ### optimize parameters
  optimized.parameters <-
        Optimizing(niter =  niter, K = K,
        N = N, D = D, S = S, pD = pD, pS = pS)
    a1 <- optimized.parameters[["a1"]]
    a2 <- optimized.parameters[["a2"]]
    k <- optimized.parameters[["k"]]
    R <- optimized.parameters[["R"]]
    Z <- optimized.parameters[["Z"]]
    ztable <- optimized.parameters[["ztable"]]

    final_peptide <- .fit_final_peptide(
        data = data,
        meta.info = meta.info,
        group = group,
        formula.str = model_formula,
        a1 = a1, a2 = a2,
        pairwise_contrasts = pairwise_contrasts,
        pD = pD, pS = pS, peptide_metadata = peptide_metadata
    )

    d = final_peptide[["d"]]
    #pD = final_peptide[["pD"]]
    p_value = final_peptide[["p_value"]]
    FDR = final_peptide[["FDR"]]
    corrected.logfc = final_peptide[["corrected.logfc"]]
    q_values = final_peptide[["q_values"]]
    BH.pvalue = final_peptide[["BH.pvalue"]]
  
    ### generate output
    LimROTS.output <- .generate_output(
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
    
    return(LimROTS.output)

}


#' Generate Bootstrap and Permutation Samples
#'
#' Generate group-preserving bootstrap samples and corresponding permuted
#' metadata, using optional correlation blocks or strata.
#'
#' @param niter_2 Integer. Number of bootstrap and permutation replicates to
#' generate.
#' @param meta.info Data frame of sample-level metadata, with sample names as
#' row names.
#' @param group Character string naming the column in `meta.info` that
#' identifies each sample's group.
#' @param correlation_block Character string or `NULL`. Name of the column in
#' `meta.info` identifying correlated samples or blocks.
#' @param strata Character vector or `NULL`. Names of character or factor
#' metadata columns used to define strata for resampling and permutation.
#'
#' @return A list containing `boot_Samples`, a matrix of bootstrap sample
#' names, and `perm_Samples`, a list of permuted metadata data frames.
.generate_bootstrap_permutations_peptide <- function(niter_2, meta.info, group, strata) {
  
  
  if(!is.null(strata)){
    message("..........")
    
    samples <- LimROTS:::bootstrapSamples_limRots(
      niter = niter_2,
      meta.info = meta.info,
      group = group, strata = strata
    )
    
    pSamples <- LimROTS:::permutationSamples_limRots(niter = niter_2, meta.info = meta.info, 
                                           group = group, strata = strata)
    
    return(list(boot_Samples = samples, perm_Samples = pSamples))
    
  }else if(is.null(strata)){
    message("..........")
    
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




.fit_final_peptide <- function(data, meta.info, group, formula.str,
        a1, a2, pairwise_contrasts, pD , pS, peptide_metadata) {
        # Fit the limma model on the original data and compute statistics, p-values, and FDR.
        fit <- fit_original_peptide(
            x = data,
            group = group,
            meta.info = meta.info,
            formula.str = formula.str,
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
                message("qvalue() failed (return NULL): ", e$message)
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