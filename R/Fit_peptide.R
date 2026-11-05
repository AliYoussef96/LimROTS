parallel_peptide <- function(BPPARAM = NULL,
                             samples,
                             data,
                             model_formula,
                             group,
                             meta.info,
                             peptide_metadata,
                             pSamples,
                             pairwise_contrasts) {
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
    model_formula = model_formula,
    group = group,
    meta.info = meta.info,
    peptide_metadata = peptide_metadata,
    pSamples = pSamples,
    pairwise_contrasts = pairwise_contrasts
  )
  for (name in names(export_vars)) {
    assign(name, export_vars[[name]])
  }
  
  
  results_list <- bplapply(seq_len(nrow(samples)), function(i) {
    
    samples.R <- split(samples[i, ], meta.info[samples[i, ], group])
    pSamples_i <- pSamples[[i]]
    d_result <- NULL
    s_result <- NULL
    pd_result <- NULL
    ps_result <- NULL
    fit <- fit_bootstrap_peptide(
      x = lapply(samples.R, function(x)
        data[, x]
      ),
      group = group,
      meta.info = meta.info,
      model_formula = model_formula,
      pairwise_contrasts = pairwise_contrasts,
      peptide_metadata = peptide_metadata
    )
    d_result <- fit$d
    s_result <- fit$s
    df1 <- data.frame(d_result = d_result, s_result = s_result)
    pFit <- fit_permutation_peptide(
      x = data,
      group = group,
      meta.info = pSamples_i,
      model_formula = model_formula,
      pairwise_contrasts = pairwise_contrasts,
      peptide_metadata = peptide_metadata
    )
    pd_result <- pFit$d
    ps_result <- pFit$s
    df2 <- data.frame(pd_result = pd_result, ps_result = ps_result)
    list(ds = df1, pdps = df2)
  }, BPPARAM = BPPARAM, 
  BPOPTIONS = bpoptions(packages = c("utils", 
                                     "stringr", "stats", 
                                     "limma" , "variancePartition")))
  return(results_list)
}

fit_original_peptide <- function(x, group, meta.info, model_formula,
                                  pairwise_contrasts, peptide_metadata){
  ### Generate combined data and covariates for model fitting
  combined_data <- x
  covariates.p <- meta.info
  
  stopifnot(identical(colnames(combined_data) , 
                      row.names(covariates.p)))
  
  covariates.p[["sample_id"]] <- row.names(covariates.p)
  
  row.names(covariates.p) <- NULL
  
  stopifnot(identical(colnames(combined_data) , 
                      covariates.p[["sample_id"]]))
  
  stopifnot(identical(row.names(combined_data) , 
                      peptide_metadata[["peptide_id"]]))
  
  
  #============================================================
  # Fast conversion from matrix to long format
  #============================================================
  
  combined_data <- as.data.table(
    combined_data,
    keep.rownames = "peptide_id"
  )
  
  long_data <- melt(
    combined_data,
    id.vars = "peptide_id",
    variable.name = "sample_id",
    value.name = "y",
    variable.factor = FALSE,
    na.rm = TRUE
  )
  
  #============================================================
  # Add peptide metadata
  #============================================================
  
  peptide_metadata <- as.data.table(peptide_metadata)
  
  if (anyDuplicated(peptide_metadata$peptide_id)) {
    stop("peptide_metadata contains duplicated peptide_id values.")
  }
  
  long_data[
    peptide_metadata,
    protein_id := i.protein_id,
    on = "peptide_id"
  ]
  
  #============================================================
  # Add all columns from covariates.p
  #============================================================
  
  covariates.p <- as.data.table(covariates.p)
  
  if (anyDuplicated(covariates.p[["sample_id"]])) {
    stop("covariates.p contains duplicated sample_id values.")
  }
  
  if (!all(long_data[["sample_id"]] %in% covariates.p[["sample_id"]])) {
    stop("Some sample IDs in long_data are missing from covariates.p.")
  }
  
  covariate_columns <- setdiff(
    colnames(covariates.p),
    "sample_id"
  )
  
  duplicated_columns <- intersect(
    covariate_columns,
    colnames(long_data)
  )
  
  if (length(duplicated_columns) > 0) {
    stop(
      paste0(
        "These columns already exist in long_data: ",
        paste(duplicated_columns, collapse = ", ")
      )
    )
  }
  
  long_data[
    covariates.p,
    (covariate_columns) := mget(
      paste0("i.", covariate_columns)
    ),
    on = "sample_id"
  ]
  
  #============================================================
  # Fit one peptide-level mixed model per protein
  #============================================================
  protein_names <- unique(long_data[["protein_id"]])
  
  coefficient <- setNames(
    numeric(length(protein_names)),
    protein_names
  )
  
  standard_error <- setNames(
    numeric(length(protein_names)),
    protein_names
  )
  
  for (current_protein in protein_names) {
    
    protein_data <- long_data[
      protein_id == current_protein
    ]
    
    protein_data[, (group) := factor(get(group))]
    protein_data[, peptide_id := factor(peptide_id)]
    protein_data[, sample_id := factor(sample_id)]
    
    model <- lmer(
      model_formula,
      data = protein_data,
      REML = TRUE
    )
    
    estimated_group_means <- emmeans(
      model,
      specs = group,
      lmer.df = "satterthwaite"
    )
    
    group_levels <- as.character(
      as.data.frame(estimated_group_means)[[group]]
    )
    
    contrast_weights <- limma::makeContrasts(
      contrasts = str_remove_all(pairwise_contrasts, fixed(group)),
      levels = group_levels
    )[, 1]
    
    
    contrast_result <- contrast(
      estimated_group_means,
      method = setNames(
        list(contrast_weights),
        pairwise_contrasts
      ),
      adjust = "none"
    )
    
    contrast_summary <- as.data.frame(
      summary(contrast_result)
    )
    
    contrast_summary <- contrast_summary[contrast_summary[["contrast"]] == pairwise_contrasts,,drop=FALSE]
    
    coefficient[current_protein] <-
      contrast_summary[["estimate"]]
    
    standard_error[current_protein] <-
      contrast_summary[["SE"]]
  }
  
  return(list(d = coefficient, s = standard_error))
                                  }


fit_bootstrap_peptide <- function(x, group, meta.info, model_formula,
                                  pairwise_contrasts, peptide_metadata) {
  
  ### Generate combined data and covariates for model fitting
  data_for_fit <- .generate_data_for_fit(x = x, meta.info = meta.info)
  combined_data <- data_for_fit$combined_data
  covariates.p <- data_for_fit$covariates.p
  
  stopifnot(identical(colnames(combined_data) , 
                      row.names(covariates.p)))
  
  covariates.p[["sample_id"]] <- row.names(covariates.p)
  
  row.names(covariates.p) <- NULL
  
  stopifnot(identical(colnames(combined_data) , 
                      covariates.p[["sample_id"]]))
  
  stopifnot(identical(row.names(combined_data) , 
                      peptide_metadata[["peptide_id"]]))
  
  
  #============================================================
  # Fast conversion from matrix to long format
  #============================================================
  
  combined_data <- as.data.table(
    combined_data,
    keep.rownames = "peptide_id"
  )
  
  long_data <- melt(
    combined_data,
    id.vars = "peptide_id",
    variable.name = "sample_id",
    value.name = "y",
    variable.factor = FALSE,
    na.rm = TRUE
  )
  
  #============================================================
  # Add peptide metadata
  #============================================================
  
  peptide_metadata <- as.data.table(peptide_metadata)
  
  if (anyDuplicated(peptide_metadata$peptide_id)) {
    stop("peptide_metadata contains duplicated peptide_id values.")
  }
  
  long_data[
    peptide_metadata,
    protein_id := i.protein_id,
    on = "peptide_id"
  ]
  
  #============================================================
  # Add all columns from covariates.p
  #============================================================
  
  covariates.p <- as.data.table(covariates.p)
  
  if (anyDuplicated(covariates.p[["sample_id"]])) {
    stop("covariates.p contains duplicated sample_id values.")
  }
  
  if (!all(long_data[["sample_id"]] %in% covariates.p[["sample_id"]])) {
    stop("Some sample IDs in long_data are missing from covariates.p.")
  }
  
  covariate_columns <- setdiff(
    colnames(covariates.p),
    "sample_id"
  )
  
  duplicated_columns <- intersect(
    covariate_columns,
    colnames(long_data)
  )
  
  if (length(duplicated_columns) > 0) {
    stop(
      paste0(
        "These columns already exist in long_data: ",
        paste(duplicated_columns, collapse = ", ")
      )
    )
  }
  
  long_data[
    covariates.p,
    (covariate_columns) := mget(
      paste0("i.", covariate_columns)
    ),
    on = "sample_id"
  ]
  
  #============================================================
  # Fit one peptide-level mixed model per protein
  #============================================================
  protein_names <- unique(long_data[["protein_id"]])
  
  coefficient <- setNames(
    numeric(length(protein_names)),
    protein_names
  )
  
  standard_error <- setNames(
    numeric(length(protein_names)),
    protein_names
  )
  
  for (current_protein in protein_names) {
    
    protein_data <- long_data[
      protein_id == current_protein
    ]
    
    protein_data[, (group) := factor(get(group))]
    protein_data[, peptide_id := factor(peptide_id)]
    protein_data[, sample_id := factor(sample_id)]
    
    model <- lmer(
      model_formula,
      data = protein_data,
      REML = TRUE
    )
    
    estimated_group_means <- emmeans(
      model,
      specs = group,
      lmer.df = "satterthwaite"
    )
    
    group_levels <- as.character(
      as.data.frame(estimated_group_means)[[group]]
    )
    
    contrast_weights <- limma::makeContrasts(
      contrasts = str_remove_all(pairwise_contrasts, fixed(group)),
      levels = group_levels
    )[, 1]
    
    
    contrast_result <- contrast(
      estimated_group_means,
      method = setNames(
        list(contrast_weights),
        pairwise_contrasts
      ),
      adjust = "none"
    )
    
    contrast_summary <- as.data.frame(
      summary(contrast_result)
    )
    
    contrast_summary <- contrast_summary[contrast_summary[["contrast"]] == pairwise_contrasts,,drop=FALSE]
    
    coefficient[current_protein] <-
      abs(contrast_summary[["estimate"]])
    
    standard_error[current_protein] <-
      contrast_summary[["SE"]]
  }
  
  return(list(d = coefficient, s = standard_error))
  
}

fit_permutation_peptide <- function(x, group, meta.info, model_formula,
                                    pairwise_contrasts, peptide_metadata) {
  
  ### Generate combined data and covariates for model fitting
  combined_data <- x
  covariates.p <- meta.info
  
  stopifnot(identical(colnames(combined_data) , 
                      row.names(covariates.p)))
  
  covariates.p[["sample_id"]] <- row.names(covariates.p)
  
  row.names(covariates.p) <- NULL
  
  stopifnot(identical(colnames(combined_data) , 
                      covariates.p[["sample_id"]]))
  
  stopifnot(identical(row.names(combined_data) , 
                      peptide_metadata[["peptide_id"]]))
  
  
  #============================================================
  # Fast conversion from matrix to long format
  #============================================================
  
  combined_data <- as.data.table(
    combined_data,
    keep.rownames = "peptide_id"
  )
  
  long_data <- melt(
    combined_data,
    id.vars = "peptide_id",
    variable.name = "sample_id",
    value.name = "y",
    variable.factor = FALSE,
    na.rm = TRUE
  )
  
  #============================================================
  # Add peptide metadata
  #============================================================
  
  peptide_metadata <- as.data.table(peptide_metadata)
  
  if (anyDuplicated(peptide_metadata$peptide_id)) {
    stop("peptide_metadata contains duplicated peptide_id values.")
  }
  
  long_data[
    peptide_metadata,
    protein_id := i.protein_id,
    on = "peptide_id"
  ]
  
  #============================================================
  # Add all columns from covariates.p
  #============================================================
  
  covariates.p <- as.data.table(covariates.p)
  
  if (anyDuplicated(covariates.p[["sample_id"]])) {
    stop("covariates.p contains duplicated sample_id values.")
  }
  
  if (!all(long_data[["sample_id"]] %in% covariates.p[["sample_id"]])) {
    stop("Some sample IDs in long_data are missing from covariates.p.")
  }
  
  covariate_columns <- setdiff(
    colnames(covariates.p),
    "sample_id"
  )
  
  duplicated_columns <- intersect(
    covariate_columns,
    colnames(long_data)
  )
  
  if (length(duplicated_columns) > 0) {
    stop(
      paste0(
        "These columns already exist in long_data: ",
        paste(duplicated_columns, collapse = ", ")
      )
    )
  }
  
  long_data[
    covariates.p,
    (covariate_columns) := mget(
      paste0("i.", covariate_columns)
    ),
    on = "sample_id"
  ]
  
  #============================================================
  # Fit one peptide-level mixed model per protein
  #============================================================
  protein_names <- unique(long_data[["protein_id"]])
  
  coefficient <- setNames(
    numeric(length(protein_names)),
    protein_names
  )
  
  standard_error <- setNames(
    numeric(length(protein_names)),
    protein_names
  )
  
  for (current_protein in protein_names) {
    
    protein_data <- long_data[
      protein_id == current_protein
    ]
    
    protein_data[, (group) := factor(get(group))]
    protein_data[, peptide_id := factor(peptide_id)]
    protein_data[, sample_id := factor(sample_id)]
    
    model <- lmer(
      model_formula,
      data = protein_data,
      REML = TRUE
    )
    
    estimated_group_means <- emmeans(
      model,
      specs = group,
      lmer.df = "satterthwaite"
    )
    
    group_levels <- as.character(
      as.data.frame(estimated_group_means)[[group]]
    )
    
    contrast_weights <- limma::makeContrasts(
      contrasts = str_remove_all(pairwise_contrasts, fixed(group)),
      levels = group_levels
    )[, 1]
    
    
    contrast_result <- contrast(
      estimated_group_means,
      method = setNames(
        list(contrast_weights),
        pairwise_contrasts
      ),
      adjust = "none"
    )
    
    contrast_summary <- as.data.frame(
      summary(contrast_result)
    )
    
    contrast_summary <- contrast_summary[contrast_summary[["contrast"]] == pairwise_contrasts,,drop=FALSE]
    
    coefficient[current_protein] <-
      abs(contrast_summary[["estimate"]])
    
    standard_error[current_protein] <-
      contrast_summary[["SE"]]
  }
  
  return(list(d = coefficient, s = standard_error))
  
  
}


#' Assemble Bootstrap Data and Matching Metadata
#'
#' Combine group-specific matrices and construct a metadata row for every
#' resampled sample, preserving alignment and separating repeated draws of a
#' correlation block.
#'
#' @param x List of feature-by-sample matrices, one per group in the bootstrap
#' replicate. Matrices must have matching feature names.
#' @param meta.info Data frame of sample-level metadata, with sample names as
#' row names.
#'
#' @return A list containing the combined feature-by-sample matrix
#' `combined_data` and corresponding sample-level metadata `covariates.p`.
.generate_data_for_fit <- function(x, meta.info) {
  ### Check if x is a list of matrices with identical rownames
  stopifnot(
    all(vapply(
      x,
      function(z) identical(rownames(z), rownames(x[[1L]])),
      logical(1)
    ))
  )
  combined_data <- do.call(cbind, x)
  
  ### create a unique colnames for the combined data to avoid duplicates
  colnames(combined_data) <-
    paste0(colnames(combined_data), ".", seq(1, ncol(combined_data)))
  
  ### create a covariates data frame for the combined data
  covariates.p <- data.frame()
  meta.info.temp <- meta.info
  meta.info.temp$sample.id <- row.names(meta.info.temp)
  for (i in colnames(combined_data)) {
    real_SampleNames <- str_split_fixed(i, fixed("."), 2)[, 1]
    df.temp <- meta.info.temp[real_SampleNames, , drop = FALSE]
    df.temp[,"sample.id"] <- i
    covariates.p <- rbind(covariates.p, df.temp)
  }
  stopifnot(identical(colnames(combined_data), covariates.p[,"sample.id"]))
  row.names(covariates.p) <- covariates.p[,"sample.id"]
  covariates.p[,"sample.id"] <- NULL
  
  stopifnot(identical(colnames(combined_data), row.names(covariates.p)))
  return(list(combined_data = combined_data, covariates.p = covariates.p))
}