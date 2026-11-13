#' Fit Peptide Bootstrap and Permutation Replicates in Parallel
#'
#' Fit the peptide-level mixed model for each bootstrap replicate and its
#' corresponding permuted metadata using a `BiocParallel` backend.
#'
#' @param BPPARAM A `BiocParallelParam` object specifying the parallel backend.
#' If `NULL`, the function selects a two-worker backend appropriate to the
#' operating system.
#' @param samples Matrix of bootstrap replicates, with one row per replicate
#' and resampled sample names in its columns.
#' @param data Feature-by-sample matrix with sample names as column names.
#' @param model_formula Formula for the peptide-level mixed model.
#' @param group Character string naming the column in `meta.info` that
#' identifies each sample's group.
#' @param meta.info Data frame of sample-level metadata, with sample names as
#' row names.
#' @param peptide_metadata Data frame mapping each `peptide_id` to a
#' `protein_id`.
#' @param pSamples List of permuted sample-level metadata data frames, one per
#' bootstrap replicate.
#' @param pairwise_contrasts Character vector of contrast expressions, such as
#' `"group1 - group2"`, to fit.
#'
#' @return A list with one element per replicate. Each element contains `ds`
#' (bootstrap statistics and standard errors) and `pdps` (permutation
#' statistics and standard errors).
#'
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
      cli_alert_info("Using SnowParam (Windows) with two workers.")
    } else {
      BPPARAM <- MulticoreParam(workers = 2)
      cli_alert_info("Using MulticoreParam (Unix-like OS) with two workers.")
    }
  } else {
    cli_alert_info("Using the provided parallel backend.")
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
                                     "data.table", "lme4", "emmeans")))
  return(results_list)
}

#' Fit the Original Peptide Data
#'
#' Fit a peptide-level linear mixed model for each protein using the original
#' expression matrix, sample metadata, and peptide-to-protein mapping.
#'
#' @param x Feature-by-sample matrix with sample names as column names and
#' peptide IDs as row names.
#' @param group Character string naming the column in `meta.info` that
#' identifies each sample's group.
#' @param meta.info Data frame of sample-level metadata, with sample names as
#' row names.
#' @param model_formula Formula for the peptide-level mixed model.
#' @param pairwise_contrasts Character vector of contrast expressions, such as
#' `"group1 - group2"`, to fit.
#' @param peptide_metadata Data frame mapping each `peptide_id` to a
#' `protein_id`.
#'
#' @return A list containing signed protein-level contrast estimates `d` and
#' their standard errors `s`.
#'
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
    cli_abort("peptide_metadata contains duplicated peptide_id values.")
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
    cli_abort("covariates.p contains duplicated sample_id values.")
  }
  
  if (!all(long_data[["sample_id"]] %in% covariates.p[["sample_id"]])) {
    cli_abort("Some sample IDs in long_data are missing from covariates.p.")
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
    cli_abort(
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
  #============================================================
  # Prepare columns only once
  #============================================================

long_data[, (group) := factor(get(group))]
long_data[, peptide_id := factor(peptide_id)]
long_data[, sample_id := factor(sample_id)]

# Avoid repeatedly scanning the complete long_data table
data.table::setkey(long_data, protein_id)

protein_names <- unique(long_data[["protein_id"]])

# Cache contrasts because most proteins have the same group levels
contrast_cache <- new.env(parent = emptyenv())

get_contrast_weights <- function(group_levels) {

  cache_key <- paste(group_levels, collapse = "\r")

  if (!exists(cache_key, envir = contrast_cache, inherits = FALSE)) {

    weights <- limma::makeContrasts(
      contrasts = stringr::str_remove_all(
        pairwise_contrasts,
        stringr::fixed(group)
      ),
      levels = group_levels
    )[, 1]

    assign(cache_key, weights, envir = contrast_cache)
  }

  get(cache_key, envir = contrast_cache, inherits = FALSE)
}

fit_one_protein <- function(current_protein) {

  protein_data <- long_data[.(current_protein)]

  # Make failures explicit rather than silently corrupting the output
  result <- tryCatch({

    present_groups <- unique(as.character(protein_data[[group]]))

    if (length(present_groups) < 2L) {
      return(c(d = NA_real_, s = NA_real_))
    }

    # droplevels is important if some proteins are missing a group
    protein_data[, (group) := droplevels(get(group))]
    protein_data[, peptide_id := droplevels(peptide_id)]
    protein_data[, sample_id := droplevels(sample_id)]

    model <- lme4::lmer(
      formula = model_formula,
      data = protein_data,
      REML = TRUE,
      control = lme4::lmerControl(
        optimizer = "nloptwrap",
        calc.derivs = FALSE
      )
    )

    estimated_group_means <- emmeans::emmeans(
      model,
      specs = group,
      lmer.df = "asymptotic"
    )

    emm_table <- as.data.frame(estimated_group_means)
    group_levels <- as.character(emm_table[[group]])

    contrast_weights <- get_contrast_weights(group_levels)

    contrast_result <- emmeans::contrast(
      estimated_group_means,
      method = stats::setNames(
        list(contrast_weights),
        pairwise_contrasts
      ),
      adjust = "none"
    )

    contrast_summary <- as.data.frame(
      summary(contrast_result)
    )

    wanted <- contrast_summary[
      contrast_summary[["contrast"]] == pairwise_contrasts,
      ,
      drop = FALSE
    ]

    if (nrow(wanted) != 1L) {
      return(c(d = NA_real_, s = NA_real_))
    }

    c(
      d = wanted[["estimate"]][1L],
      s = wanted[["SE"]][1L]
    )

  }, error = function(e) {

    # Replace with a structured log if you need failure diagnostics
    c(d = NA_real_, s = NA_real_)
  })

  result
}

# Sequential within a dataset.
# Use this when datasets are already processed in parallel.
fit_results <- lapply(protein_names, fit_one_protein)

fit_matrix <- do.call(rbind, fit_results)

coefficient <- stats::setNames(
  fit_matrix[, "d"],
  protein_names
)

standard_error <- stats::setNames(
  fit_matrix[, "s"],
  protein_names
)

return(list(
  d = coefficient,
  s = standard_error
))
}


#' Fit Peptide Bootstrap Samples
#'
#' Combine group-specific bootstrap matrices and fit peptide-level mixed
#' models, returning a statistic and standard error for each protein.
#'
#' @param x List of feature-by-sample matrices, one per group in the bootstrap
#' replicate. Matrices must have matching peptide IDs as row names.
#' @param group Character string naming the column in `meta.info` that
#' identifies each sample's group.
#' @param meta.info Data frame of sample-level metadata, with sample names as
#' row names.
#' @param model_formula Formula for the peptide-level mixed model.
#' @param pairwise_contrasts Character vector of contrast expressions, such as
#' `"group1 - group2"`, to fit.
#' @param peptide_metadata Data frame mapping each `peptide_id` to a
#' `protein_id`.
#'
#' @return A list containing absolute protein-level statistics `d` and their
#' standard errors `s`.
#'
fit_bootstrap_peptide <- function(x, group, meta.info, model_formula,
                                  pairwise_contrasts, peptide_metadata) {
  
  ### Generate combined data and covariates for model fitting
  data_for_fit <- .generate_data_for_fit_peptide(x = x, meta.info = meta.info)
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
    cli_abort("peptide_metadata contains duplicated peptide_id values.")
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
    cli_abort("covariates.p contains duplicated sample_id values.")
  }
  
  if (!all(long_data[["sample_id"]] %in% covariates.p[["sample_id"]])) {
    cli_abort("Some sample IDs in long_data are missing from covariates.p.")
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
    cli_abort(
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
  #============================================================
# Prepare columns only once
#============================================================

long_data[, (group) := factor(get(group))]
long_data[, peptide_id := factor(peptide_id)]
long_data[, sample_id := factor(sample_id)]

# Avoid repeatedly scanning the complete long_data table
data.table::setkey(long_data, protein_id)

protein_names <- unique(long_data[["protein_id"]])

# Cache contrasts because most proteins have the same group levels
contrast_cache <- new.env(parent = emptyenv())

get_contrast_weights <- function(group_levels) {

  cache_key <- paste(group_levels, collapse = "\r")

  if (!exists(cache_key, envir = contrast_cache, inherits = FALSE)) {

    weights <- limma::makeContrasts(
      contrasts = stringr::str_remove_all(
        pairwise_contrasts,
        stringr::fixed(group)
      ),
      levels = group_levels
    )[, 1]

    assign(cache_key, weights, envir = contrast_cache)
  }

  get(cache_key, envir = contrast_cache, inherits = FALSE)
}

fit_one_protein <- function(current_protein) {

  protein_data <- long_data[.(current_protein)]

  # Make failures explicit rather than silently corrupting the output
  result <- tryCatch({

    present_groups <- unique(as.character(protein_data[[group]]))

    if (length(present_groups) < 2L) {
      return(c(d = NA_real_, s = NA_real_))
    }

    # droplevels is important if some proteins are missing a group
    protein_data[, (group) := droplevels(get(group))]
    protein_data[, peptide_id := droplevels(peptide_id)]
    protein_data[, sample_id := droplevels(sample_id)]

    model <- lme4::lmer(
      formula = model_formula,
      data = protein_data,
      REML = TRUE,
      control = lme4::lmerControl(
        optimizer = "nloptwrap",
        calc.derivs = FALSE
      )
    )

    estimated_group_means <- emmeans::emmeans(
      model,
      specs = group,
      lmer.df = "asymptotic"
    )

    emm_table <- as.data.frame(estimated_group_means)
    group_levels <- as.character(emm_table[[group]])

    contrast_weights <- get_contrast_weights(group_levels)

    contrast_result <- emmeans::contrast(
      estimated_group_means,
      method = stats::setNames(
        list(contrast_weights),
        pairwise_contrasts
      ),
      adjust = "none"
    )

    contrast_summary <- as.data.frame(
      summary(contrast_result)
    )

    wanted <- contrast_summary[
      contrast_summary[["contrast"]] == pairwise_contrasts,
      ,
      drop = FALSE
    ]

    if (nrow(wanted) != 1L) {
      return(c(d = NA_real_, s = NA_real_))
    }

    c(
      d = abs(wanted[["estimate"]][1L]),
      s = wanted[["SE"]][1L]
    )

  }, error = function(e) {

    # Replace with a structured log if you need failure diagnostics
    c(d = NA_real_, s = NA_real_)
  })

  result
}

# Sequential within a dataset.
# Use this when datasets are already processed in parallel.
fit_results <- lapply(protein_names, fit_one_protein)

fit_matrix <- do.call(rbind, fit_results)

coefficient <- stats::setNames(
  fit_matrix[, "d"],
  protein_names
)

standard_error <- stats::setNames(
  fit_matrix[, "s"],
  protein_names
)

return(list(
  d = coefficient,
  s = standard_error
))

}

#' Fit Peptide Models with Permuted Metadata
#'
#' Fit peptide-level mixed models to the original expression matrix using
#' permuted sample metadata to obtain null statistics.
#'
#' @param x Feature-by-sample matrix with sample names as column names and
#' peptide IDs as row names.
#' @param group Character string naming the column in `meta.info` that
#' identifies each sample's group.
#' @param meta.info Data frame of permuted sample-level metadata, with sample
#' names as row names.
#' @param model_formula Formula for the peptide-level mixed model.
#' @param pairwise_contrasts Character vector of contrast expressions, such as
#' `"group1 - group2"`, to fit.
#' @param peptide_metadata Data frame mapping each `peptide_id` to a
#' `protein_id`.
#'
#' @return A list containing absolute protein-level null statistics `d` and
#' their standard errors `s`.
#'
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
    cli_abort("peptide_metadata contains duplicated peptide_id values.")
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
    cli_abort("covariates.p contains duplicated sample_id values.")
  }
  
  if (!all(long_data[["sample_id"]] %in% covariates.p[["sample_id"]])) {
    cli_abort("Some sample IDs in long_data are missing from covariates.p.")
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
    cli_abort(
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
  #============================================================
# Prepare columns only once
#============================================================

long_data[, (group) := factor(get(group))]
long_data[, peptide_id := factor(peptide_id)]
long_data[, sample_id := factor(sample_id)]

# Avoid repeatedly scanning the complete long_data table
data.table::setkey(long_data, protein_id)

protein_names <- unique(long_data[["protein_id"]])

# Cache contrasts because most proteins have the same group levels
contrast_cache <- new.env(parent = emptyenv())

get_contrast_weights <- function(group_levels) {

  cache_key <- paste(group_levels, collapse = "\r")

  if (!exists(cache_key, envir = contrast_cache, inherits = FALSE)) {

    weights <- limma::makeContrasts(
      contrasts = stringr::str_remove_all(
        pairwise_contrasts,
        stringr::fixed(group)
      ),
      levels = group_levels
    )[, 1]

    assign(cache_key, weights, envir = contrast_cache)
  }

  get(cache_key, envir = contrast_cache, inherits = FALSE)
}

fit_one_protein <- function(current_protein) {

  protein_data <- long_data[.(current_protein)]

  # Make failures explicit rather than silently corrupting the output
  result <- tryCatch({

    present_groups <- unique(as.character(protein_data[[group]]))

    if (length(present_groups) < 2L) {
      return(c(d = NA_real_, s = NA_real_))
    }

    # droplevels is important if some proteins are missing a group
    protein_data[, (group) := droplevels(get(group))]
    protein_data[, peptide_id := droplevels(peptide_id)]
    protein_data[, sample_id := droplevels(sample_id)]

    model <- lme4::lmer(
      formula = model_formula,
      data = protein_data,
      REML = TRUE,
      control = lme4::lmerControl(
        optimizer = "nloptwrap",
        calc.derivs = FALSE
      )
    )

    estimated_group_means <- emmeans::emmeans(
      model,
      specs = group,
      lmer.df = "asymptotic"
    )

    emm_table <- as.data.frame(estimated_group_means)
    group_levels <- as.character(emm_table[[group]])

    contrast_weights <- get_contrast_weights(group_levels)

    contrast_result <- emmeans::contrast(
      estimated_group_means,
      method = stats::setNames(
        list(contrast_weights),
        pairwise_contrasts
      ),
      adjust = "none"
    )

    contrast_summary <- as.data.frame(
      summary(contrast_result)
    )

    wanted <- contrast_summary[
      contrast_summary[["contrast"]] == pairwise_contrasts,
      ,
      drop = FALSE
    ]

    if (nrow(wanted) != 1L) {
      return(c(d = NA_real_, s = NA_real_))
    }

    c(
      d = abs(wanted[["estimate"]][1L]),
      s = wanted[["SE"]][1L]
    )

  }, error = function(e) {

    # Replace with a structured log if you need failure diagnostics
    c(d = NA_real_, s = NA_real_)
  })

  result
}

# Sequential within a dataset.
# Use this when datasets are already processed in parallel.
fit_results <- lapply(protein_names, fit_one_protein)

fit_matrix <- do.call(rbind, fit_results)

coefficient <- stats::setNames(
  fit_matrix[, "d"],
  protein_names
)

standard_error <- stats::setNames(
  fit_matrix[, "s"],
  protein_names
)

return(list(
  d = coefficient,
  s = standard_error
))

  
}


#' Assemble Bootstrap Data and Matching Metadata
#'
#' Combine group-specific matrices and construct a metadata row for each
#' resampled sample while preserving their alignment.
#'
#' @param x List of feature-by-sample matrices, one per group in the bootstrap
#' replicate. Matrices must have matching peptide IDs as row names.
#' @param meta.info Data frame of sample-level metadata, with sample names as
#' row names.
#'
#' @return A list containing the combined feature-by-sample matrix
#' `combined_data` and corresponding sample-level metadata `covariates.p`.
#'
.generate_data_for_fit_peptide <- function(x, meta.info) {
  ### Check if x is a list of matrices with identical rownames
  stopifnot(
    all(vapply(
      x,
      function(z) identical(rownames(z), rownames(x[[1L]])),
      logical(1)
    ))
  )
  combined_data <- do.call(cbind, unname(x))
  
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