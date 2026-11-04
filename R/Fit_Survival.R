# Copyright (C) 2025 Ali Youssef.
# This program is licensed under GPL (>= 2).
#
#' Fit Survival Bootstrap and Permutation Replicates in Parallel
#'
#' Fit a survival model to each bootstrap replicate and its corresponding
#' permuted metadata using a `BiocParallel` backend.
#'
#' @param BPPARAM A `BiocParallelParam` object specifying the parallel backend.
#' If `NULL`, the function selects a two-worker backend appropriate to the
#' operating system.
#' @param samples Matrix of bootstrap replicates, with one row per replicate
#' and resampled sample names in its columns.
#' @param data Feature-by-sample matrix with sample names as column names.
#' @param formula.str Character string containing the survival model formula.
#' It must include `Surv(time,event)~y`; additional covariates must be columns
#' in `meta.info`.
#' @param meta.info Data frame of sample-level metadata, with sample names as
#' row names.
#' @param pSamples List of permuted sample-level metadata data frames, one per
#' bootstrap replicate.
#' @param competing_risks Logical. If `TRUE`, fit a Fine-Gray competing-risks
#' model with `cmprsk::crr()`; otherwise fit a Cox proportional hazards model
#' with `survival::coxph()`.
#'
#' @return A list with one element per replicate. Each element contains `ds`
#' (bootstrap statistics and standard errors) and `pdps` (permutation
#' statistics and standard errors).
#'
#' @importFrom BiocParallel SnowParam MulticoreParam bplapply
#' @importFrom stats model.matrix formula p.adjust
#'


Boot_parallel_survival <- function(BPPARAM = NULL,
        samples,
        data,
        formula.str,
        meta.info,
        pSamples,
        competing_risks) {
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
        formula.str = formula.str,
        meta.info = meta.info,
        pSamples = pSamples,
        competing_risks = competing_risks
    )
    export_funcs <- list(
        bootstrap_survival = bootstrap_survival,
        permutating_survival = permutating_survival
    )
    for (name in names(export_vars)) {
        assign(name, export_vars[[name]])
    }
    for (name in names(export_funcs)) {
        assign(name, export_funcs[[name]])
    }
    
    results_list <- bplapply(seq_len(nrow(samples)), function(i) {
        samples.R <- split(samples[i, ], rep(1, length(samples[i, ])))
        pSamples_i <- pSamples[[i]]
        d_result <- s_result <- pd_result <- ps_result <- NULL
        fit <- bootstrap_survival(
            x = lapply(samples.R, function(x)
                data[, x]
            ),
            meta.info = meta.info,
            formula.str = formula.str,
            competing_risks = competing_risks
        )
        d_result <- fit$d
        s_result <- fit$s
        df1 <- data.frame(d_result = d_result, s_result = s_result)
        pFit <- permutating_survival(
            x = data,
            meta.info = pSamples_i,
            formula.str = formula.str,
            competing_risks = competing_risks
        )
        pd_result <- pFit$d
        ps_result <- pFit$s
        df2 <- data.frame(pd_result = pd_result, ps_result = ps_result)
        list(ds = df1, pdps = df2)
    }, BPPARAM = BPPARAM,
    BPOPTIONS = bpoptions(packages = c(
        "utils", "stringr",
        "stats", "survival",
        "cmprsk")))
    return(results_list)
}


#' Fit Survival Models to the Original Data
#'
#' Fit a survival model for every feature using the original data and observed
#' sample metadata.
#'
#' @param x Feature-by-sample matrix with sample names as column names.
#' @param meta.info Data frame of sample-level metadata, with sample names as
#' row names. It must contain `time`, `event`, and any covariates in
#' `formula.str`.
#' @param formula.str Character string containing the survival model formula.
#' It must include `Surv(time,event)~y`; additional covariates must be columns
#' in `meta.info`.
#' @param competing_risks Logical. If `TRUE`, fit a Fine-Gray competing-risks
#' model with `cmprsk::crr()`; otherwise fit a Cox proportional hazards model
#' with `survival::coxph()`.
#'
#' @details
#' For each feature (row), the function appends the feature expression values
#' as `y` to the sample metadata and fits either a Cox proportional
#' hazards model (`competing_risks = FALSE`) or a Fine-Gray
#' subdistribution hazard model (`competing_risks = TRUE`). The
#' coefficient, its standard error, and the exponentiated coefficient
#' (hazard ratio) for the feature term \code{y} are extracted.
#'
#' Unlike `permutating_survival` and `bootstrap_survival`, this function
#' operates on the observed (non-permuted, non-resampled) data to
#' produce the final statistics used for feature ranking.
#'
#' @return A list containing the following elements:
#' \item{d}{A numeric vector of absolute coefficients
#' (\eqn{|\beta|}) for each feature.}
#' \item{s}{A numeric vector of coefficient standard errors for each feature.}
#' \item{exp_coef}{A numeric vector of exponentiated coefficients (hazard
#' ratios, \eqn{e^{\beta}}) for each feature.}
#'
#' @seealso \code{\link[survival]{coxph}}, \code{\link[cmprsk]{crr}}
#'
#' @importFrom stats model.matrix formula setNames
#' @importFrom dplyr bind_cols
#' @import survival
#' @importFrom utils combn
#' @importFrom cmprsk crr
#'
#'
#'
fit_survival <- function(x, meta.info, formula.str,
competing_risks
) {
    combined_data <- data.frame(x, check.rows = FALSE, check.names = FALSE)
    covariates.p <- meta.info
    covariates.p$sample.id <- NULL
    row.names(covariates.p) <- NULL
    
    d  <- numeric(nrow(combined_data))
    s    <- numeric(nrow(combined_data))
    taxon_names <- character(nrow(combined_data))
    exp_coef <- character(nrow(combined_data))
    
    if(competing_risks == FALSE){
        for (j in seq_len(nrow(combined_data))) {
            
            df.temp <- covariates.p
            df.temp$y <- as.numeric(combined_data[j,])
            
            fit <- coxph(
                formula = formula(formula.str),
                ties = "breslow", data = df.temp
            )
            
            # Extract taxon coefficient only
            d[j] <- summary(fit)$coefficients["y", "coef"]

            s[j]   <- summary(fit)$coefficients["y", "se(coef)"]

            exp_coef[j] <- summary(fit)$coefficients["y", "exp(coef)"]

            taxon_names[j] <- row.names(combined_data)[j]
            
        }
        
        names(d) <- taxon_names
        names(s) <- taxon_names
        names(exp_coef) <- taxon_names
        
        
        return(list(d = abs(d), s = s , exp_coef = as.numeric(exp_coef)))
    }else if(competing_risks == TRUE){
        for (j in seq_len(nrow(combined_data))) {

            df.temp <- covariates.p
            df.temp$y <- as.numeric(combined_data[j,])
            
            formula.str.temp <- str_remove(formula.str , fixed("Surv(time,event)~"))
            cov_matrix <- model.matrix(as.formula(paste0("~0+",formula.str.temp)), data = df.temp)
            
            fit <- crr(
                ftime   = df.temp$time,
                fstatus = df.temp$event,
                cov1    = cov_matrix,
                failcode = 1,
                cencode  = 0,
                variance = TRUE
            )


            d[j] <- summary(fit)$coef["y", "coef"]

            s[j]   <- summary(fit)$coef["y", "se(coef)"]

            exp_coef[j] <- summary(fit)$coef["y", "exp(coef)"]

            taxon_names[j] <- row.names(combined_data)[j]

        }

        names(d) <- taxon_names
        names(s) <- taxon_names
        names(exp_coef) <- taxon_names


        return(list(d = abs(d), s = s , exp_coef = as.numeric(exp_coef)))

        }
}



#' Fit Survival Models to Bootstrap Resamples
#'
#' Fit a survival model for every feature using bootstrap-resampled expression
#' values and the matching sample metadata.
#'
#' @param x List of feature-by-sample matrices, one per group in the bootstrap
#' replicate. Matrices must have matching feature names.
#' @param meta.info Data frame of sample-level metadata, with sample names as
#' row names. It must contain `time`, `event`, and any covariates in
#' `formula.str`.
#' @param formula.str Character string containing the survival model formula.
#' It must include `Surv(time,event)~y`; additional covariates must be columns
#' in `meta.info`.
#' @param competing_risks Logical. If `TRUE`, fit a Fine-Gray competing-risks
#' model with `cmprsk::crr()`; otherwise fit a Cox proportional hazards model
#' with `survival::coxph()`.
#'
#' @details
#' For each feature (row), the function constructs a temporary data frame
#' combining the bootstrap-resampled expression values (`y`) with the
#' sample metadata. It then fits either a Cox proportional hazards model
#' (`competing_risks = FALSE`) or a Fine-Gray subdistribution
#' hazard model (`competing_risks = TRUE`) per feature, extracting
#' the coefficient and its standard error for the feature term \code{y}.
#'
#' @return A list containing the following elements:
#' \item{d}{A numeric vector of absolute coefficients
#' (\eqn{|\beta|}) for each feature.}
#' \item{s}{A numeric vector of coefficient standard errors for each feature.}
#' Assemble Bootstrap Data and Matching Survival Metadata
#'
#' Combine group-specific bootstrap matrices and construct a metadata row for
#' each resampled sample, preserving the alignment between expression values
#' and sample-level survival covariates.
#'
#' @param x List of feature-by-sample matrices, one per group in the bootstrap
#' replicate. Matrices must have matching feature names.
#' @param meta.info Data frame of sample-level metadata, with sample names as
#' row names. It must contain the survival outcomes and model covariates.
#'
#' @return A list containing the combined feature-by-sample matrix
#' `combined_data` and corresponding sample-level metadata `covariates.p`.
#' @seealso \code{\link[survival]{coxph}}, \code{\link[cmprsk]{crr}}
#'
#' @importFrom stats model.matrix formula
#' @importFrom dplyr bind_cols
#' @importFrom stringr str_split_fixed fixed
#' @import survival
#' @importFrom utils combn
#' @importFrom cmprsk crr
#'

bootstrap_survival <-
  function(x, meta.info, formula.str, competing_risks) {
    
    ### Generate combined data and covariates for model fitting
    data_for_fit <- .generate_data_for_fit(x = x, meta.info = meta.info)
    combined_data <- data_for_fit$combined_data
    covariates.p <- data_for_fit$covariates.p
    
    stopifnot(identical(colnames(combined_data), row.names(covariates.p)))
    
    d  <- numeric(nrow(combined_data))
    s    <- numeric(nrow(combined_data))
    taxon_names <- character(nrow(combined_data))
    
    if(competing_risks == FALSE){
      for (j in seq_len(nrow(combined_data))) {
        df.temp <- covariates.p
        df.temp$y <- as.numeric(combined_data[j,])
        
        fit <- coxph(
          formula = formula(formula.str),
          ties = "breslow", data = df.temp
        )
        
        # Extract taxon coefficient only
        d[j] <- summary(fit)$coefficients["y", "coef"]

        s[j]   <- summary(fit)$coefficients["y", "se(coef)"]

        taxon_names[j] <- row.names(combined_data)[j]
      }
      
      names(d) <- taxon_names
      names(s) <- taxon_names
      
      return(list(d = abs(d), s = s))
      
    }else if(competing_risks == TRUE){
      for (j in seq_len(nrow(combined_data))) {
        df.temp <- covariates.p
        df.temp$y <- as.numeric(combined_data[j,])
        
        formula.str.temp <- str_remove(formula.str , fixed("Surv(time,event)~"))
        cov_matrix <- model.matrix(as.formula(paste0("~0+",formula.str.temp)), data = df.temp)
        
        fit <- crr(
          ftime   = df.temp$time,
          fstatus = df.temp$event,
          cov1    = cov_matrix,
          failcode = 1,
          cencode  = 0,
          variance = TRUE
        )
        
        d[j] <- summary(fit)$coef["y", "coef"]

        s[j]   <- summary(fit)$coef["y", "se(coef)"]

        taxon_names[j] <- row.names(combined_data)[j]
      }
      
      names(d) <- taxon_names
      names(s) <- taxon_names
      
      return(list(d = abs(d), s = s))
    }
  }


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



#' Fit Survival Models to Permuted Metadata
#'
#' Fit a survival model for every feature using the original data and permuted
#' sample metadata to generate null statistics.
#'
#' @param x Feature-by-sample matrix with sample names as column names.
#' @param meta.info Data frame of sample-level metadata, with sample names as
#' row names. It must contain `time`, `event`, and any covariates in
#' `formula.str`; the group or event labels are permuted for the null fit.
#' @param formula.str Character string containing the survival model formula.
#' It must include `Surv(time,event)~y`; additional covariates must be columns
#' in `meta.info`.
#' @param competing_risks Logical. If `TRUE`, fit a Fine-Gray competing-risks
#' model with `cmprsk::crr()`; otherwise fit a Cox proportional hazards model
#' with `survival::coxph()`.
#'
#' @details
#' For each feature (row), the function appends the feature expression values
#' as `y` to the permuted metadata and fits either a Cox proportional
#' hazards model (`competing_risks = FALSE`) or a Fine-Gray
#' subdistribution hazard model (`competing_risks = TRUE`). Because
#' the metadata is permuted, the resulting statistics form the null distribution
#' used to compute empirical p-values and FDR.
#'
#' Unlike `bootstrap_survival`, this function operates directly on the
#' full data matrix without group splitting or resampling.
#'
#' @return A list containing the following elements:
#' \item{d}{A numeric vector of absolute model coefficients
#' (\eqn{|\beta|}) for each feature.}
#' \item{s}{A numeric vector of coefficient standard errors for each feature.}
#'
#' @seealso \code{\link[survival]{coxph}}, \code{\link[cmprsk]{crr}}
#'
#' @importFrom stats model.matrix formula
#' @importFrom dplyr bind_cols
#' @import survival
#' @importFrom utils combn
#' @importFrom cmprsk crr
#'
#'
#'
permutating_survival <- function(x, meta.info, formula.str,
                                 competing_risks
) {
  combined_data <- data.frame(x, check.rows = FALSE, check.names = FALSE)
  covariates.p <- meta.info
  row.names(covariates.p) <- NULL
  
  d  <- numeric(nrow(combined_data))
  s    <- numeric(nrow(combined_data))
  taxon_names <- character(nrow(combined_data))
  
  if(competing_risks == FALSE){
    for (j in seq_len(nrow(combined_data))) {
      df.temp <- covariates.p
      df.temp$y <- as.numeric(combined_data[j,])
      
      fit <- coxph(
        formula = formula(formula.str),
        ties = "breslow", data = df.temp
      )
      
      # Extract taxon coefficient only
      d[j] <- summary(fit)$coefficients["y", "coef"]

      s[j]   <- summary(fit)$coefficients["y", "se(coef)"]

      taxon_names[j] <- row.names(combined_data)[j]
    }
    
    names(d) <- taxon_names
    names(s) <- taxon_names
    
    return(list(d = abs(d), s = s))
  }else if(competing_risks == TRUE){
    for (j in seq_len(nrow(combined_data))) {
      
      df.temp <- covariates.p
      df.temp$y <- as.numeric(combined_data[j,])
      
      formula.str.temp <- str_remove(formula.str , fixed("Surv(time,event)~"))
      cov_matrix <- model.matrix(as.formula(paste0("~0+",formula.str.temp)), data = df.temp)
      
      fit <- crr(
        ftime   = df.temp$time,
        fstatus = df.temp$event,
        cov1    = cov_matrix,
        failcode = 1,
        cencode  = 0,
        variance = TRUE
      )
      
      d[j] <- summary(fit)$coef["y", "coef"]

      s[j]   <- summary(fit)$coef["y", "se(coef)"]
      
      taxon_names[j] <- row.names(combined_data)[j]
      
    }
    
    names(d) <- taxon_names
    names(s) <- taxon_names
    
    return(list(d = abs(d), s = s))
  }
  
}