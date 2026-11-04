# Copyright (C) 2025 Ali Youssef.
# This program is licensed under GPL (>= 2).
#
#' Parallel processing handling function for LimROTS survival
#'
#' @param data A \code{SummarizedExperiment} object or a matrix where rows
#' represent features (e.g., genes, proteins) and columns represent samples.
#' The values should be log-transformed.
#' @param a1 Optional numeric value used in the optimization process.
#' If defined by the user, no optimization occurs.
#' @param a2 Optional numeric value used in the optimization process.
#' If defined by the user, no optimization occurs.
#' @param meta.info A data frame containing sample-level metadata, where each
#' row corresponds to a sample. It should include the grouping variable
#' specified in \code{group}. If \code{x} is a \code{SummarizedExperiment}
#' object, \code{meta.info} must be a vector of the metadata needed for the
#' model to run and can be retrieved using \code{colData()}.
#' @param BPPARAM A parallel BPPARAM object for distributed computation.
#' @param formula.str A formula string used when covariates are present in meta.
#' info for modeling. It should include "~ 0 + ..." to exclude the
#' intercept from the model.
#' @param samples bootstrapped samples matrix
#' @param pSamples a permutated list of samples
#' @param competing_risks Logical. If \code{TRUE}, the competing
#' risks model via \code{crr} from \code{cmprsk} is used instead of the
#' standard Cox proportional hazards model.
#'
#'
#' @return A list containing: \code{D, S, pD, pS} for bootstrapped data and
#'  for permuted data.
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


# Copyright (C) 2025 Ali Youssef.
# This program is licensed under GPL (>= 2).
#
#' Final Per-Feature Survival Model Fit on the Full Dataset
#'
#' This function fits a per-feature survival model to the full (non-resampled)
#' data matrix using the observed sample metadata, producing the final
#' statistics used for ranking features.
#'
#' @param x A data matrix where rows represent features (e.g., proteins,
#' metabolites) and columns represent samples.
#' @param meta.info A data frame containing the metadata for the samples.
#' Must include \code{time}, \code{event}, and any additional covariates
#' used in \code{formula.str}.
#' @param formula.str A string specifying the formula to be used in model
#' fitting. Must include a \code{Surv(time, event)} term. The per-feature
#' coefficient term (\code{y}) is prepended automatically.
#' @param competing_risks Logical. If \code{FALSE} (default), a Cox
#' proportional hazards model is fitted per feature using \code{coxph}.
#' If \code{TRUE}, a competing risks model is fitted per
#' feature using \code{crr} from the \code{cmprsk} package.
#'
#' @details
#' For each feature (row), the function appends the feature expression values
#' as \code{y} to the sample metadata and fits either a Cox proportional
#' hazards model (\code{competing_risks = FALSE}) or a
#' subdistribution hazard model (\code{competing_risks = TRUE}). The
#' coefficient, its standard error, and the exponentiated coefficient
#' (hazard ratio) for the feature term \code{y} are extracted.
#'
#' Unlike \code{permutating_survival} and \code{bootstrap_survival}, this
#' function operates on the observed (non-permuted, non-resampled) data to
#' produce the final statistics used for feature ranking.
#'
#' @return A list containing the following elements:
#' \item{d}{A numeric vector of absolute coefficients
#' (\eqn{|\beta|}) for each feature.}
#' \item{s}{A numeric vector of standard errors of the coefficients for
#' each feature.}
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
            beta <- fit$coef[1]
            var  <- fit$var[1, 1]
            
            d[j] <- beta
            s[j]   <- sqrt(var)
            exp_coef[j] <- exp(fit$coef[1])
            taxon_names[j] <- row.names(combined_data[j,])
            
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

            beta <- fit$coef["y"]
            var  <- fit$var[1, 1]

            d[j] <- beta
            s[j]   <- sqrt(var)
            exp_coef[j] <- exp(beta)
            taxon_names[j] <- row.names(combined_data[j,])

        }

        names(d) <- taxon_names
        names(s) <- taxon_names
        names(exp_coef) <- taxon_names


        return(list(d = abs(d), s = s , exp_coef = as.numeric(exp_coef)))

        }
}



# Copyright (C) 2025 Ali Youssef.
# This program is licensed under GPL (>= 2).
#
#' Perform Per-Feature Survival Modeling on Bootstrap Resamples
#'
#'
#' @param x A list of data matrices (bootstrap resample), where each element
#' corresponds to a resampled group. Rows represent features (e.g., proteins,
#' metabolites) and columns represent samples.
#' @param meta.info A data frame containing the metadata for the samples,
#' including \code{time}, \code{event}, and any additional covariates used
#' in \code{formula.str}.
#' @param formula.str A string specifying the formula to be used in model
#' fitting. Must include a \code{Surv(time, event)} term. The per-feature
#' coefficient term (\code{y}) is prepended automatically.
#' @param competing_risks Logical. If \code{FALSE} (default), a Cox
#' proportional hazards model is fitted per feature using \code{coxph}.
#' If \code{TRUE}, a competing risks model is fitted per
#' feature using \code{crr} from the \code{cmprsk} package.
#'
#' @details
#' For each feature (row), the function constructs a temporary data frame
#' combining the bootstrap-resampled expression values (\code{y}) with the
#' sample metadata. It then fits either a Cox proportional hazards model
#' (\code{competing_risks = FALSE}) or a Fine subdistribution
#' hazard model (\code{competing_risks = TRUE}) per feature, extracting
#' the coefficient and its standard error for the feature term \code{y}.
#'
#' @return A list containing the following elements:
#' \item{d}{A numeric vector of absolute coefficients
#' (\eqn{|\beta|}) for each feature.}
#' \item{s}{A numeric vector of standard errors of the coefficients for
#' each feature.}
#'
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
        beta <- fit$coef[1]
        var  <- fit$var[1, 1]
        
        d[j] <- beta
        s[j]   <- sqrt(var)
        taxon_names[j] <- row.names(combined_data[j,])
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
        
        beta <- fit$coef["y"]
        var  <- fit$var[1, 1]
        
        d[j] <- beta
        s[j]   <- sqrt(var)
        taxon_names[j] <- row.names(combined_data[j,])
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



# Copyright (C) 2025 Ali Youssef.
# This program is licensed under GPL (>= 2).
#
#' Perform Per-Feature Survival Modeling on Permuted Data
#'
#' This function fits a per-feature survival model to the original data matrix
#' using permuted sample metadata, generating a null distribution of test
#' statistics for downstream FDR and p-value estimation.
#'
#' @param x A data matrix where rows represent features (e.g., proteins,
#' metabolites) and columns represent samples.
#' @param meta.info A data frame of permuted sample metadata, where each row
#' corresponds to a sample. Must include \code{time}, \code{event}, and any
#' additional covariates used in \code{formula.str}.
#' @param formula.str A string specifying the formula to be used in model
#' fitting. Must include a \code{Surv(time, event)} term. The per-feature
#' coefficient term (\code{y}) is prepended automatically.
#' @param competing_risks Logical. If \code{FALSE} (default), a Cox
#' proportional hazards model is fitted per feature using \code{coxph}.
#' If \code{TRUE}, a competing risks model is fitted per
#' feature using \code{crr} from the \code{cmprsk} package.
#'
#' @details
#' For each feature (row), the function appends the feature expression values
#' as \code{y} to the permuted metadata and fits either a Cox proportional
#' hazards model (\code{competing_risks = FALSE}) or a
#' subdistribution hazard model (\code{competing_risks = TRUE}). Because
#' the metadata is permuted, the resulting statistics form the null distribution
#' used to compute empirical p-values and FDR.
#'
#' Unlike \code{bootstrap_survival}, this function operates directly on the
#' full data matrix without group splitting or resampling.
#'
#' @return A list containing the following elements:
#' \item{d}{A numeric vector of absolute Cox coefficients
#' (\eqn{|\beta|}) for each feature.}
#' \item{s}{A numeric vector of standard errors of the coefficients for
#' each feature.}
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
      beta <- fit$coef[1]
      var  <- fit$var[1, 1]
      
      d[j] <- beta
      s[j]   <- sqrt(var)
      taxon_names[j] <- row.names(combined_data[j,])
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
      
      beta <- fit$coef["y"]
      var  <- fit$var[1, 1]
      
      d[j] <- beta
      s[j]   <- sqrt(var)
      
      taxon_names[j] <- row.names(combined_data[j,])
      
    }
    
    names(d) <- taxon_names
    names(s) <- taxon_names
    
    return(list(d = abs(d), s = s))
  }
  
}