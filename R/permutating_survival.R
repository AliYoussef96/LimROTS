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
    covariates.p$sample.id <- NULL
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
    }else{
        for (j in seq_len(nrow(combined_data))) {

            df.temp <- covariates.p
            df.temp$y <- as.numeric(combined_data[j,])

            cov_matrix <- as.matrix(df.temp[, "y", drop = FALSE])
            extra_terms <- setdiff(
                colnames(df.temp),
                c("time", "event", "y")
            )
            if (length(extra_terms) > 0) {
                cov_matrix <- cbind(
                    cov_matrix,
                    as.matrix(df.temp[, extra_terms, drop = FALSE])
                )
            }

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