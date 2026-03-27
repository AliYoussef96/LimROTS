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
#' @importFrom stats model.matrix formula
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
    } else {
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
            exp_coef[j] <- exp(beta)
            taxon_names[j] <- row.names(combined_data[j,])

        }

        names(d) <- taxon_names
        names(s) <- taxon_names
        names(exp_coef) <- taxon_names


        return(list(d = abs(d), s = s , exp_coef = as.numeric(exp_coef)))

        }
}