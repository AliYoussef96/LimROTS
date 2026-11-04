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