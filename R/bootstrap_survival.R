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
#' If \code{TRUE}, a Fine\u2013Gray competing risks model is fitted per
#' feature using \code{crr} from the \code{cmprsk} package.
#'
#' @details
#' For each feature (row), the function constructs a temporary data frame
#' combining the bootstrap-resampled expression values (\code{y}) with the
#' sample metadata. It then fits either a Cox proportional hazards model
#' (\code{competing_risks = FALSE}) or a Fine\u2013Gray subdistribution
#' hazard model (\code{competing_risks = TRUE}) per feature, extracting
#' the coefficient and its standard error for the feature term \code{y}.
#'
#' @return A list containing the following elements:
#' \item{d}{A numeric vector of absolute Cox/Fine\u2013Gray coefficients
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
        data <- x
        combined_data <- data.frame(
            check.rows = FALSE,
            check.names = FALSE,
            none = rep("none", nrow(data[[1]]))
        )
        for (k.list in names(data)) {
            combined_data <- cbind(
                combined_data,
                data.frame(data[[k.list]], check.rows = FALSE, 
                        check.names = FALSE)
            )
        }
        combined_data <- combined_data[, -1]
        colnames(combined_data) <-
            paste0(colnames(combined_data), ".", seq(1, ncol(combined_data)))
        covariates.p <- data.frame()
        meta.info.temp <- meta.info
        meta.info.temp$sample.id <- row.names(meta.info.temp)
        for (i in colnames(combined_data)) {
            real_SampleNames <- str_split_fixed(i, fixed("."), 2)[, 1]
            df.temp <-
                meta.info.temp[row.names(meta.info.temp) %in% real_SampleNames
                            , ]
            df.temp$sample.id <- i
            covariates.p <- rbind(covariates.p, df.temp)
        }
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

