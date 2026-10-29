# Copyright (C) 2025 Ali Youssef.
# This program is licensed under GPL (>= 2).
#
#' Perform Linear Modeling with Covariates using Limma
#'
#' This function performs linear modeling using the Limma package, incorporating
#'  covariates in the model fitting process. It is designed to handle both
#'  two-group comparisons and multi-group settings with covariates.
#'
#' @param x A list containing two or more data matrices where rows represent
#' features (e.g., genes, proteins) and columns represent samples. The list
#' should contain at least two matrices for pairwise group comparison.
#' @param group A character string indicating the name of the group
#' variable in`meta.info` to be used in the analysis.
#' @param meta.info A data frame containing the metadata for the samples.
#' This includes sample grouping and any covariates to be included in the model.
#' @param formula.str A string specifying the formula to be used in model
#' fitting. It should follow the standard R formula syntax
#' (e.g., `~ covariate1 + covariate2`). If the formula contains a random
#' effects term (i.e., includes `|`, e.g., `~ 0 + group + (1|subject)`),
#' the function automatically routes to DREAM analysis via
#' \code{\link[variancePartition]{dream}}. The DREAM path is restricted to
#' two-group comparisons; passing more than two groups with a random term
#' raises an error.
#' @param trend A logical value indicating whether to allow for an
#' intensity-dependent trend in the prior variance.
#' @param robust A logical value indicating whether to use a robust fitting
#' procedure to protect against outliers.
#' @param correlation_block Character or NULL. The name of a column in
#' \code{meta.info} that defines correlation blocks. When not NULL,
#' within-block correlation is estimated via \code{duplicateCorrelation} and
#' accounted for in \code{lmFit}. If NULL, standard independent fitting is used.
#' @param pairwise_contrasts A character vector specifying the pairwise contrasts to be computed. 
#' Each contrast should be in the format "group1 - group2". for multi-group analysis, 
#' this will be generated automatically.
#' @param model A character string specifying the modeling approach to be used. 
#' limma for fixed effects and dream for random effects. Default is "limma".
#'
#' @details
#' This function combines the data matrices from different groups and fits a
#' linear model using covariates provided in the `meta.info`. For two-group
#' comparisons, the function computes contrasts between the two groups and
#' applies empirical Bayes moderation. For multi-group analysis with a single
#' covariate, pairwise contrasts are computed, and the moderated F-statistic is
#' calculated for each feature.
#'
#' When `formula.str` contains a `|` character (random effects term), the
#' function bypasses the standard limma workflow and instead uses
#' \code{\link[variancePartition]{dream}} to fit a linear mixed model. This
#' path is only supported for two-group comparisons. Passing more than two
#' groups with a random-effects formula raises an informative error directing
#' the user to use the \code{correlation_block} argument with a fixed-effects
#' formula instead.
#'
#' @return A list containing the following elements:
#' \item{d}{A vector of the test statistics (log-fold changes or F-statistics)
#' for each feature.}
#' \item{s}{A vector of the standard deviations for each feature, adjusted by t
#' he empirical Bayes procedure.}
#' \item{corrected.logfc}{The log-fold changes for each feature after fitting
#' the model.}
#'
#' @seealso \code{\link[limma]{lmFit}},
#' \code{\link[limma]{eBayes}},
#' \code{\link[limma]{topTable}},
#' \code{\link[limma]{makeContrasts}},
#' \code{\link[variancePartition]{dream}},
#' \code{\link[variancePartition]{makeContrastsDream}}
#'
#'
#'
#' @importFrom stats model.matrix formula
#' @importFrom dplyr bind_cols
#' @importFrom utils combn
#' @importFrom limma makeContrasts lmFit contrasts.fit
#' @importFrom limma duplicateCorrelation
#' @importFrom variancePartition dream makeContrastsDream
#' @importFrom BiocParallel SerialParam
#'


fit_original <- function(x, group, meta.info, formula.str, trend, robust,
                          correlation_block = NULL, model = "limma",
                          pairwise_contrasts) {
    if(model == "limma"){
        return(.Limma_fit_original(x, group, meta.info, formula.str, 
        correlation_block, pairwise_contrasts))
    }else if (model == "dream") {
        return(.Dream_fit_original(x, group, meta.info, formula.str, 
        correlation_block, pairwise_contrasts))
    }else{
        stop("Invalid model specified. Choose either 'limma' or 'dream'.")
    }
}

.Limma_fit_original <- function(x, group, meta.info, formula.str, 
            trend, robust, correlation_block = NULL, pairwise_contrasts) {
        ### Check if x is a list of matrices with identical rownames
        combined_data <- x

        if(!identical(colnames(combined_data), row.names(meta.info))){
            stop("Column names of the data matrix must match row names of meta.info.")
        }

        ### Create the design matrix based on the provided formula and covariates
        design.matrix <- model.matrix(formula(formula.str), data = meta.info)
        colnames(design.matrix) <- make.names(colnames(design.matrix))
        stopifnot(identical(colnames(combined_data) , row.names(design.matrix)))

        ### Fit the linear model using Limma, accounting for correlation blocks if specified
        if (!is.null(correlation_block)) {
            corfit <- duplicateCorrelation(combined_data, design.matrix,
                                        block = meta.info[, correlation_block])
            fit <- lmFit(combined_data, design.matrix,
                        block = meta.info[, correlation_block],
                        correlation = corfit$consensus.correlation)
        } else {
            fit <- lmFit(combined_data, design.matrix)
        }

        ### For two-group comparisons, compute the specified pairwise contrasts
        cont_matrix <- makeContrasts(
            contrasts = pairwise_contrasts, levels =
                design.matrix
        )

        ### Fit the contrasts and apply empirical Bayes moderation
        fit2 <- limma::contrasts.fit(fit, cont_matrix)
        fit.ebayes <- limma::eBayes(fit2, trend = trend, robust = robust)

        d_values <- setNames(
            abs(as.numeric(fit.ebayes[["coefficients"]])),
            rownames(fit.ebayes[["coefficients"]]))

        corrected.logfc <- setNames(
            as.numeric(fit.ebayes[["coefficients"]]),
            rownames(fit.ebayes[["coefficients"]]))

        stopifnot(identical(names(fit.ebayes[["s2.post"]]) ,  
                        row.names(fit.ebayes[["stdev.unscaled"]])))

        s_values <- setNames(
            as.numeric(sqrt(fit.ebayes[["s2.post"]]) *
                            fit.ebayes[["stdev.unscaled"]]),
            names(fit.ebayes[["s2.post"]]))

        stopifnot(identical(names(s_values), names(d_values) ) )
        stopifnot(identical(names(corrected.logfc), names(d_values) ) )
        stopifnot(identical(names(s_values), names(corrected.logfc) ) )

        return(list(
            d = d_values,
            s = s_values,
            corrected.logfc = corrected.logfc
        ))
            
    }


.Dream_fit_original <- function(x, group, meta.info, formula.str, trend, robust,
            correlation_block = NULL, pairwise_contrasts) {

        combined_data <- x

        if(!identical(colnames(combined_data), row.names(meta.info))){
            stop("Column names of the data matrix must match row names of meta.info.")
        }

        ### Create the design matrix based on the provided formula and covariates
        L <- variancePartition::makeContrastsDream(
            formula(formula.str), meta.info,
            contrasts = pairwise_contrasts
        )

        ### Fit the model using DREAM
        fit <- variancePartition::dream(
            combined_data, formula(formula.str), meta.info,
            L = L,
            BPPARAM = BiocParallel::SerialParam()
        )
        fit.ebayes <- variancePartition::eBayes(
            fit, trend = trend, robust = robust
        )
        d_values <- variancePartition::topTable(
            fit.ebayes,
            coef = pairwise_contrasts,
            number = "Inf",
            sort.by = "none"
        )
        corrected.logfc <- d_values[["logFC"]]
        d_values <- abs(d_values[["logFC"]])
        s_values <-
            as.numeric(sqrt(fit.ebayes[["s2.post"]]) *
                fit.ebayes[["stdev.unscaled"]][, 1])
        return(list(
            d = d_values,
            s = s_values,
            corrected.logfc = corrected.logfc
        ))
}