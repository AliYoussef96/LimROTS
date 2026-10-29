# Copyright (C) 2025 Ali Youssef.
# This program is licensed under GPL (>= 2).
#
#' Perform Permutation-Based Linear Modeling with Covariates using Limma
#'
#' This function performs linear modeling using the Limma package with
#' permutation of the covariates to evaluate the test statistics under random
#' assignments. It handles two-group comparisons and multi-group settings.
#'
#' @param x A data matrices where rows represent
#' features (e.g., genes, proteins) and columns represent samples. The list
#' should contain at least two matrices for pairwise group comparison.
#' @param group A character string indicating the name of the group
#' variable in `meta.info` to be used in the analysis.
#' @param meta.info A data frame containing the metadata for the samples.
#' This includes sample grouping and any covariates to be included in the model.
#' @param formula.str A string specifying the formula to be used in model
#' fitting. It should follow the standard R formula syntax
#' (e.g., `~ covariate1 + covariate2`).
#' @param correlation_block Character or NULL. The name of a column in
#' \code{meta.info} that defines correlation blocks. When not NULL,
#' within-block correlation is estimated via \code{duplicateCorrelation} and
#' accounted for in \code{lmFit}. If NULL, standard independent fitting is used.
#'
#' @details
#' This function combines the data matrices from different groups and permutes
#' the covariates from `meta.info`before fitting a linear model using Limma.
#' Permutation helps assess how the covariates behave under random conditions,
#' providing a null distribution of the test statistics. For two-group
#' comparisons, the function computes contrasts between the two groups and
#' applies empirical Bayes moderation. For multi-group analysis with a single
#' covariate, pairwise contrasts are computed, and the moderated F-statistic is
#' calculated for each feature.
#'
#' @return A list containing the following elements:
#' \item{d}{A vector of the test statistics (log-fold changes or F-statistics)
#' for each feature.}
#' \item{s}{A vector of the standard deviations for each feature, adjusted by
#' the empirical Bayes procedure.}
#'
#' @seealso \code{\link[limma]{lmFit}}, \code{\link[limma]{eBayes}},
#' \code{\link[limma]{topTable}},
#' \code{\link[limma]{makeContrasts}}
#'
#' @importFrom stats model.matrix formula
#' @importFrom dplyr bind_cols
#' @importFrom limma makeContrasts lmFit contrasts.fit 
#' @importFrom limma duplicateCorrelation
#' @importFrom stringr str_split_fixed fixed
#' @importFrom utils combn
#' @importFrom variancePartition dream makeContrastsDream
#' @importFrom BiocParallel SerialParam
#'
#' 


fit_permutation <- function(x, group, meta.info, formula.str,
                             correlation_block = NULL, model = "limma",
                             pairwise_contrasts) {
    if(model == "limma"){
        return(.Limma_permutating(x, group, meta.info, formula.str, 
        correlation_block, pairwise_contrasts))
    }else if (model == "dream") {
        return(.Dream_permutating(x, group, meta.info, formula.str, 
        correlation_block, pairwise_contrasts))
    }else{
        stop("Invalid model specified. Choose either 'limma' or 'dream'.")
    }
}

.Limma_permutating <- function(x, group, meta.info, formula.str,
                             correlation_block = NULL, pairwise_contrasts) {

    ### Generate combined data and covariates for model fitting
    combined_data <- x
    covariates.p <- meta.info
    row.names(covariates.p) <- NULL
    
    ### Create the design matrix based on the provided formula and covariates
    design.matrix <-
        model.matrix(formula(formula.str), data = covariates.p)
    colnames(design.matrix) <-
        make.names(colnames(design.matrix))

    ### Fit the linear model using Limma, accounting for correlation blocks if specified
    if (!is.null(correlation_block)) {
        corfit <- duplicateCorrelation(combined_data, design.matrix,
                                       block = 
                                         covariates.p[, correlation_block])
        fit <- lmFit(combined_data, design.matrix,
                     block = covariates.p[, correlation_block],
                     correlation = corfit$consensus.correlation)
    } else {
        fit <- lmFit(combined_data, design.matrix)
    }

    ### Compute contrasts and apply empirical Bayes moderation based on the number of groups

    cont_matrix <- makeContrasts(
        contrasts = pairwise_contrasts,
        levels = design.matrix
    )

    ### Fit the contrasts and apply empirical Bayes moderation
    fit2 <- limma::contrasts.fit(fit, cont_matrix)
    fit.ebayes <- limma::eBayes(fit2, trend = FALSE, robust = FALSE)

    d_values <- setNames(
        abs(as.numeric(fit.ebayes[["coefficients"]])),
        rownames(fit.ebayes[["coefficients"]]))

    stopifnot(identical(names(fit.ebayes[["s2.post"]]) ,  
                    row.names(fit.ebayes[["stdev.unscaled"]])))

    s_values <- setNames(
        as.numeric(sqrt(fit.ebayes[["s2.post"]]) *
                        fit.ebayes[["stdev.unscaled"]]),
        names(fit.ebayes[["s2.post"]]))

    stopifnot(identical(names(s_values), names(d_values) ) )

    return(list(d = d_values, s = s_values))
    
    }


.Dream_permutating <- function(x, group, meta.info, formula.str,
            correlation_block = NULL, pairwise_contrasts) {

    ### Generate combined data and covariates for model fitting
    combined_data <- x
    covariates.p <- meta.info
    row.names(covariates.p) <- NULL
    
    L <- variancePartition::makeContrastsDream(
        formula(formula.str), covariates.p,
        contrasts = pairwise_contrasts
    )
    fit <- variancePartition::dream(
        combined_data, formula(formula.str), covariates.p,
        L = L,
        BPPARAM = BiocParallel::SerialParam()
    )
    fit.ebayes <- variancePartition::eBayes(
        fit, trend = FALSE, robust = FALSE
    )
    d_values <- variancePartition::topTable(fit.ebayes,
        coef = pairwise_contrasts,
        number = "Inf", sort.by = "none"
    )
    d_values <- abs(d_values[["logFC"]])
    s_values <- as.numeric(sqrt(fit.ebayes[["s2.post"]]) *
        fit.ebayes[["stdev.unscaled"]][, 1])

    return(list(d = d_values, s = s_values))

}