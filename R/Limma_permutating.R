#' Perform Permutation-Based Linear Modeling with Covariates using Limma
#'
#' This function performs linear modeling using the Limma package with
#' permutation of the covariates to evaluate the test statistics under random
#' assignments. It handles two-group comparisons and multi-group settings.
#'
#' @param x A data matrices where rows represent
#' features (e.g., genes, proteins) and columns represent samples. The list
#' should contain at least two matrices for pairwise group comparison.
#' @param group.name A character string indicating the name of the group
#' variable in `meta.info` to be used in the analysis.
#' @param meta.info A data frame containing the metadata for the samples.
#' This includes sample grouping and any covariates to be included in the model.
#' @param formula.str A string specifying the formula to be used in model
#' fitting. It should follow the standard R formula syntax
#' (e.g., `~ covariate1 + covariate2`).
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
#' @importFrom limma makeContrasts lmFit contrasts.fit eBayes topTable
#' @importFrom stringr str_split_fixed fixed
#' @importFrom utils combn
#'
#'
#'
Limma_permutating <- function(x, group.name, meta.info, formula.str
                                        ) {
    combined_data <- x
    covariates.p <- meta.info
    covariates.p$sample.id <- NULL
    row.names(covariates.p) <- NULL
    design.matrix <-
        model.matrix(formula(formula.str), data = covariates.p)
    colnames(design.matrix) <-
        make.names(colnames(design.matrix))
    fit <- lmFit(combined_data, design.matrix)
    if (length( unique( covariates.p[,group.name] ) ) == 2) {
        pairwise_contrasts <-
            paste0(group.name, unique(covariates.p[, group.name]))
        pairwise_contrasts <- combn(pairwise_contrasts, 2, function(x) {
            paste(x[1], "-", x[2])
        })
        cont_matrix <- makeContrasts(
            contrasts = pairwise_contrasts,
            levels = design.matrix
        )
        fit2 <- contrasts.fit(fit, cont_matrix)
        fit.ebayes <- eBayes(fit2, trend = FALSE, robust = FALSE)
        d_values <- topTable(fit.ebayes,
            coef = pairwise_contrasts,
            number = "Inf", sort.by = "none"
        )
        d_values <- abs(d_values$logFC)
        s_values <- as.numeric(sqrt(fit.ebayes$s2.post) *
            fit.ebayes$stdev.unscaled[, 1])
        return(list(d = d_values, s = s_values))
    } else if (length( unique( covariates.p[,group.name] ) ) > 2 & 
                                                    ncol(covariates.p) == 1) {
        pairwise_contrasts <-
            paste0(group.name, unique(covariates.p[, group.name]))
        pairwise_contrasts <- combn(pairwise_contrasts, 2, function(x) {
            paste(x[1], "-", x[2])
        })
        cont_matrix <- makeContrasts(
            contrasts = pairwise_contrasts,
            levels = design.matrix
        )
        fit2 <- contrasts.fit(fit, cont_matrix)
        fit.ebayes <- eBayes(fit2, trend = FALSE, robust = FALSE)
        msr <- fit.ebayes$F * fit.ebayes$s2.post
        return(list(d = msr, s = fit.ebayes$s2.post))
    }
}
