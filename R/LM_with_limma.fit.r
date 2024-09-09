#' Perform Linear Modeling with Covariates using Limma
#'
#' This function performs linear modeling using the Limma package, incorporating covariates
#' in the model fitting process. It is designed to handle both two-group comparisons and
#' multi-group settings with covariates.
#'
#' @param data A list containing two or more data matrices where rows represent features
#' (e.g., genes, proteins) and columns represent samples. The list should contain at least
#' two matrices for pairwise group comparison.
#' @param group.name A character string indicating the name of the group variable in
#' `meta.info` to be used in the analysis.
#' @param meta.info A data frame containing the metadata for the samples. This includes
#' sample grouping and any covariates to be included in the model.
#' @param formula.str A string specifying the formula to be used in model fitting. It should
#' follow the standard R formula syntax (e.g., `~ covariate1 + covariate2`).
#' @param trend A logical value indicating whether to allow for an intensity-dependent trend in
#' the prior variance.
#' @param robust A logical value indicating whether to use a robust fitting procedure to protect
#' against outliers.
#'
#' @details
#' This function combines the data matrices from different groups and fits a linear model
#' using covariates provided in the `meta.info`. For two-group comparisons, the function
#' computes contrasts between the two groups and applies empirical Bayes moderation. For
#' multi-group analysis with a single covariate, pairwise contrasts are computed, and the
#' moderated F-statistic is calculated for each feature.
#'
#' @return A list containing the following elements:
#' \item{d}{A vector of the test statistics (log-fold changes or F-statistics) for each feature.}
#' \item{s}{A vector of the standard deviations for each feature, adjusted by the empirical Bayes
#' procedure.}
#' \item{corrected.logfc}{The log-fold changes for each feature after fitting the model.}
#'
#' @seealso \code{\link[limma]{lmFit}}, \code{\link[limma]{eBayes}}, \code{\link[limma]{topTable}},
#' \code{\link[limma]{makeContrasts}}
#'


testStatistic_with_covariates_Fit <- function(data, group.name, meta.info , formula.str,
                                                      trend, robust) {



  combined_data <- cbind(data[[1]], data[[2]])


  covariates.p <- meta.info


  covariates.p$sample.id <- NULL


  design.matrix <- model.matrix(formula(formula.str), data = covariates.p)

  colnames(design.matrix) <- make.names(colnames(design.matrix) )


  fit <- limma::lmFit(combined_data, design.matrix)

  if(length(data) == 2){

  cont_matrix <- limma::makeContrasts("groups1-groups2",  levels=design.matrix)

  fit2 <- limma::contrasts.fit(fit, cont_matrix)

  fit.ebayes <- limma::eBayes(fit2, trend=trend, robust=robust)

  d_values <- limma::topTable(fit.ebayes, coef="groups1-groups2" , number = "Inf", sort.by = 'none')

  corrected.logfc <- d_values$logFC

  d_values <- abs(d_values$logFC)

  s_values <- as.numeric( sqrt(fit.ebayes$s2.post)  * fit.ebayes$stdev.unscaled[,1] )



  return(list(d = d_values, s = s_values , corrected.logfc = corrected.logfc))

  }else if(length(data) > 2 & ncol(covariates.p) == 1){

    pairwise_contrasts <- combn(colnames(design.matrix), 2, function(x) paste(x[1], "-", x[2]))

    cont_matrix <- limma::makeContrasts(contrasts = pairwise_contrasts, levels = design.matrix)

    fit2 <- limma::contrasts.fit(fit, cont_matrix)

    fit.ebayes <- limma::eBayes(fit2, trend=trend, robust=robust)

    msr <- fit.ebayes$F * fit.ebayes$s2.post

    corrected.logfc <- fit.ebayes$coefficients

    return(list(d = msr, s = fit.ebayes$s2.post , corrected.logfc = corrected.logfc))

  }

}



