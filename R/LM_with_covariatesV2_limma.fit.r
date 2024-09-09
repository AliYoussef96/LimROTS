#' testStatistic_with_covariates_fit
#'
#' Performs statistical testing with covariates adjustment using linear modeling and permutation.
#'
#' @param data A list of numeric matrices containing expression data for two groups.
#' @param group.name Column name in covariates indicating group labels (optional).
#' @param covariates Data frame containing covariates (Confounding Variables).
#' @param formula.str Formula string for covariate adjustment.
#' @param trend Logical, whether to adjust for trend in the linear model (default = TRUE).
#' @param robust Logical, whether to use robust fitting in the linear model (default = TRUE).
#'
#' @return A list with three components:
#'   \item{d}{A numeric vector of absolute log fold changes.}
#'   \item{s}{A numeric vector of standard errors.}
#'   \item{corrected.logfc}{A numeric vector of log fold changes adjusted for covariates.}
#' @export
#' @importFrom limma lmFit makeContrasts contrasts.fit eBayes topTable
#' @importFrom stats formula model.matrix


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



