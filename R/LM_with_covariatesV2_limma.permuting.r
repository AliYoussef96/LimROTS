#' testStatistic_with_covariates
#' 
#' Performs statistical testing with covariates adjustment using linear modeling.
#'
#' @param paired Logical, indicating whether the data is paired.
#' @param data A list of numeric matrices containing expression data for two groups.
#' @param group.name Column name in covariates indicating group labels (optional).
#' @param covariates Data frame containing covariates (Confounding Variables).
#' @param formula.str Formula string for covariate adjustment.
#' @param trend Logical, whether to adjust for trend in the linear model (default = TRUE).
#' @param robust Logical, whether to use robust fitting in the linear model (default = TRUE).
#'
#' @return A list with two components:
#'   \item{d}{A numeric vector of absolute log fold changes.}
#'   \item{s}{A numeric vector of standard errors.}
#' @export
#' @importFrom limma lmFit makeContrasts contrasts.fit eBayes topTable
#' @importFrom utils str_split_fixed
#' @importFrom stats formula model.matrix


testStatistic_with_covariates_permutating <- function(paired, data, group.name, covariates, formula.str ,
                                          trend, robust) {
  # Combine the two data sets
  combined_data <- cbind(data[[1]], data[[2]])
  
  
  covariates.p <- data.frame()
  
  for(i in colnames(combined_data)){
    
    real_SampleNames <-  str_split_fixed(i , fixed(".") , 2)[,1]
    
    df.temp <- covariates[covariates$sample.id %in% real_SampleNames,]
    
    df.temp$sample.id <- i
    
    covariates.p <- rbind(covariates.p , df.temp)
  }
  
  covariates.p <- covariates.p[sample(nrow(covariates.p)),]
  
  covariates.p$sample.id <- NULL
  
  row.names(covariates.p) <- NULL
  
  
  
  #covariates.p$groups <- sample(covariates.p$groups, nrow(covariates.p))
  #covariates.p <- as.data.frame(apply(covariates.p, 2, function(x) sample(x)))
  
  
  # formula.str <-  paste0(colnames(covariates.p), collapse = " + ")
  # 
  # formula.str <- paste0("~ 0 +"  ,  formula.str )
  
  design.matrix <- model.matrix(formula(formula.str), data = covariates.p)
  
  colnames(design.matrix) <- make.names(colnames(design.matrix) )
  
  # Fit the linear model
  fit <- limma::lmFit(combined_data, design.matrix)
  
  cont_matrix <- limma::makeContrasts("groups1-groups2",  levels=design.matrix)
  
  fit2 <- limma::contrasts.fit(fit, cont_matrix)
  
  fit.ebayes <- limma::eBayes(fit2, trend=trend, robust=robust)
  
  
  d_values <- limma::topTable(fit.ebayes, coef="groups1-groups2" , number = "Inf" , sort.by = 'none')
  
  d_values <- abs(d_values$logFC) 
  
  s_values <- as.numeric( sqrt(fit.ebayes$s2.post)  * fit.ebayes$stdev.unscaled[,1] )
  
  
  
  return(list(d = d_values, s = s_values))
}



