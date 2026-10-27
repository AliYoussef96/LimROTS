# Copyright (C) 2025 Ali Youssef.
# This program is licensed under GPL (>= 2).
#
#' Sanity Check for Input Data and Parameters
#'
#' This function performs a series of checks and initial setups for input data,
#' metadata, and parameters, ensuring everything is correctly formatted for
#' downstream analysis.
#'
#' @param x A matrix-like object or a `SummarizedExperiment` containing the
#' data to be analyzed.
#' @param assay.type A character string or numeric index specifying the assay
#' to use if `x` is a `SummarizedExperiment`. Default is `NULL`.
#' @param niter Integer. Number of bootstrap samples or resampling iterations.
#' Default is 1000.
#' @param K Integer. Top list size. If NULL, it will be set to a quarter of
#' the number of rows in the data matrix. Default is NULL.
#' @param meta.info Data frame. Metadata associated with the samples
#' (columns of `x`). If `x` is a `SummarizedExperiment`,
#' `meta.info` can be a vector of `colData` column names to use.
#' @param group Character. Column name in `meta.info` that defines the
#' groups or conditions for comparison.
#' @param formula.str Optional character string representing the formula for
#' the model.
#' @param survival Logical, indicating whether the analysis is survival 
#' analysis.
#' Default is \code{FALSE}. If \code{TRUE}, \code{meta.info} must contain
#' \code{time} and \code{event} columns.
#'
#' @details
#' This function checks whether the input data and metadata are in the correct
#' format, processes metadata from a `SummarizedExperiment` object if provided,
#' and ensures that group information is correctly specified. If no top list
#' size (`K`) is provided, it defaults to a quarter of the number of rows in
#' the data.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{meta.info}: Processed metadata.
#'   \item \code{data}: Processed data matrix.
#'   \item \code{groups}: Numeric or factor vector indicating group assignments.
#'   \item \code{K}: Top list size to be used in the analysis.
#' }
#'
#'
SanityChecK <- function(x,
        assay.type = NULL,
        niter = 1000, K = NULL,
        meta.info, group,
        formula.str, survival = FALSE) {

    
    ### Check if x is a SummarizedExperiment and extract data and metadata
    ### and extract the assay if specified
    Check_SExp <- Check_SummarizedExperiment(
        x = x,
        assay.type = assay.type,
        meta.info = meta.info,
        group = group,
        survival = survival
    )
    data <- Check_SExp[["data"]]
    groups <- Check_SExp[["groups"]]
    meta.info <- Check_SExp[["meta.info"]]

    Check_meta_info(
        meta.info = meta.info,
        data = data
    )
  
    ### Check if survival analysis is requested and perform additional checks
    if(survival){
        message("Survival analysis mode activated. Performing additional checks.")
        survival_check <- .generate_survival(meta.info = meta.info, 
        data = data, 
        K = K, 
        formula.str = formula.str)
        return(survival_check)
      }

    ### Final checks and sorting of meta.info and data based on the group column

    meta.info <- meta.info[order(meta.info[,group]),]
    data <- data[, row.names(meta.info)]
    if(!identical(colnames(data) , row.names(meta.info)) ){
        stop("After sorting, column names of data do not match row names of
             meta info in colData")
    }

    if(!inherits(meta.info[, group], "character")) {
        stop("The group column in meta info from colData must be of type character.")
    }

    ### Check if K is NULL and set it to a quarter of the number of rows in the data
    if (is.null(K)) {
        K <- floor(nrow(data) / 4)
        message("Top list size (K) not provided. Setting K to a quarter of the number of rows in the data: ", K)
    }
    return(list(
        meta.info = meta.info,
        data = data,
        groups = groups,
        K = K
    ))
}



.generate_survival <- function(meta.info, data, K, formula.str) {
      
      if(!all(c("time", "event") %in% colnames(meta.info))){
        stop("For survival analysis, 'meta.info' must contain 'time' and 
             'event' columns.")
      }
      
      if(!grepl("Surv\\(time, event\\)", formula.str)){
        stop("For survival analysis, 'formula.str' must include a 
             'Surv(time, event)' term.")
      }
      
      if (is.null(K)) {
        K <- floor(nrow(data) / 4)
      }
      
      formula.str <- gsub("\\s+", "", formula.str)
      
      terms <- str_split_fixed(formula.str, fixed("~"), 2)
      rhs <- terms[, 2]
      surv_part <- regmatches(rhs, regexpr("Surv\\([^)]+\\)", rhs))
      rhs_covariates <- gsub("Surv\\([^)]+\\)\\+?|\\+?Surv\\([^)]+\\)",
                             "", rhs)
      rhs_covariates <- gsub("^\\+|\\+$", "", rhs_covariates)
      if (rhs_covariates == "") {
        formula.str <- paste0(surv_part, "~y")
      } else {
        formula.str <- paste0(surv_part, "~y+", rhs_covariates)
      }
      
    meta.info <- meta.info[order(meta.info[,event]),]
    data <- data[, row.names(meta.info)]
    if(!identical(colnames(data) , row.names(meta.info)) ){
        stop("After sorting, column names of data do not match row names of
             meta info in colData")
    }
      
      if(!identical(colnames(data), row.names(meta.info))){
        stop("After sorting, column names of data do not match row names of
             meta.info.")
      }
      message("Sanity check completed successfully!")
      
      return(list(
        meta.info = meta.info,
        data = data,
        K = K,
        formula.str = formula.str
      ))
    
}