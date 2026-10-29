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
#' group or conditions for comparison.
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
#'   \item \code{group}: Character vector indicating group assignments.
#'   \item \code{K}: Top list size to be used in the analysis.
#' }
#'
#'
SanityChecK <- function(x,
        assay.type,
        niter, K,
        meta.info, group, model,
        formula.str, survival, pairwise_contrasts,
        correlation_block) {

    ### Check if x is a SummarizedExperiment and extract data and metadata
    ### and extract the assay if specified
    Check_SExp <- .Check_SummarizedExperiment(
        x = x,
        assay.type = assay.type,
        meta.info = meta.info,
        group = group,
        survival = survival,
        K = K
    )
    data <- Check_SExp[["data"]]
    group <- Check_SExp[["group"]]
    meta.info <- Check_SExp[["meta.info"]]
    K <- Check_SExp[["K"]]
    
    .Check_meta_info(
        meta.info = meta.info,
        data = data,
        group = group,
        K = K
    )

    ### Check if survival analysis is requested and perform additional checks
    if(survival){
        message("Survival analysis mode activated. Performing additional checks.")
        survival_check <- .check_survival(meta.info = meta.info, 
        data = data, 
        K = K, 
        formula.str = formula.str)
        return(survival_check)
      }


    ### Checks for the model type and random effects in the formula
    .check_random_effects(formula.str = formula.str, model = model, correlation_block = correlation_block)

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

    .Check_meta_info(
        meta.info = meta.info,
        data = data,
        group = group,
        K = K
    )

    return(list(
        meta.info = meta.info,
        data = data,
        group = group,
        K = K
    ))
}


# Copyright (C) 2025 Ali Youssef.
# This program is licensed under GPL (>= 2).
#
#' Check if SummarizedExperiment or data is correct
#'
#' @param x A matrix-like object or a `SummarizedExperiment` containing
#' the data to be analyzed.
#' @param assay.type A character string or numeric index specifying the assay
#' to use if `x` is a `SummarizedExperiment`. Default is `NULL`
#' @param meta.info Data frame. Metadata associated with the samples
#' (columns of `x`). If `x` is a `SummarizedExperiment`,
#' @param group Character. Column name in `meta.info` that defines the
#' group or conditions for comparison.
#' @param survival Logical, indicating whether the analysis is survival
#' analysis. Default is \code{FALSE}.
#'
#' @import SummarizedExperiment
#'
#' @return a list of `data` , `group` and `meta.info`
#'
#'
#'

.Check_SummarizedExperiment <- function(x,
                                       assay.type, meta.info,
                                       group, survival, K) {

    
    if(!inherits(x, "SummarizedExperiment")){
        stop("Input 'x' must be a SummarizedExperiment object.")
    }

    if (inherits(x, "SummarizedExperiment")) {
        message("Data is SummarizedExperiment object")
        if (is.null(meta.info)) {
            stop("meta.info should be a vector of colData names to be used")
        } else {
            meta.info.colnames <- meta.info
            meta.info <- data.frame(
                colData(x)[, meta.info ,  drop = FALSE],
                check.names = FALSE,
                row.names = row.names(colData(x))
            )
            if (!identical(colnames(meta.info) , meta.info.colnames)) {
                stop("meta.info should be a vector of colData names to be used")
            }
        }
        if(!survival){
            if (!group %in% colnames(meta.info)) {
                stop(
                    "group should be a string specifying",
                    " the column in `meta.info` that",
                    " represents the group or conditions",
                    " for comparison."
                )
            }
        }
        if (is.null(assay.type)) {
            assay.type <- assayNames(x)[1]
        }
        message(sprintf("Assay: %s will be used", assay.type))
        data <- assay(x, assay.type)
    } else {
        stop("Input 'x' must be a SummarizedExperiment object.")
    }

    ### Check if K is NULL and set it to a quarter of the number of rows in the data
    if (is.null(K)) {
        K <- floor(nrow(data) / 4)
        message("Top list size (K) not provided. 
        Setting K to a quarter of the number of rows in the data: ", K)
    }

    return(list(
        data = data, group = group,
        meta.info = meta.info,
        K = K
    ))
}

# Copyright (C) 2025 Ali Youssef.
# This program is licensed under GPL (>= 2).
#
#' Check if meta info is correct
#'
#' @param meta.info Data frame. Metadata associated with the samples
#' (columns of `x`). If `x` is a `SummarizedExperiment`,
#' @param data A matrix-like object or a `SummarizedExperiment` containing the
#' data to be analyzed.
#'
#' @return Logical
#'
#'

.Check_meta_info <- function(meta.info, data , group, K) {
    if (any(!row.names(meta.info) %in% colnames(data))) {
        stop("rownames
            for colData should match the data colnames (samples names)")
    }
    if (any(grepl(".", colnames(data), fixed = TRUE))) {
        stop("Sample names should contains no '.', please remove it if any")
    }
    if (nrow(meta.info) != ncol(data)) {
        stop("Number of samples in the data does not match the group")
    }
    if(!identical(sort(rownames(meta.info)), sort(colnames(data)))){
        stop("Sample names in the data and meta info in colData do not match")
    }
    if (is.null(row.names(meta.info)) ||
        anyDuplicated(row.names(meta.info))) {
        stop("meta.info must have unique sample names as row names.")
    }
    if(length(unique(meta.info[, group])) != 2){
        stop("LimROTS only supports two-group comparisons.")
    }
    if (anyNA(meta.info[, group])) {
        stop("The group column contains missing values.")
    }
    if(K < 1 || K > nrow(data)){
        stop("K must be between 1 and the number of rows in the data matrix.")
    }
    return(TRUE)
}

# Copyright (C) 2025 Ali Youssef.
# This program is licensed under GPL (>= 2).
#' Check if survival analysis is correct
#' 
#' @param meta.info Data frame. Metadata associated with the samples
#' (columns of `x`). If `x` is a `SummarizedExperiment
#' @param data A matrix-like object or a `SummarizedExperiment` containing the
#' data to be analyzed.
#' @param K Integer. Top list size. If NULL, it will be set to a quarter of
#' the number of rows in the data matrix. Default is NULL.
#' @param formula.str Optional character string representing the formula for
#' the model.
#' 
#' @return A list containing:
#' \itemize{meta.info: Processed metadata.
#' \itemize{data: Processed data matrix.
#' \itemize{K: Top list size to be used in the analysis.
#' \itemize{formula.str: Processed formula string for the model.
#' 
#' 


.check_survival <- function(meta.info, data, K, formula.str) {
      
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
      
    meta.info <- meta.info[order(meta.info[,event,  drop = FALSE]), ]
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

.check_random_effects <- function(formula.str, model, correlation_block) {

        #### Check if correlation_block is provided and is a single character string
        if(!is.null(correlation_block)){
            if(length(correlation_block) != 1){
                stop("correlation_block should be a single character string specifying the 
                column name in meta.info that defines correlation blocks.")
            }
            if (!correlation_block %in% colnames(meta.info)) {
            stop(
                "The correlation-block column '",
                correlation_block,
                "' was not found in meta.info."
            )
        }
        }

        ### Check for random effects in the formula
        random_terms <- lme4::findbars(formula(formula.str))
        if(model == "limma"){
            if (!is.null(random_terms)) { 
                stop("random terms analysis with random effects is not supported in Limma.
                Please use model='dream' for random effects analysis.")
            }
        }else if(model == "dream"){
            if (is.null(random_terms)) { 
                stop("dream need a random term(s) in the formula.str")
            }
            if(is.null(correlation_block)){
                stop("Please specify which random effect term to use via correlation_block")
            }
            any(vapply(random_terms, function(x) {
                rlang::as_string(rlang::get_expr(x[[3]])) == correlation_block
            }, logical(1)))
        }
}
