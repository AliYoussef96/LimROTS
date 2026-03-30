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
#' (columns of `data.exp`). If `data.exp` is a `SummarizedExperiment`,
#' `meta.info` can be a vector of `colData` column names to use.
#' @param group.name Character. Column name in `meta.info` that defines the
#' groups or conditions for comparison.
#' @param formula.str Optional character string representing the formula for
#' the model.
#' @param verbose Logical, indicating whether to display messages during the
#' function's execution. Default is TRUE.
#' @param log Logical, indicating whether the data is already log-transformed.
#' Default is TRUE.
#' @param survival Logical, indicating whether the analysis is survival analysis.
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
        meta.info, group.name,
        formula.str, verbose = TRUE,
        log = TRUE, survival = FALSE) {
    data.exp <- x
    Check_SExp <- Check_SummarizedExperiment(
        data.exp = x,
        assay.type = assay.type,
        meta.info = meta.info,
        group.name = group.name,
        survival = survival
    )
    data <- Check_SExp$data
    groups <- Check_SExp$groups
    meta.info <- Check_SExp$meta.info

    Check_meta_info(
        meta.info = meta.info,
        data = data, log = log
    )
    
    if(survival){
      
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
        if (verbose) {
          message(sprintf("No top list size K given, using %s", K))
        }
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
      
      
      sort.df <- data.frame(
        sample.id = colnames(data),
        event = meta.info[, "event"]
      )
      sort.df <- sort.df[order(sort.df$event), ]
      data <- data[, sort.df$sample.id]
      meta.info$temp <- row.names(meta.info)
      meta.info <- data.frame(meta.info[colnames(data), ],
                              check.rows = FALSE,
                              check.names = FALSE
      )
      meta.info$temp <- NULL
      
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

    sort.df <- data.frame(
        sample.id = colnames(data),
        groups = meta.info[, group.name]
    )
    sort.df <- sort.df[order(sort.df$groups), ]
    data <- data[, sort.df$sample.id]
    meta.info$temp <- row.names(meta.info)
    meta.info <- data.frame(
        meta.info[colnames(data), ],
        check.rows = FALSE,
        check.names = FALSE
    )
    meta.info$temp <- NULL
    if (inherits(meta.info[, group.name], "character")) {
        meta.info[, group.name] <- as.factor(
            meta.info[, group.name]
        )
        message(paste(
            "Group Level: ",
            levels(meta.info[, group.name]),
            collapse = " & "
        ))
        meta.info[, group.name] <- as.numeric(
            meta.info[, group.name]
        )
        groups <- as.numeric(meta.info[, group.name])
        meta.info[, group.name] <- as.factor(
            meta.info[, group.name]
        )
    } else if (inherits(
        meta.info[, group.name], "factor"
    )) {
        groups <- as.numeric(meta.info[, group.name])
        message(paste(
            "Group Level: ",
            levels(meta.info[, group.name]),
            collapse = " & "
        ))
        meta.info[, group.name] <- as.numeric(
            meta.info[, group.name]
        )
        meta.info[, group.name] <- as.factor(
            meta.info[, group.name]
        )
    } else {
        meta.info[, group.name] <- as.factor(
            meta.info[, group.name]
        )
        message(paste(
            "Group Level: ",
            levels(meta.info[, group.name]),
            collapse = " & "
        ))
        meta.info[, group.name] <- as.numeric(
            meta.info[, group.name]
        )
        groups <- as.numeric(meta.info[, group.name])
        meta.info[, group.name] <- as.factor(
            meta.info[, group.name]
        )
    }
    groups <- groups + (1 - min(groups))
    if (is.null(K)) {
        K <- floor(nrow(data) / 4)
        if (verbose) {
            message(sprintf(
                "No top list size K given, using %s", K
            ))
        }
    }
    return(list(
        meta.info = meta.info,
        data = data,
        groups = groups,
        K = K
    ))
}