# Copyright (C) 2025 Ali Youssef.
# This program is licensed under GPL (>= 2).
#
#' Check if SummarizedExperiment or data is correct
#'
#' @param data.exp A matrix-like object or a `SummarizedExperiment` containing
#' the data to be analyzed.
#' @param assay.type A character string or numeric index specifying the assay
#' to use if `data.exp` is a `SummarizedExperiment`. Default is `NULL`
#' @param meta.info Data frame. Metadata associated with the samples
#' (columns of `data.exp`). If `data.exp` is a `SummarizedExperiment`,
#' @param group.name Character. Column name in `meta.info` that defines the
#' groups or conditions for comparison.
#'
#' @import SummarizedExperiment
#'
#' @return a list of `data` , `groups` and `meta.info`
#'
#'
#'

Check_SummarizedExperiment <- function(data.exp,
                                       assay.type = NULL, meta.info, group.name) {
  if (inherits(data.exp, "SummarizedExperiment")) {
    message("Data is SummarizedExperiment object")
    
    if (is.null(meta.info)) {
      stop("meta.info should be a vector of colData names to be used")
    } else {
      meta.info.colnames <- meta.info
      meta.info <- data.frame(
        colData(data.exp)[, meta.info],
        check.names = FALSE,
        row.names = row.names(colData(data.exp))
      )
      if (length(meta.info) != length(meta.info.colnames)) {
        stop("meta.info should be a vector of colData names to be used")
      } else {
        colnames(meta.info) <- meta.info.colnames
      }
    }
    if (!group.name %in% colnames(meta.info)) {
      stop(
        "group.name should be a string specifying the column in
                `meta.info` that represents the groups or conditions
                for comparison."
      )
    }
    if (is.null(assay.type)) {
      assay.type <- assayNames(data.exp)[1]
    }
    message(sprintf("Assay: %s will be used", assay.type))
    data <- assay(data.exp, assay.type)
    groups <- NULL
  } else {
    data <- data.exp
    groups <- meta.info[, group.name]
    meta.info <- meta.info
  }
  return(list(data = data, groups = groups, meta.info = meta.info))
}