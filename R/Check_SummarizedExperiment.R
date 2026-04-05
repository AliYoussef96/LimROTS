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
#' groups or conditions for comparison.
#' @param survival Logical, indicating whether the analysis is survival
#' analysis. Default is \code{FALSE}.
#'
#' @import SummarizedExperiment
#'
#' @return a list of `data` , `groups` and `meta.info`
#'
#'
#'

Check_SummarizedExperiment <- function(x,
        assay.type = NULL, meta.info,
        group, survival = FALSE) {
    if (inherits(x, "SummarizedExperiment")) {
        message("Data is SummarizedExperiment object")

        if (is.null(meta.info)) {
            stop("meta.info should be a vector of colData names to be used")
        } else {
            meta.info.colnames <- meta.info
            meta.info <- data.frame(
                colData(x)[, meta.info],
                check.names = FALSE,
                row.names = row.names(colData(x))
            )
            if (length(meta.info) != length(meta.info.colnames)) {
                stop("meta.info should be a vector of colData names to be used")
            } else {
                colnames(meta.info) <- meta.info.colnames
            }
        }
        if(!survival){
          if (!group %in% colnames(meta.info)) {
              stop(
                  "group should be a string specifying",
                  " the column in `meta.info` that",
                  " represents the groups or conditions",
                  " for comparison."
              )
          }
        }
        if (is.null(assay.type)) {
            assay.type <- assayNames(x)[1]
        }
          message(sprintf("Assay: %s will be used", assay.type))
          data <- assay(x, assay.type)
          groups <- NULL
      } else {
          data <- x
          if (!survival) {
              groups <- meta.info[, group]
          } else {
              groups <- NULL
          }
          meta.info <- meta.info
      }
      return(list(
          data = data, groups = groups,
          meta.info = meta.info
      ))
  }