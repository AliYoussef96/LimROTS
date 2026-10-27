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

Check_meta_info <- function(meta.info, data) {
    if (any(!row.names(meta.info) %in% colnames(data))) {
        stop("rownames
            for colData should match the data colnames (samples names)")
    }
    if (any(grepl(".", colnames(data), fixed = TRUE))) {
        stop("Sample names should contains no '.', please remove it if any")
    }
    if (nrow(meta.info) != ncol(data)) {
        stop("Number of samples in the data does not match the groups")
    }
    if(!identical(sort(rownames(meta.info)), sort(colnames(data)))){
        stop("Sample names in the data and meta info in colData do not match")
    }

    return(TRUE)
}
