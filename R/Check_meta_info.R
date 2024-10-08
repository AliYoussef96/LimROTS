#'Check if meta info is correct
#'
#' @param meta.info Data frame. Metadata associated with the samples
#' (columns of `data.exp`). If `data.exp` is a `SummarizedExperiment`,
#' @param data A matrix-like object or a `SummarizedExperiment` containing the
#' data to be analyzed.
#' @param log Logical, indicating whether the data is already log-transformed.
#' Default is TRUE.
#'
#' @return Logical
#'
#'

Check_meta_info <- function(meta.info, data, log) {
    if (any(!row.names(meta.info) %in% colnames(data))) {
        stop("rownames
            for meta.info should match the data colnames (samples names)")
    }
    if (any(grepl(".", colnames(data), fixed = TRUE))) {
        stop("Sample names should contains no '.', please remove it if any")
    }
    if (nrow(meta.info) != ncol(data)) {
        stop("Number of samples in the data does not match the groups.")
    }
    if (log == FALSE) {
        stop(
            "log is FALSE, Unlogged values, Please log2 the data and
                set log to TRUE.
                If the data is log2 values then only set log to TRUE."
        )
    }

    return(TRUE)
}
