#' Check if SummarizedExperiment or data is correct
#'
#' @param meta.info Data frame. Metadata associated with the samples
#' (columns of `data.exp`). If `data.exp` is a `SummarizedExperiment`,
#' @param data.exp A matrix-like object or a `SummarizedExperiment` containing
#' the data to be analyzed.
#' @param group.name Character. Column name in `meta.info` that defines the
#' groups or conditions for comparison.
#'

Check_SummarizedExperiment <- function(data.exp, meta.info, group.name) {
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
        message(sprintf("Assay: %s will be used", assayNames(data.exp)[1]))
        data <- assay(data.exp, assayNames(data.exp)[1])
    } else {
        data <- data.exp
        groups <- meta.info[, group.name]
    }
    return(list(data = data, groups = groups))
}
