#' Sanity Check for Input Data and Parameters
#'
#' This function performs a series of checks and initial setups for input data, metadata, and parameters, ensuring everything is correctly formatted for downstream analysis.
#'
#' @param data.exp A matrix-like object or a `SummarizedExperiment` containing the data to be analyzed.
#' @param B Integer. Number of bootstrap samples or resampling iterations. Default is 1000.
#' @param K Integer. Top list size. If NULL, it will be set to a quarter of the number of rows in the data matrix. Default is NULL.
#' @param a1,a2 Optional numeric parameters related to optimization.
#' @param meta.info Data frame. Metadata associated with the samples (columns of `data.exp`). If `data.exp` is a `SummarizedExperiment`, `meta.info` can be a vector of `colData` column names to use.
#' @param group.name Character. Column name in `meta.info` that defines the groups or conditions for comparison.
#' @param formula.str Optional character string representing the formula for the model.
#' @param survival Logical. If TRUE, survival analysis is used, requiring `time` and `event` columns in `meta.info`. Default is FALSE.
#' @param paired Logical. If TRUE, indicates paired test setup. Default is FALSE.
#' @param n.ROTS Logical. If TRUE, uses the ROTS method instead of LimROTS. Default is FALSE.
#' @param verbose Logical, indicating whether to display messages during the function's execution. Default is TRUE.
#' @param log Logical, indicating whether the data is already log-transformed. Default is TRUE.
#'
#' @details
#' This function checks whether the input data and metadata are in the correct format, processes metadata from a `SummarizedExperiment` object if provided, and ensures that group information is correctly specified. If no top list size (`K`) is provided, it defaults to a quarter of the number of rows in the data.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{meta.info}: Processed metadata.
#'   \item \code{data}: Processed data matrix.
#'   \item \code{groups}: Numeric or factor vector indicating group assignments.
#'   \item \code{event}: Event data for survival analysis (if applicable).
#'   \item \code{K}: Top list size to be used in the analysis.
#' }
#'
#'



SanityChecK <- function(data.exp,
                        B = 1000,
                        K = NULL,
                        a1 = NULL,
                        a2 = NULL,
                        meta.info = NULL,
                        group.name = NULL ,
                        formula.str = NULL,
                        survival = FALSE,
                        paired = FALSE,
                        n.ROTS = FALSE,
                        verbose = TRUE,
                        log = TRUE) {
    if (survival == TRUE) {
        if (n.ROTS == FALSE) {
            stop(
                "survival is TRUE, Survival analysis is only available through the old ROTS implementation. To enable this, please set n.ROTS to TRUE."
            )
        }
    }

    if (paired == TRUE) {
        if (n.ROTS == FALSE) {
            stop(
                "paired if TRUE, Paired analysis is only available through the old ROTS implementation. To enable this, please set n.ROTS to TRUE."
            )
        }
    }

    if (log == FALSE) {
        if (n.ROTS == FALSE) {
            stop(
                "log is FALSE, Unlogged values can only be handled through the old ROTS implementation. To enable this, please set n.ROTS to TRUE."
            )

        }
    }

    ### SummarizedExperiment
    if (inherits(data.exp, "SummarizedExperiment")) {
        message("Data is SummarizedExperiment object")

        if (is.null(meta.info)) {
            stop("meta.info should be a vector of colData names to be used")
        } else{
            meta.info.colnames <- meta.info
            meta.info <- data.frame(
                colData(data.exp)[, meta.info],
                check.names = FALSE,
                row.names = row.names(colData(data.exp))
            )
            if (length(meta.info) != length(meta.info.colnames)) {
                stop("meta.info should be a vector of colData names to be used")
            } else{
                colnames(meta.info) <- meta.info.colnames
            }

        }
        if (!group.name %in% colnames(meta.info)) {
            stop(
                "group.name should be a string specifying the column in `meta.info` that represents the groups or conditions for comparison."
            )
        }
        message(sprintf("Assay: %s will be used" , assayNames(data.exp)[1]))
        data <- assay(data.exp , assayNames(data.exp)[1])
    } else{
        data <- data.exp
        groups <- meta.info[, group.name]
    }
    ### meta.info
    if (any(!row.names(meta.info) %in% colnames(data))) {
        stop("rownames for meta.info should match the data colnames (samples names)")
    }
    if (any(grepl("." , colnames(data) , fixed = TRUE))) {
        stop("Sample names should contains no '.', please remove it if any")
    }
    if (!is.null(meta.info) & n.ROTS == FALSE) {
        if (ncol(meta.info) == 1) {
            message(
                "A meta.info table is provided with only group infomration >>> LimROTS with no covariates will be used"
            )
            if (is.null(formula.str)) {
                stop("formula.str should by provided for the model")
            }
        } else{
            message(
                "A meta.info table is provided with covariates >>> LimROTS with covariates will be used"
            )
            if (is.null(formula.str)) {
                stop("formula.str should by provided for the model")
            }
        }
    } else{
        message("n.ROTS is TRUE >>> ROTS will be used")
    }

    ### Sort
    if (nrow(meta.info) != ncol(data)) {
        stop("Number of samples in the data does not match the groups.")
    }
    sort.df <- data.frame(sample.id = colnames(data), groups = meta.info[, group.name])
    sort.df <- sort.df[order(sort.df$groups), ]
    data <- data[, sort.df$sample.id]
    meta.info$temp <- row.names(meta.info)
    meta.info <- data.frame(meta.info[colnames(data), ],
                            check.rows = FALSE,
                            check.names = FALSE)
    meta.info$temp <- NULL
    ### Groups
    if (inherits(meta.info[, group.name], "character")) {
        meta.info[, group.name] <- factor(meta.info[, group.name])
        groups.levels <- levels(groups)
        meta.info[, group.name] <- as.numeric(meta.info[, group.name])
    } else if (inherits(meta.info[, group.name], "factor")) {
        groups <- as.numeric(meta.info[, group.name])
    } else{
        meta.info[, group.name] <- as.numeric(meta.info[, group.name])
        groups <- as.numeric(meta.info[, group.name])
    }
    if (survival  == TRUE) {
        if (all(c("time", "event") %in% colnames(meta.info))) {
            stop(
                "meta.info must have two columns time and event. Also, group.name must be time"
            )
        }
        event <- meta.info[, "event"]
        groups <- meta.info[, "time"]
    } else{
        event <- NULL
        time <- NULL
    }
    groups <- groups + (1 - min(groups))
    ### paired
    if (paired) {
        for (i in unique(groups)[-1]) {
            if (length(which(groups == 1)) != length(which(groups ==
                                                           i)))
                stop("Uneven number of samples for paired test.")
        }
    }
    if (is.null(K)) {
        K <- floor(nrow(data) / 4)
        if (verbose)
            message(sprintf("No top list size K given, using %s", K))
    }
    return(list(
        meta.info = meta.info,
        data = data,
        groups = groups,
        event = event,
        K = K
    ))
}
