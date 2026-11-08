LimROTS_microb <- function(
    assay.type = NULL,
    niter = 1000,
    K = NULL,
    meta.info,
    BPPARAM  = NULL,
    group,
    formula.str,
    robust = TRUE,
    trend = TRUE,
    correlation_block = NULL,
    pairwise_contrasts,
    model,
    strata){
    ### Perform sanity checks on the input parameters and data.
    SanityChecK.list <- SanityChecK(
        x,
        assay.type = assay.type,
        niter = niter,
        K = K,
        meta.info = meta.info,
        group = group,
        model = model,
        formula.str = formula.str,
        pairwise_contrasts = pairwise_contrasts,
        correlation_block = correlation_block,
        strata = strata,
        survival = FALSE,
        competing_risks = FALSE
    )

    meta.info <- SanityChecK.list[["meta.info"]]
    data <- SanityChecK.list[["data"]]
    group <- SanityChecK.list[["group"]]
    K <- SanityChecK.list[["K"]]

    if (anyDuplicated(rownames(feature_data)) ||
        anyDuplicated(colnames(feature_data)) ||
        anyDuplicated(rownames(metadata))) {
        stop("Feature names and sample names must be unique.")
    }

    if (anyNA(feature_data) || any(!is.finite(feature_data))) {
        stop("feature_data contains missing or non-finite values.")
    }

    if (any(feature_data < 0)) {
        stop("feature_data cannot contain negative values.")
    }

    if (any(feature_data != round(feature_data))) {
        stop("feature_data must contain original integer counts.")
    }

    if (!all(colnames(feature_data) %in% rownames(metadata))) {
        stop(
            "Samples missing from metadata: ",
            paste(
                setdiff(colnames(feature_data), rownames(metadata)),
                collapse = ", "
            )
        )
    }



}