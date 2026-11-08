parallel_microb <- function(BPPARAM = NULL,
                          samples,
                          data,
                          formula.str,
                          group,
                          meta.info,
                          pSamples,
                          pairwise_contrasts){

    if (is.null(BPPARAM)) {
        if (.Platform$OS.type == "windows") {
            BPPARAM <- SnowParam(workers = 2)
            message("Using SnowParam (Windows) with two workers.")
        } else {
            BPPARAM <- MulticoreParam(workers = 2)
            message("Using MulticoreParam (Unix-like OS) with two workers.")
        }
    } else {
        message("Using provided parallel backend BPPARAM.")
    }
   
    results_list <- bplapply(seq_len(nrow(samples)), function(i) {
        samples.R <- split(samples[i, ], meta.info[samples[i, ], group])
        pSamples_i <- pSamples[[i]]
        d_result <- NULL
        s_result <- NULL
        pd_result <- NULL
        ps_result <- NULL
        fit <- fit_bootstrap_microb(
            x = lapply(samples.R, function(x)
                data[, x]
            ),
            group = group,
            meta.info = meta.info,
            formula.str = formula.str,
            pairwise_contrasts = pairwise_contrasts
        )
        d_result <- fit$d
        s_result <- fit$s
        df1 <- data.frame(d_result = d_result, s_result = s_result)
        pFit <- fit_permutation_microb(
            x = data,
            group = group,
            meta.info = pSamples_i,
            formula.str = formula.str,
            pairwise_contrasts = pairwise_contrasts
        )
        pd_result <- pFit$d
        ps_result <- pFit$s
        df2 <- data.frame(pd_result = pd_result, ps_result = ps_result)
        list(ds = df1, pdps = df2)
    }, BPPARAM = BPPARAM, 
    BPOPTIONS = bpoptions(packages = c("utils", 
                                       "stringr", "stats", 
                                       "limma")))
    return(results_list)
}


fit_original_microb <- function(x, group, meta.info, formula.str, 
            trend, robust, pairwise_contrasts){
    combined_data <- x
    if(!identical(colnames(combined_data), row.names(meta.info))){
        stop("Column names of the data matrix must match row names of meta.info.")
    }
    ### Create the design matrix based on the provided formula and covariates
    design.matrix <- model.matrix(formula(formula.str), data = meta.info)
    colnames(design.matrix) <- make.names(colnames(design.matrix))
    stopifnot(identical(colnames(combined_data) , row.names(design.matrix)))

    ### Fit the linear model using Limma, accounting for correlation blocks if specified
    fit <- lmFit(combined_data, design.matrix)
    cont_matrix <- makeContrasts(
            contrasts = pairwise_contrasts, levels =
                design.matrix
        )

    ### Fit the contrasts and apply empirical Bayes moderation
    fit2 <- limma::contrasts.fit(fit, cont_matrix)


    beta <- fit$coefficients[, 1]
    bias <- .estimate_bias(
        beta,
        ncol(processed_counts)
    )


}


.estimate_bias <- function(coefficient, n_samples) {
    tryCatch(
        as.numeric(
            modeest::mlv(
                sqrt(n_samples) * coefficient,
                method = "meanshift",
                kernel = "gaussian"
            )
        ) / sqrt(n_samples),
        error = function(e) median(coefficient, na.rm = TRUE)
    )
}