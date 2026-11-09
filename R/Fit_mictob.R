parallel_microb <- function(BPPARAM = NULL,
                          samples,
                          data,
                          formula.str,
                          group,
                          meta.info,
                          pSamples,
                          pairwise_contrasts,
                          trend,
                          robust) {

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
            pairwise_contrasts = pairwise_contrasts,
            trend = trend,
            robust = robust
        )
        d_result <- fit$d
        s_result <- fit$s
        df1 <- data.frame(d_result = d_result, s_result = s_result)
        pFit <- fit_permutation_microb(
            x = data,
            group = group,
            meta.info = pSamples_i,
            formula.str = formula.str,
            pairwise_contrasts = pairwise_contrasts,
            trend = trend,
            robust = robust
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

    ### Fit the contrasts 
    fit2 <- contrasts.fit(fit, cont_matrix)
    fit.ebayes <- eBayes(fit2, trend = trend, robust = robust)

    beta <- as.numeric(fit.ebayes[["coefficients"]])
    bias <- .estimate_bias(
        beta,
        ncol(processed_counts)
    )
    coef <- beta - bias

    d_values <- setNames(abs(coef), rownames(fit.ebayes[["coefficients"]]))
    corrected.logfc <- setNames(coef , rownames(fit.ebayes[["coefficients"]]))

    stopifnot(identical(names(fit.ebayes[["s2.post"]]) ,  
                    row.names(fit.ebayes[["stdev.unscaled"]])))

    s_values <- setNames(
        as.numeric(
            sqrt(fit.ebayes[["s2.post"]]) *
                        fit.ebayes[["stdev.unscaled"]]
                        ),
        row.names(fit.ebayes[["stdev.unscaled"]]))

    stopifnot(identical(names(s_values), names(d_values) ) )
    stopifnot(identical(names(corrected.logfc), names(d_values) ) )
    stopifnot(identical(names(s_values), names(corrected.logfc) ) )

    return(list(d = d_values, s = s_values, corrected.logfc = corrected.logfc))

}

fit_bootstrap_microb <- function(x, group, meta.info, formula.str,
 trend, robust, pairwise_contrasts){
                            
    ### Generate combined data and covariates for model fitting
    data_for_fit <- .generate_data_for_fit(x = x, meta.info = meta.info)
    combined_data <- data_for_fit$combined_data
    covariates.p <- data_for_fit$covariates.p

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

    ### Fit the contrasts 
    fit2 <- contrasts.fit(fit, cont_matrix)
    fit.ebayes <- eBayes(fit2, trend = trend, robust = robust)

    beta <- as.numeric(fit.ebayes[["coefficients"]])
    bias <- .estimate_bias(
        beta,
        ncol(processed_counts)
    )
    coef <- beta - bias

    d_values <- setNames(abs(coef), rownames(fit.ebayes[["coefficients"]]))

    stopifnot(identical(names(fit.ebayes[["s2.post"]]) ,  
                    row.names(fit.ebayes[["stdev.unscaled"]])))

    s_values <- setNames(
        as.numeric(
            sqrt(fit.ebayes[["s2.post"]]) *
                        fit.ebayes[["stdev.unscaled"]]
                        ),
        row.names(fit.ebayes[["stdev.unscaled"]]))

    stopifnot(identical(names(s_values), names(d_values) ) )

    return(list(d = d_values, s = s_values))
}

fit_permutation_microb <- function(x, group, meta.info, formula.str,
 trend, robust, pairwise_contrasts){
    
    ### Generate combined data and covariates for model fitting
    combined_data <- x
    covariates.p <- meta.info
    row.names(covariates.p) <- NULL
    
    ### Create the design matrix based on the provided formula and covariates
    design.matrix <-
        model.matrix(formula(formula.str), data = covariates.p)
    colnames(design.matrix) <-
        make.names(colnames(design.matrix))
 

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

    ### Fit the contrasts 
    fit2 <- contrasts.fit(fit, cont_matrix)
    fit.ebayes <- eBayes(fit2, trend = trend, robust = robust)

    beta <- as.numeric(fit.ebayes[["coefficients"]])
    bias <- .estimate_bias(
        beta,
        ncol(processed_counts)
    )
    coef <- beta - bias

    d_values <- setNames(abs(coef), rownames(fit.ebayes[["coefficients"]]))

    stopifnot(identical(names(fit.ebayes[["s2.post"]]) ,  
                    row.names(fit.ebayes[["stdev.unscaled"]])))

    s_values <- setNames(
        as.numeric(
            sqrt(fit.ebayes[["s2.post"]]) *
                        fit.ebayes[["stdev.unscaled"]]
                        ),
        row.names(fit.ebayes[["stdev.unscaled"]]))

    stopifnot(identical(names(s_values), names(d_values) ) )

    return(list(d = d_values, s = s_values))

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


.generate_data_for_fit <- function(x, meta.info, correlation_block) {
    ### Check if x is a list of matrices with identical rownames
    stopifnot(
        all(vapply(
            x,
            function(z) identical(rownames(z), rownames(x[[1L]])),
            logical(1)
        ))
    )
    combined_data <- do.call(cbind, x)
    
    ### create a unique colnames for the combined data to avoid duplicates
    colnames(combined_data) <-
        paste0(colnames(combined_data), ".", seq(1, ncol(combined_data)))

    ### create a covariates data frame for the combined data
    covariates.p <- data.frame()
    meta.info.temp <- meta.info
    meta.info.temp$sample.id <- row.names(meta.info.temp)
    for (i in colnames(combined_data)) {
        real_SampleNames <- str_split_fixed(i, fixed("."), 2)[, 1]
        df.temp <- meta.info.temp[real_SampleNames, , drop = FALSE]
        df.temp[,"sample.id"] <- i
        covariates.p <- rbind(covariates.p, df.temp)
    }
    stopifnot(identical(colnames(combined_data), covariates.p[,"sample.id"]))
    row.names(covariates.p) <- covariates.p[,"sample.id"]
    covariates.p[,"sample.id"] <- NULL

    stopifnot(identical(colnames(combined_data), row.names(covariates.p)))
    return(list(combined_data = combined_data, covariates.p = covariates.p))
}