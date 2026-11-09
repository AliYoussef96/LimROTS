# Copyright (C) 2025 Ali Youssef.
# This program is licensed under GPL (>= 2).
#
#' Fit Bootstrap and Permutation Replicates in Parallel
#'
#' Apply the selected linear-model fitting method to each bootstrap replicate
#' and its corresponding permuted metadata in parallel.
#'
#' @param BPPARAM A `BiocParallelParam` object specifying the parallel backend.
#' If `NULL`, the function selects a two-worker backend appropriate to the
#' operating system.
#' @param samples Matrix of bootstrap replicates, with sample names in rows and
#' resampled sample names in columns.
#' @param data Feature-by-sample matrix with sample names as column names.
#' @param formula.str Character string containing the model formula.
#' @param group Character string naming the column in `meta.info` that
#' identifies each sample's group.
#' @param meta.info Data frame of sample-level metadata, with sample names as
#' row names.
#' @param pSamples List of permuted sample-level metadata data frames, one per
#' bootstrap replicate.
#' @param correlation_block Character string or `NULL`. Name of the column in
#' `meta.info` identifying correlated samples or blocks.
#' @param model Character string specifying the fitting method: `"limma"` for
#' fixed effects or `"dream"` for random effects.
#' @param pairwise_contrasts Character vector of contrast expressions, such as
#' `"group1 - group2"`, to fit for each replicate.
#'
#' @return A list with one element per replicate. Each element contains `ds`
#' (bootstrap statistics and standard errors) and `pdps` (permutation
#' statistics and standard errors).
#'
#' @importFrom BiocParallel SnowParam MulticoreParam bplapply
#' @importFrom stats model.matrix formula p.adjust
#'


Boot_parallel <- function(BPPARAM = NULL,
                          samples,
                          data,
                          formula.str,
                          group,
                          meta.info,
                          pSamples,
                          correlation_block,
                          model,
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
    if (inherits(BPPARAM, "SnowParam")) {
        BPPARAM$exportglobals <- FALSE
    }
    export_vars <- list(
        samples = samples,
        data = data,
        formula.str = formula.str,
        group = group,
        meta.info = meta.info,
        pSamples = pSamples,
        correlation_block = correlation_block,
        model = model,
        pairwise_contrasts = pairwise_contrasts,
        trend = trend,
        robust = robust
    )

    for (name in names(export_vars)) {
        assign(name, export_vars[[name]])
    }

    results_list <- bplapply(seq_len(nrow(samples)), function(i) {
        
        samples.R <- split(samples[i, ], meta.info[samples[i, ], group])
        pSamples_i <- pSamples[[i]]
        d_result <- NULL
        s_result <- NULL
        pd_result <- NULL
        ps_result <- NULL
        fit <- fit_bootstrap(
            x = lapply(samples.R, function(x)
                data[, x]
            ),
            group = group,
            meta.info = meta.info,
            formula.str = formula.str,
            correlation_block = correlation_block,
            model = model,
            pairwise_contrasts = pairwise_contrasts,
            trend = trend,
            robust = robust
        )
        d_result <- fit$d
        s_result <- fit$s
        df1 <- data.frame(d_result = d_result, s_result = s_result)
        pFit <- fit_permutation(
            x = data,
            group = group,
            meta.info = pSamples_i,
            formula.str = formula.str,
            correlation_block = correlation_block,
            model = model,
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
                                       "limma" , "variancePartition")))
    return(results_list)
}


#' Fit the Original Data with Limma or DREAM
#'
#' Fit the selected linear model to the original feature-by-sample matrix and
#' return feature statistics and standard errors.
#'
#' @param x Feature-by-sample matrix with sample names as column names.
#' @param group Character string naming the column in `meta.info` that
#' identifies each sample's group.
#' @param meta.info Data frame of sample-level metadata, with sample names as
#' row names.
#' @param formula.str Character string containing the model formula.
#' @param trend Logical indicating whether to fit an intensity-dependent trend
#' in the prior variance.
#' @param robust Logical indicating whether to use robust empirical Bayes
#' moderation.
#' @param correlation_block Character string or `NULL`. Name of the column in
#' `meta.info` identifying correlated samples or blocks.
#' @param model Character string specifying the fitting method: `"limma"` for
#' fixed effects or `"dream"` for random effects.
#' @param pairwise_contrasts Character vector of contrast expressions, such as
#' `"group1 - group2"`, to fit.
#'
#' @details
#' For `"limma"`, fits a linear model and applies empirical Bayes moderation.
#' For `"dream"`, fits a linear mixed model and applies empirical Bayes moderation. 
#' The specified contrasts areapplied to the original data.
#'
#' @return A list containing the following elements:
#' \item{d}{Absolute model statistics for each feature.}
#' \item{s}{Model standard errors for each feature after empirical Bayes
#' moderation.}
#' \item{corrected.logfc}{Signed model coefficients for each feature.}
#'
#' @seealso \code{\link[limma]{lmFit}},
#' \code{\link[limma]{eBayes}},
#' \code{\link[limma]{topTable}},
#' \code{\link[limma]{makeContrasts}},
#' \code{\link[variancePartition]{dream}},
#' \code{\link[variancePartition]{makeContrastsDream}}
#' @importFrom stats model.matrix formula
#' @importFrom dplyr bind_cols
#' @importFrom utils combn
#' @importFrom limma makeContrasts lmFit contrasts.fit
#' @importFrom limma duplicateCorrelation
#' @importFrom variancePartition dream makeContrastsDream
#' @importFrom BiocParallel SerialParam
#'


fit_original <- function(x, group, meta.info, formula.str, trend, robust,
                          correlation_block, model,
                          pairwise_contrasts) {
    if(model == "limma"){
        return(.Limma_fit_original(x = x, group = group, 
        meta.info = meta.info, formula.str = formula.str, 
        correlation_block = correlation_block, 
        pairwise_contrasts = pairwise_contrasts,
        trend = trend, robust = robust ))

    }else if (model == "dream") {
        return(.Dream_fit_original(x = x, group = group, 
        meta.info = meta.info, formula.str = formula.str, 
        correlation_block = correlation_block, 
        pairwise_contrasts = pairwise_contrasts,
        trend = trend, robust = robust ))
    }else{
        stop("Invalid model specified. Choose either 'limma' or 'dream'.")
    }
}

#' Fit Original Data with Limma
#'
#' Construct the design matrix, fit the original data with Limma, apply the
#' requested contrasts and empirical Bayes moderation, and return feature-level
#' statistics.
#'
#' @param x Feature-by-sample matrix with sample names as column names.
#' @param group Character string naming the column in `meta.info` that
#' identifies each sample's group.
#' @param meta.info Data frame of sample-level metadata, with sample names as
#' row names.
#' @param formula.str Character string containing the model formula.
#' @param trend Logical indicating whether to fit an intensity-dependent trend
#' in the prior variance.
#' @param robust Logical indicating whether to use robust empirical Bayes
#' moderation.
#' @param correlation_block Character string or `NULL`. Name of the column in
#' `meta.info` identifying correlated samples or blocks.
#' @param pairwise_contrasts Character vector of contrast expressions, such as
#' `"group1 - group2"`, to fit.
#'
#' @return A list containing absolute feature statistics `d`, standard errors
#' `s`, and signed model coefficients `corrected.logfc`.
.Limma_fit_original <- function(x, group, meta.info, formula.str, 
            trend, robust, correlation_block, pairwise_contrasts) {
        ### Check if x is a list of matrices with identical rownames
        combined_data <- x

        if(!identical(colnames(combined_data), row.names(meta.info))){
            stop("Column names of the data matrix must match row names of meta.info.")
        }

        ### Create the design matrix based on the provided formula and covariates
        design.matrix <- model.matrix(formula(formula.str), data = meta.info)
        colnames(design.matrix) <- make.names(colnames(design.matrix))
        stopifnot(identical(colnames(combined_data) , row.names(design.matrix)))

        ### Fit the linear model using Limma, accounting for correlation blocks if specified
        if (!is.null(correlation_block)) {
            corfit <- duplicateCorrelation(combined_data, design.matrix,
                                        block = meta.info[, correlation_block])
            fit <- lmFit(combined_data, design.matrix,
                        block = meta.info[, correlation_block],
                        correlation = corfit$consensus.correlation)
        } else {
            fit <- lmFit(combined_data, design.matrix)
        }

        ### For two-group comparisons, compute the specified pairwise contrasts
        cont_matrix <- makeContrasts(
            contrasts = pairwise_contrasts, levels =
                design.matrix
        )

        ### Fit the contrasts and apply empirical Bayes moderation
        fit2 <- limma::contrasts.fit(fit, cont_matrix)
        fit.ebayes <- limma::eBayes(fit2, trend = trend, robust = robust)

        d_values <- setNames(
            abs(as.numeric(fit.ebayes[["coefficients"]])),
            rownames(fit.ebayes[["coefficients"]]))

        corrected.logfc <- setNames(
            as.numeric(fit.ebayes[["coefficients"]]),
            rownames(fit.ebayes[["coefficients"]]))

        if (is.null(correlation_block)) {
        stopifnot(identical(names(fit.ebayes[["s2.post"]]) ,  
                        row.names(fit.ebayes[["stdev.unscaled"]])))
        }

        s_values <- setNames(
            as.numeric(sqrt(fit.ebayes[["s2.post"]]) *
                            fit.ebayes[["stdev.unscaled"]]),
            row.names(fit.ebayes[["stdev.unscaled"]]))

        stopifnot(identical(names(s_values), names(d_values) ) )
        stopifnot(identical(names(corrected.logfc), names(d_values) ) )
        stopifnot(identical(names(s_values), names(corrected.logfc) ) )

        return(list(
            d = d_values,
            s = s_values,
            corrected.logfc = corrected.logfc
        ))
            
    }


#' Fit Original Data with DREAM
#'
#' Fit a linear mixed model to the original data with DREAM, apply the requested
#' contrasts and empirical Bayes moderation, and return feature-level statistics.
#'
#' @param x Feature-by-sample matrix with sample names as column names.
#' @param group Character string naming the column in `meta.info` that
#' identifies each sample's group.
#' @param meta.info Data frame of sample-level metadata, with sample names as
#' row names.
#' @param formula.str Character string containing the model formula.
#' @param trend Logical indicating whether to fit an intensity-dependent trend
#' in the prior variance.
#' @param robust Logical indicating whether to use robust empirical Bayes
#' moderation.
#' @param correlation_block Character string or `NULL`. Name of the column in
#' `meta.info` identifying correlated samples or blocks.
#' @param pairwise_contrasts Character vector of contrast expressions, such as
#' `"group1 - group2"`, to fit.
#'
#' @return A list containing absolute feature statistics `d`, standard errors
#' `s`, and signed model coefficients `corrected.logfc`.
.Dream_fit_original <- function(x, group, meta.info, formula.str, trend, robust,
            correlation_block, pairwise_contrasts) {

        combined_data <- x

        if(!identical(colnames(combined_data), row.names(meta.info))){
            stop("Column names of the data matrix must match row names of meta.info.")
        }

        ### Create the design matrix based on the provided formula and covariates
        L <- variancePartition::makeContrastsDream(
            formula(formula.str), meta.info,
            contrasts = pairwise_contrasts
        )

        ### Fit the model using DREAM
        fit <- variancePartition::dream(
            combined_data, formula(formula.str), meta.info,
            L = L,
            BPPARAM = BiocParallel::SerialParam()
        )
        fit.ebayes <- variancePartition::eBayes(
            fit, trend = trend, robust = robust
        )
        d_values <- variancePartition::topTable(
            fit.ebayes,
            coef = pairwise_contrasts,
            number = "Inf",
            sort.by = "none"
        )
        corrected.logfc <- d_values[["logFC"]]
        d_values <- abs(d_values[["logFC"]])
        s_values <-
            as.numeric(sqrt(fit.ebayes[["s2.post"]]) *
                fit.ebayes[["stdev.unscaled"]][, 1])
        return(list(
            d = d_values,
            s = s_values,
            corrected.logfc = corrected.logfc
        ))
}


#' Fit Bootstrap Samples with Limma or DREAM
#'
#' Fit the selected linear model to a bootstrap replicate assembled from
#' group-specific feature-by-sample matrices.
#'
#' @param x List of feature-by-sample matrices, one per group in the bootstrap
#' replicate. Matrices must have matching feature names.
#' @param group Character string naming the column in `meta.info` that
#' identifies each sample's group.
#' @param meta.info Data frame of sample-level metadata, with sample names as
#' row names.
#' @param formula.str Character string containing the model formula.
#' @param correlation_block Character string or `NULL`. Name of the column in
#' `meta.info` identifying correlated samples or blocks.
#' @param model Character string specifying the fitting method: `"limma"` for
#' fixed effects or `"dream"` for random effects.
#' @param pairwise_contrasts Character vector of contrast expressions, such as
#' `"group1 - group2"`, to fit.
#'
#' @details
#' Dispatches to the Limma or DREAM fitting implementation. Both implementations
#' combine the group matrices and use the specified formula and contrasts.
#'
#' @return A list containing the following elements:
#' \item{d}{Absolute model statistics for each feature.}
#' \item{s}{Model standard errors for each feature after empirical Bayes
#' moderation.}
#'
#' @seealso \code{\link[limma]{lmFit}}, \code{\link[limma]{eBayes}},
#' \code{\link[limma]{topTable}},
#' \code{\link[limma]{makeContrasts}}
#'
#' @importFrom stats model.matrix formula
#' @importFrom dplyr bind_cols
#' @importFrom stringr str_split_fixed fixed
#' @importFrom limma makeContrasts lmFit contrasts.fit
#' @importFrom limma duplicateCorrelation
#' @importFrom utils combn
#' @importFrom variancePartition dream makeContrastsDream
#' @importFrom BiocParallel SerialParam
#'


fit_bootstrap <- function(x, group, meta.info, formula.str, 
                            correlation_block, model,
                            trend , robust, pairwise_contrasts) {
    if(model == "limma"){
        return(.Limma_bootstrap(x = x, group = group, 
        meta.info = meta.info, formula.str = formula.str, 
        correlation_block = correlation_block, 
        pairwise_contrasts = pairwise_contrasts, 
        trend = trend, robust = robust))

    }else if (model == "dream") {
        return(.Dream_bootstrap(x = x, group = group, 
        meta.info = meta.info, formula.str = formula.str, 
        correlation_block = correlation_block, 
        pairwise_contrasts = pairwise_contrasts , 
        trend = trend, robust = robust))
    }else{
        stop("Invalid model specified. Choose either 'limma' or 'dream'.")
    }
}

#' Fit Bootstrap Samples with Limma
#'
#' Combine group-specific bootstrap matrices, construct the design matrix,
#' account for optional correlation blocks, and fit contrasts with Limma.
#'
#' @param x List of feature-by-sample matrices, one per group in the bootstrap
#' replicate. Matrices must have matching feature names.
#' @param group Character string naming the column in `meta.info` that
#' identifies each sample's group.
#' @param meta.info Data frame of sample-level metadata, with sample names as
#' row names.
#' @param formula.str Character string containing the model formula.
#' @param correlation_block Character string or `NULL`. Name of the column in
#' `meta.info` identifying correlated samples or blocks.
#' @param pairwise_contrasts Character vector of contrast expressions, such as
#' `"group1 - group2"`, to fit.
#'
#' @return A list containing absolute feature statistics `d` and standard
#' errors `s`.
.Limma_bootstrap <-
    function(x, group, meta.info, formula.str,
    trend , robust, correlation_block , pairwise_contrasts) {

        
        ### Generate combined data and covariates for model fitting
        data_for_fit <- .generate_data_for_fit(x = x, meta.info = meta.info, 
        correlation_block = correlation_block)
        combined_data <- data_for_fit$combined_data
        covariates.p <- data_for_fit$covariates.p

        ### Create the design matrix based on the provided formula and covariates
        design.matrix <- model.matrix(formula(formula.str), data = covariates.p)
        colnames(design.matrix) <- make.names(colnames(design.matrix))
        stopifnot(identical(colnames(combined_data) , row.names(design.matrix)))

        ### Fit the linear model using Limma, accounting for correlation blocks if specified
        if (!is.null(correlation_block)) {
            corfit <- duplicateCorrelation(combined_data, design.matrix,
                                           block = 
                                             covariates.p[, correlation_block])
            fit <- lmFit(combined_data, design.matrix,
                         block = covariates.p[, correlation_block],
                         correlation = corfit$consensus.correlation)
        } else {
            fit <- lmFit(combined_data, design.matrix)
        }

        ### Compute contrasts and apply empirical Bayes moderation based on the number of groups
        cont_matrix <- makeContrasts(
            contrasts = pairwise_contrasts, levels =
                design.matrix
        )

        ### Fit the contrasts and apply empirical Bayes moderation
        fit2 <- limma::contrasts.fit(fit, cont_matrix)
        fit.ebayes <- limma::eBayes(fit2, trend = trend, robust = robust)

        d_values <- setNames(
            abs(as.numeric(fit.ebayes[["coefficients"]])),
            rownames(fit.ebayes[["coefficients"]]))
    
        if (is.null(correlation_block)) {
        stopifnot(identical(names(fit.ebayes[["s2.post"]]) ,  
                        row.names(fit.ebayes[["stdev.unscaled"]])))
        }

        s_values <- setNames(
            as.numeric(sqrt(fit.ebayes[["s2.post"]]) *
                            fit.ebayes[["stdev.unscaled"]]),
            row.names(fit.ebayes[["stdev.unscaled"]]))
    
        stopifnot(identical(names(s_values), names(d_values) ) )

        return(list(d = d_values, s = s_values))

        } 
             
        

#' Fit Bootstrap Samples with DREAM
#'
#' Combine group-specific bootstrap matrices, fit the specified mixed model
#' with DREAM, and return feature statistics and standard errors.
#'
#' @param x List of feature-by-sample matrices, one per group in the bootstrap
#' replicate. Matrices must have matching feature names.
#' @param group Character string naming the column in `meta.info` that
#' identifies each sample's group.
#' @param meta.info Data frame of sample-level metadata, with sample names as
#' row names.
#' @param formula.str Character string containing the model formula.
#' @param correlation_block Character string or `NULL`. Name of the column in
#' `meta.info` identifying correlated samples or blocks.
#' @param pairwise_contrasts Character vector of contrast expressions, such as
#' `"group1 - group2"`, to fit.
#'
#' @return A list containing absolute feature statistics `d` and standard
#' errors `s`.
.Dream_bootstrap <- function(x, group, meta.info, formula.str,
             trend , robust, correlation_block, pairwise_contrasts) {

        ### Generate combined data and covariates for model fitting
        data_for_fit <- .generate_data_for_fit(x = x, meta.info = meta.info, correlation_block = correlation_block)
        combined_data <- data_for_fit$combined_data
        covariates.p <- data_for_fit$covariates.p

        ### Create the design matrix based on the provided formula and covariates
        L <- variancePartition::makeContrastsDream(
            formula(formula.str), covariates.p,
            contrasts = pairwise_contrasts
        )

        ### Fit the model using DREAM
        fit <- variancePartition::dream(
            combined_data, formula(formula.str), covariates.p,
            L = L,
            BPPARAM = BiocParallel::SerialParam()
        )
        fit.ebayes <- variancePartition::eBayes(
            fit, trend = trend, robust = robust
        )
        d_values <- variancePartition::topTable(
            fit.ebayes,
            coef = pairwise_contrasts,
            number = "Inf",
            sort.by = "none"
        )
        d_values <- abs(d_values[["logFC"]])
        s_values <-
            as.numeric(sqrt(fit.ebayes[["s2.post"]]) *
                fit.ebayes[["stdev.unscaled"]][, 1])
        return(list(d = d_values, s = s_values))
            }


#' Assemble Bootstrap Data and Matching Metadata
#'
#' Combine group-specific matrices and construct a metadata row for every
#' resampled sample, preserving alignment and separating repeated draws of a
#' correlation block.
#'
#' @param x List of feature-by-sample matrices, one per group in the bootstrap
#' replicate. Matrices must have matching feature names.
#' @param meta.info Data frame of sample-level metadata, with sample names as
#' row names.
#' @param correlation_block Character string or `NULL`. Name of the column in
#' `meta.info` identifying correlated samples or blocks.
#'
#' @return A list containing the combined feature-by-sample matrix
#' `combined_data` and corresponding sample-level metadata `covariates.p`.
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

    ### if correlation_block is specified, modify the block IDs for duplicate samples

    if (!is.null(correlation_block)) {
        # A block drawn more than once is a separate unit each time; give
        # the n-th copy of each sample the n-th draw's block ID so that
        # duplicateCorrelation does not pool the copies into one block.
        real_names <- str_split_fixed(colnames(combined_data),
                                    fixed("."), 2)[, 1]
        draw <- stats::ave(seq_along(real_names), real_names,
                        FUN = seq_along)
        covariates.p[, correlation_block] <- paste(
            covariates.p[, correlation_block], draw, sep = "_")
    }

    stopifnot(identical(colnames(combined_data), row.names(covariates.p)))
    return(list(combined_data = combined_data, covariates.p = covariates.p))
}


#' Fit Permuted Metadata with Limma or DREAM
#'
#' Fit the selected linear model to the data matrix and permuted sample
#' metadata to obtain statistics for the permutation null distribution.
#'
#' @param x Feature-by-sample matrix with sample names as column names.
#' @param group Character string naming the column in `meta.info` that
#' identifies each sample's group.
#' @param meta.info Data frame of sample-level metadata, with sample names as
#' row names.
#' @param formula.str Character string containing the model formula.
#' @param correlation_block Character string or `NULL`. Name of the column in
#' `meta.info` identifying correlated samples or blocks.
#' @param model Character string specifying the fitting method: `"limma"` for
#' fixed effects or `"dream"` for random effects.
#' @param pairwise_contrasts Character vector of contrast expressions, such as
#' `"group1 - group2"`, to fit.
#'
#' @details
#' Dispatches to the Limma or DREAM fitting implementation, using the supplied
#' permuted metadata and contrast expressions.
#'
#' @return A list containing the following elements:
#' \item{d}{Absolute model statistics for each feature.}
#' \item{s}{Model standard errors for each feature after empirical Bayes
#' moderation.}
#'
#' @seealso \code{\link[limma]{lmFit}}, \code{\link[limma]{eBayes}},
#' \code{\link[limma]{topTable}},
#' \code{\link[limma]{makeContrasts}}
#'
#' @importFrom stats model.matrix formula
#' @importFrom dplyr bind_cols
#' @importFrom limma makeContrasts lmFit contrasts.fit
#' @importFrom limma duplicateCorrelation
#' @importFrom stringr str_split_fixed fixed
#' @importFrom utils combn
#' @importFrom variancePartition dream makeContrastsDream
#' @importFrom BiocParallel SerialParam
#'
#' 


fit_permutation <- function(x, group, meta.info, formula.str,
                             correlation_block, model,
                             pairwise_contrasts , trend , robust) {
    if(model == "limma"){
        return(.Limma_permutating(x = x, group = group, 
        meta.info = meta.info, formula.str = formula.str, 
        correlation_block = correlation_block, 
        pairwise_contrasts = pairwise_contrasts, 
        trend = trend, robust = robust))

    }else if (model == "dream") {
        return(.Dream_permutating(x = x, group = group, 
        meta.info = meta.info, formula.str = formula.str, 
        correlation_block = correlation_block, 
        pairwise_contrasts = pairwise_contrasts, 
        trend = trend, robust = robust))
    }else{
        stop("Invalid model specified. Choose either 'limma' or 'dream'.")
    }
}

#' Fit Permuted Metadata with Limma
#'
#' Construct the design matrix from the permuted metadata, fit the data with
#' Limma, account for optional correlation blocks, and apply the requested
#' contrasts.
#'
#' @param x Feature-by-sample matrix with sample names as column names.
#' @param group Character string naming the column in `meta.info` that
#' identifies each sample's group.
#' @param meta.info Data frame of sample-level metadata, with sample names as
#' row names.
#' @param formula.str Character string containing the model formula.
#' @param correlation_block Character string or `NULL`. Name of the column in
#' `meta.info` identifying correlated samples or blocks.
#' @param pairwise_contrasts Character vector of contrast expressions, such as
#' `"group1 - group2"`, to fit.
#'
#' @return A list containing absolute feature statistics `d` and standard
#' errors `s`.
.Limma_permutating <- function(x, group, meta.info, formula.str,
                             trend , robust, correlation_block, pairwise_contrasts) {

    ### Generate combined data and covariates for model fitting
    combined_data <- x
    covariates.p <- meta.info
    row.names(covariates.p) <- NULL
    
    ### Create the design matrix based on the provided formula and covariates
    design.matrix <-
        model.matrix(formula(formula.str), data = covariates.p)
    colnames(design.matrix) <-
        make.names(colnames(design.matrix))

    ### Fit the linear model using Limma, accounting for correlation blocks if specified
    if (!is.null(correlation_block)) {
        corfit <- duplicateCorrelation(combined_data, design.matrix,
                                       block = 
                                         covariates.p[, correlation_block])
        fit <- lmFit(combined_data, design.matrix,
                     block = covariates.p[, correlation_block],
                     correlation = corfit$consensus.correlation)
    } else {
        fit <- lmFit(combined_data, design.matrix)
    }

    ### Compute contrasts and apply empirical Bayes moderation based on the number of groups

    cont_matrix <- makeContrasts(
        contrasts = pairwise_contrasts,
        levels = design.matrix
    )

    ### Fit the contrasts and apply empirical Bayes moderation
    fit2 <- limma::contrasts.fit(fit, cont_matrix)
    fit.ebayes <- limma::eBayes(fit2, trend = trend, robust = robust)

    d_values <- setNames(
        abs(as.numeric(fit.ebayes[["coefficients"]])),
        rownames(fit.ebayes[["coefficients"]]))

    if (is.null(correlation_block)) {
    stopifnot(identical(names(fit.ebayes[["s2.post"]]) ,  
                    row.names(fit.ebayes[["stdev.unscaled"]])))
    }

    s_values <- setNames(
        as.numeric(sqrt(fit.ebayes[["s2.post"]]) *
                        fit.ebayes[["stdev.unscaled"]]),
        row.names(fit.ebayes[["stdev.unscaled"]]))

    stopifnot(identical(names(s_values), names(d_values) ) )

    return(list(d = d_values, s = s_values))
    
    }

#' Fit Permuted Metadata with DREAM
#'
#' Fit a linear mixed model to the data and permuted metadata with DREAM, then
#' return feature statistics and standard errors.
#'
#' @param x Feature-by-sample matrix with sample names as column names.
#' @param group Character string naming the column in `meta.info` that
#' identifies each sample's group.
#' @param meta.info Data frame of sample-level metadata, with sample names as
#' row names.
#' @param formula.str Character string containing the model formula.
#' @param correlation_block Character string or `NULL`. Name of the column in
#' `meta.info` identifying correlated samples or blocks.
#' @param pairwise_contrasts Character vector of contrast expressions, such as
#' `"group1 - group2"`, to fit.
#'
#' @return A list containing absolute feature statistics `d` and standard
#' errors `s`.
.Dream_permutating <- function(x, group, meta.info, formula.str,
            correlation_block, pairwise_contrasts, trend, robust) {

    ### Generate combined data and covariates for model fitting
    combined_data <- x
    covariates.p <- meta.info
    row.names(covariates.p) <- NULL
    
    L <- variancePartition::makeContrastsDream(
        formula(formula.str), covariates.p,
        contrasts = pairwise_contrasts
    )
    fit <- variancePartition::dream(
        combined_data, formula(formula.str), covariates.p,
        L = L,
        BPPARAM = BiocParallel::SerialParam()
    )
    fit.ebayes <- variancePartition::eBayes(
        fit, trend = trend, robust = robust
    )
    d_values <- variancePartition::topTable(fit.ebayes,
        coef = pairwise_contrasts,
        number = "Inf", sort.by = "none"
    )
    d_values <- abs(d_values[["logFC"]])
    s_values <- as.numeric(sqrt(fit.ebayes[["s2.post"]]) *
        fit.ebayes[["stdev.unscaled"]][, 1])

    return(list(d = d_values, s = s_values))

}
