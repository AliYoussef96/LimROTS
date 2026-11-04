# Copyright (C) 2025 Ali Youssef.
# This program is licensed under GPL (>= 2).
#
#' Parallel processing handling function
#'
#' @param data A \code{SummarizedExperiment} object or a matrix where rows
#' represent features (e.g., genes, proteins) and columns represent samples.
#' The values should be log-transformed.
#' @param meta.info A data frame containing sample-level metadata, where each
#' row corresponds to a sample. It should include the grouping variable
#' specified in \code{group}. If \code{x} is a \code{SummarizedExperiment}
#' object, \code{meta.info} must be a vector of the metadata needed for the
#' model to run and can be retrieved using \code{colData()}.
#' @param group A string specifying the column in \code{meta.info} that
#' represents the groups or conditions for comparison.
#' @param BPPARAM A parallel BPPARAM object for distributed computation.
#' @param formula.str A formula string used when covariates are present in meta.
#' info for modeling. It should include "~ 0 + ..." to exclude the
#' intercept from the model.
#' @param samples bootstrapped samples matrix
#' @param groups groups information from `meta.info`
#' @param pSamples a permutated list of samples
#' @param correlation_block Character or NULL. The name of a column in
#' \code{meta.info} that defines correlation blocks. Passed to
#' \code{Limma_bootstrap} and \code{Limma_permutating} to account for
#' within-block correlation during model fitting via
#' \code{duplicateCorrelation}. If NULL, standard independent fitting is used.
#' @param model A character string specifying the modeling approach to be used. 
#' limma for fixed effects and dream for random effects. Default is "limma".
#' @param pairwise_contrasts A character vector specifying the pairwise contrasts to be computed. 
#' Each contrast should be in the format "group1 - group2". for multi-group analysis, 
#' this will be generated automatically.
#' 
#' @return A list containing: \code{D, S, pD, pS} for bootstrapped data and
#'  for permuted data.
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
                          pairwise_contrasts) {
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
        pairwise_contrasts = pairwise_contrasts
    )
    export_funcs <- list(fit_bootstrap = fit_bootstrap,
                         fit_permutation = fit_permutation,
                         .Limma_bootstrap = .Limma_bootstrap,
                         .Dream_bootstrap  = .Dream_bootstrap,
                         .generate_data_for_fit = .generate_data_for_fit,
                         .Limma_permutating = .Limma_permutating,
                         .Dream_permutating = .Dream_permutating)
    for (name in names(export_vars)) {
        assign(name, export_vars[[name]])
    }
    for (name in names(export_funcs)) {
        assign(name, export_funcs[[name]])
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
            pairwise_contrasts = pairwise_contrasts
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
            pairwise_contrasts = pairwise_contrasts
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



# Copyright (C) 2025 Ali Youssef.
# This program is licensed under GPL (>= 2).
#
#' Perform Linear Modeling with Covariates using Limma
#'
#' This function performs linear modeling using the Limma package, incorporating
#'  covariates in the model fitting process. It is designed to handle both
#'  two-group comparisons and multi-group settings with covariates.
#'
#' @param x A list containing two or more data matrices where rows represent
#' features (e.g., genes, proteins) and columns represent samples. The list
#' should contain at least two matrices for pairwise group comparison.
#' @param group A character string indicating the name of the group
#' variable in`meta.info` to be used in the analysis.
#' @param meta.info A data frame containing the metadata for the samples.
#' This includes sample grouping and any covariates to be included in the model.
#' @param formula.str A string specifying the formula to be used in model
#' fitting. It should follow the standard R formula syntax
#' (e.g., `~ covariate1 + covariate2`). If the formula contains a random
#' effects term (i.e., includes `|`, e.g., `~ 0 + group + (1|subject)`),
#' the function automatically routes to DREAM analysis via
#' \code{\link[variancePartition]{dream}}. The DREAM path is restricted to
#' two-group comparisons; passing more than two groups with a random term
#' raises an error.
#' @param trend A logical value indicating whether to allow for an
#' intensity-dependent trend in the prior variance.
#' @param robust A logical value indicating whether to use a robust fitting
#' procedure to protect against outliers.
#' @param correlation_block Character or NULL. The name of a column in
#' \code{meta.info} that defines correlation blocks. When not NULL,
#' within-block correlation is estimated via \code{duplicateCorrelation} and
#' accounted for in \code{lmFit}. If NULL, standard independent fitting is used.
#' @param pairwise_contrasts A character vector specifying the pairwise contrasts to be computed. 
#' Each contrast should be in the format "group1 - group2". for multi-group analysis, 
#' this will be generated automatically.
#' @param model A character string specifying the modeling approach to be used. 
#' limma for fixed effects and dream for random effects. Default is "limma".
#'
#' @details
#' This function combines the data matrices from different groups and fits a
#' linear model using covariates provided in the `meta.info`. For two-group
#' comparisons, the function computes contrasts between the two groups and
#' applies empirical Bayes moderation. For multi-group analysis with a single
#' covariate, pairwise contrasts are computed, and the moderated F-statistic is
#' calculated for each feature.
#'
#' When `formula.str` contains a `|` character (random effects term), the
#' function bypasses the standard limma workflow and instead uses
#' \code{\link[variancePartition]{dream}} to fit a linear mixed model. This
#' path is only supported for two-group comparisons. Passing more than two
#' groups with a random-effects formula raises an informative error directing
#' the user to use the \code{correlation_block} argument with a fixed-effects
#' formula instead.
#'
#' @return A list containing the following elements:
#' \item{d}{A vector of the test statistics (log-fold changes or F-statistics)
#' for each feature.}
#' \item{s}{A vector of the standard deviations for each feature, adjusted by t
#' he empirical Bayes procedure.}
#' \item{corrected.logfc}{The log-fold changes for each feature after fitting
#' the model.}
#'
#' @seealso \code{\link[limma]{lmFit}},
#' \code{\link[limma]{eBayes}},
#' \code{\link[limma]{topTable}},
#' \code{\link[limma]{makeContrasts}},
#' \code{\link[variancePartition]{dream}},
#' \code{\link[variancePartition]{makeContrastsDream}}
#'
#'
#'
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
        return(.Limma_fit_original(x = x, group = group, meta.info = meta.info, formula.str = formula.str, 
        correlation_block = correlation_block, pairwise_contrasts = pairwise_contrasts,
        trend = trend, robust = robust ))
    }else if (model == "dream") {
        return(.Dream_fit_original(x = x, group = group, meta.info = meta.info, formula.str = formula.str, 
        correlation_block = correlation_block, pairwise_contrasts = pairwise_contrasts,
        trend = trend, robust = robust ))
    }else{
        stop("Invalid model specified. Choose either 'limma' or 'dream'.")
    }
}

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


# Copyright (C) 2025 Ali Youssef.
# This program is licensed under GPL (>= 2).
#
#' Perform Linear Modeling with Covariates using Limma
#'
#' This function performs linear modeling using the Limma package while
#' accounting for covariates specified
#' in the `meta.info`. It supports two-group comparisons and multi-group
#' analysis, incorporating covariates
#' through a design matrix.
#'
#' @param x A list containing two or more data matrices where rows represent
#' features (e.g., genes, proteins) and columns represent samples. The list
#' should contain at least two matrices for pairwise group comparison.
#' @param group A character string indicating the name of the group
#' variable in `meta.info` to be used in the analysis.
#' @param meta.info A data frame containing the metadata for the samples.
#' This includes sample grouping and any covariates to be included in the model.
#' @param formula.str A string specifying the formula to be used in model
#' fitting. It should follow the standard R formula syntax
#' (e.g., `~ covariate1 + covariate2`).
#' @param correlation_block Character or NULL. The name of a column in
#' \code{meta.info} that defines correlation blocks. When not NULL,
#' within-block correlation is estimated via \code{duplicateCorrelation} and
#' accounted for in \code{lmFit}. If NULL, standard independent fitting is used.
#' @param pairwise_contrasts A character vector specifying the pairwise contrasts to be computed. 
#' Each contrast should be in the format "group1 - group2". for multi-group analysis, 
#' this will be generated automatically.
#' @param model A character string specifying the modeling approach to be used. 
#' limma for fixed effects and dream for random effects. Default is "limma".
#'
#' @details
#' This function first combines the data matrices from different groups and
#' prepares a design matrix based on the covariates specified in `meta.info`
#' using the provided formula. It fits a linear model using Limma,
#' computes contrasts between groups, and applies empirical Bayes moderation.
#' For two-group comparisons, the function returns log-fold changes and
#' associated statistics. In multi-group settings with a single covariate,
#' it calculates pairwise contrasts and moderated F-statistics.
#'
#' @return A list containing the following elements:
#' \item{d}{A vector of the test statistics (log-fold changes or F-statistics)
#' for each feature.}
#' \item{s}{A vector of the standard deviations for each feature, adjusted by
#' the empirical Bayes procedure.}
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
                            pairwise_contrasts) {
    if(model == "limma"){
        return(.Limma_bootstrap(x, group, meta.info, formula.str, 
        correlation_block, pairwise_contrasts))
    }else if (model == "dream") {
        return(.Dream_bootstrap(x, group, meta.info, formula.str, 
        correlation_block, pairwise_contrasts))
    }else{
        stop("Invalid model specified. Choose either 'limma' or 'dream'.")
    }
}

.Limma_bootstrap <-
    function(x, group, meta.info, formula.str,
             correlation_block , pairwise_contrasts) {

        
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
        fit.ebayes <- limma::eBayes(fit2, trend = FALSE, robust = FALSE)

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
             
        

.Dream_bootstrap <- function(x, group, meta.info, formula.str,
             correlation_block, pairwise_contrasts) {

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
            fit, trend = FALSE, robust = FALSE
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


# Copyright (C) 2025 Ali Youssef.
# This program is licensed under GPL (>= 2).
#
#' Perform Permutation-Based Linear Modeling with Covariates using Limma
#'
#' This function performs linear modeling using the Limma package with
#' permutation of the covariates to evaluate the test statistics under random
#' assignments. It handles two-group comparisons and multi-group settings.
#'
#' @param x A data matrices where rows represent
#' features (e.g., genes, proteins) and columns represent samples. The list
#' should contain at least two matrices for pairwise group comparison.
#' @param group A character string indicating the name of the group
#' variable in `meta.info` to be used in the analysis.
#' @param meta.info A data frame containing the metadata for the samples.
#' This includes sample grouping and any covariates to be included in the model.
#' @param formula.str A string specifying the formula to be used in model
#' fitting. It should follow the standard R formula syntax
#' (e.g., `~ covariate1 + covariate2`).
#' @param correlation_block Character or NULL. The name of a column in
#' \code{meta.info} that defines correlation blocks. When not NULL,
#' within-block correlation is estimated via \code{duplicateCorrelation} and
#' accounted for in \code{lmFit}. If NULL, standard independent fitting is used.
#'
#' @details
#' This function combines the data matrices from different groups and permutes
#' the covariates from `meta.info`before fitting a linear model using Limma.
#' Permutation helps assess how the covariates behave under random conditions,
#' providing a null distribution of the test statistics. For two-group
#' comparisons, the function computes contrasts between the two groups and
#' applies empirical Bayes moderation. For multi-group analysis with a single
#' covariate, pairwise contrasts are computed, and the moderated F-statistic is
#' calculated for each feature.
#'
#' @return A list containing the following elements:
#' \item{d}{A vector of the test statistics (log-fold changes or F-statistics)
#' for each feature.}
#' \item{s}{A vector of the standard deviations for each feature, adjusted by
#' the empirical Bayes procedure.}
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
                             pairwise_contrasts) {
    if(model == "limma"){
        return(.Limma_permutating(x, group, meta.info, formula.str, 
        correlation_block, pairwise_contrasts))
    }else if (model == "dream") {
        return(.Dream_permutating(x, group, meta.info, formula.str, 
        correlation_block, pairwise_contrasts))
    }else{
        stop("Invalid model specified. Choose either 'limma' or 'dream'.")
    }
}

.Limma_permutating <- function(x, group, meta.info, formula.str,
                             correlation_block, pairwise_contrasts) {

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
    fit.ebayes <- limma::eBayes(fit2, trend = FALSE, robust = FALSE)

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


.Dream_permutating <- function(x, group, meta.info, formula.str,
            correlation_block, pairwise_contrasts) {

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
        fit, trend = FALSE, robust = FALSE
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

