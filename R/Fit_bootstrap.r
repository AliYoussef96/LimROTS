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
        data_for_fit <- .generate_data_for_fit(x = x, meta.info = meta.info, correlation_block = correlation_block)
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
    
        stopifnot(identical(names(fit.ebayes[["s2.post"]]) ,  
                        row.names(fit.ebayes[["stdev.unscaled"]])))

        s_values <- setNames(
            as.numeric(sqrt(fit.ebayes[["s2.post"]]) *
                            fit.ebayes[["stdev.unscaled"]]),
            names(fit.ebayes[["s2.post"]]))
    
        stopifnot(identical(names(s_values), names(d_values) ) )

        return(list(d = d_values, s = s_values))

        } 
             
        

.Dream_bootstrap <- function(x, group, meta.info, formula.str,
             correlation_block) {

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