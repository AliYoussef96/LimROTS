#' Generate Bootstrap Samples
#'
#' This function generates bootstrap samples from the input metadata. It samples
#'  with replacement within each group defined in the metadata, and optionally
#'  adjusts for paired groups.
#'
#' @param niter Integer. The number of bootstrap samples to generate.
#' @param meta.info Data frame. Metadata containing sample information, where
#' each row corresponds to a sample.
#' @param group.name Character. The name of the column in `meta.info` that
#' defines the grouping variable for the samples.
#' 
#' @details
#' The function works by resampling the row names of the metadata for each group
#'  separately.
#'
#' @return A matrix of dimension \code{niter} x \code{n}, where \code{n} is the
#' number of samples. Each row corresponds to a bootstrap sample, and each
#' entry is a resampled row name from the metadata.
#'
#' @examples
#' # Example usage:
#' set.seed(123)
#' meta.info <- data.frame(
#'     group = rep(c("A", "B"), each = 5),
#'     row.names = paste0("Sample", 1:10)
#' )
#' bootstrapS(
#'     niter = 10, meta.info = meta.info, group.name = "group"
#' )
#'
bootstrapS <- function(niter, meta.info, group.name) {
    groups <- meta.info[, group.name]
    bootsamples <- matrix(nrow = niter, ncol = length(groups))
    for (i in seq_len(niter)) {
        for (g in unique(groups)) {
            g.names <- row.names(meta.info)[which(groups == g)]
            bootsamples[i, which(groups == g)] <-
                sample(g.names, length(g.names), replace = TRUE)
        }
    }
    return(bootsamples)
}


#' Generate Stratified Bootstrap Samples for limRots
#'
#' This function generates stratified bootstrap samples based on the groupings
#' and additional factors in the metadata. The function ensures that samples
#' are drawn proportionally based on strata defined by the interaction of
#' factor columns in the metadata.
#'
#' @param niter Integer. The number of bootstrap samples to generate.
#' @param meta.info Data frame. Metadata containing sample information,
#' where each row corresponds to a sample. Factor columns in `meta.info`
#' are used to define strata for sampling.
#' @param group.name Character. The name of the column in `meta.info` that
#' defines the grouping variable for the samples.
#'
#' @details
#' The function works by first identifying the factors in the `meta.info` data
#' frame that are used to create strata for sampling. Within each group defined
#' by `group.name`, the function samples according to the strata proportions,
#' ensuring that samples are drawn from the correct groups and strata in a
#' proportional manner.
#'
#' @return A matrix of dimension \code{niter} x \code{n}, where \code{n} is the
#' number of samples. Each row corresponds to a bootstrap sample, and each
#' entry is a resampled row name from the metadata, stratified by group and
#' additional factors.
#'
#'
#'
bootstrapSamples_limRots <- function(niter, meta.info, group.name) {
    labels <- as.numeric(meta.info[, group.name])
    samples <- matrix(nrow = niter, ncol = length(labels))
    for (i in seq_len(niter)) {
        for (label in unique(labels)) {
            pos <- which(labels == label)
            meta.info.pos <- meta.info[meta.info[, group.name] == label, ]
            meta.info.factors <- c()
            for (j in seq_len(ncol(meta.info))) {
                if (is.factor(meta.info.pos[, j])) {
                    meta.info.factors <-
                        c(meta.info.factors, colnames(meta.info.pos)[j])
                }
            }
            if (is.null(meta.info.factors)) {
                samples <- bootstrapS(
                    niter = niter,
                    meta.info = meta.info,
                    group.name = group.name
                )
                return(samples)
            }
            meta.info.factors <-
                meta.info.factors[meta.info.factors != group.name]
            meta.info.pos$stratum <-
                interaction(meta.info.pos[, meta.info.factors])
            stratum_sizes <- table(meta.info.pos$stratum)
            stratum_samples <- round(length(pos) * prop.table(stratum_sizes))
            sampled_indices <-
                unlist(lapply(names(stratum_samples), function(stratum) {
                    stratum_indices <-
                        row.names(meta.info.pos)[which(meta.info.pos$stratum ==
                            stratum)]
                    sample(stratum_indices, stratum_samples[stratum],
                        replace = TRUE
                    )
                }))
            samples[i, pos] <- sampled_indices
        }
    }
    return(samples)
}
