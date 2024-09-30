#' Generate Bootstrap Samples with Optional Pairing
#'
#' This function generates bootstrap samples from the input metadata. It samples with replacement
#' within each group defined in the metadata, and optionally adjusts for paired groups.
#'
#' @param B Integer. The number of bootstrap samples to generate.
#' @param meta.info Data frame. Metadata containing sample information, where each row corresponds to a sample.
#' @param group.name Character. The name of the column in `meta.info` that defines the grouping variable for the samples.
#' @param paired Logical. If `TRUE`, the function ensures the bootstrap samples are paired between two groups.
#'
#' @details
#' The function works by resampling the row names of the metadata for each group separately. If `paired` is `TRUE`,
#' it assumes there are exactly two groups and samples the second group based on the positions of the first group to maintain pairing.
#'
#' @return A matrix of dimension \code{B} x \code{n}, where \code{n} is the number of samples. Each row corresponds
#' to a bootstrap sample, and each entry is a resampled row name from the metadata.
#'
#' @export
#' @examples
#' # Example usage:
#' set.seed(123)
#' meta.info <- data.frame(group = rep(c("A", "B"), each = 5), row.names = paste0("Sample", 1:10))
#' bootstrapS(B = 10, meta.info = meta.info, group.name = "group", paired = FALSE)
#'
#' # Paired bootstrap sampling
#' bootstrapS(B = 10, meta.info = meta.info, group.name = "group", paired = TRUE)
bootstrapS <- function(B, meta.info, group.name, paired) {
    groups <- meta.info[, group.name]
    bootsamples <- matrix(nrow = B, ncol = length(groups))
    for (i in seq_len(B)) {
        for (g in unique(groups)) {
            g.names <- row.names(meta.info)[which(groups == g)]
            bootsamples[i, which(groups == g)] <- sample(g.names, length(g.names), replace = TRUE)
        }
    }
    if (paired) {
        for (i in seq_len(B)) {
            g.names1 <- bootsamples[i, which(groups == unique(groups)[1])]
            g.names2 <- match(g.names1, row.names(meta.info)) + length(g.names1)
            bootsamples[i, which(groups == unique(groups)[2])] <- row.names(meta.info)[g.names2]
        }
    }
    return(bootsamples)
}


#' Generate Permutated Samples
#'
#' This function generates permuted samples by shuffling the row names of the metadata.
#'
#' @param meta.info Data frame. Metadata containing sample information, where each row corresponds to a sample.
#' @param B Integer. The number of permutations to generate.
#'
#' @details
#' The function creates a matrix where each row is a permuted version of the row names from `meta.info`.
#' This can be used to generate null distributions or perform randomization-based tests.
#'
#' @return A matrix of dimension \code{B} x \code{n}, where \code{n} is the number of samples (i.e., rows in `meta.info`).
#' Each row is a permutation of the row names of the metadata.
#'
#' @export
#' @examples
#' # Example usage:
#' set.seed(123)
#' meta.info <- data.frame(group = rep(c("A", "B"), each = 5), row.names = paste0("Sample", 1:10))
#' permutatedS(meta.info = meta.info, B = 10)
permutatedS <- function(meta.info, B)
{
    persamples <- matrix(nrow = B, ncol = nrow(meta.info))
    for (i in seq_len(B)) {
        persamples[i, ] <- sample(row.names(meta.info))
    }
    return(persamples)
}



#' Generate Stratified Bootstrap Samples for limRots
#'
#' This function generates stratified bootstrap samples based on the groupings and additional factors in the metadata.
#' The function ensures that samples are drawn proportionally based on strata defined by the interaction of factor columns in the metadata.
#'
#' @param B Integer. The number of bootstrap samples to generate.
#' @param meta.info Data frame. Metadata containing sample information, where each row corresponds to a sample. Factor columns in `meta.info` are used to define strata for sampling.
#' @param group.name Character. The name of the column in `meta.info` that defines the grouping variable for the samples.
#'
#' @details
#' The function works by first identifying the factors in the `meta.info` data frame that are used to create strata for sampling.
#' Within each group defined by `group.name`, the function samples according to the strata proportions, ensuring that samples are drawn
#' from the correct groups and strata in a proportional manner.
#'
#' @return A matrix of dimension \code{B} x \code{n}, where \code{n} is the number of samples. Each row corresponds
#' to a bootstrap sample, and each entry is a resampled row name from the metadata, stratified by group and additional factors.
#'
#' @export
#' @examples
#' # Example usage:
#' set.seed(123)
#' meta.info <- data.frame(group = rep(c(1, 2), each = 5),
#'     batch = rep(c("A", "B"), 5),
#'     row.names = paste0("Sample", 1:10))
#' meta.info$batch <- as.factor(meta.info$batch)
#' bootstrapSamples.limRots(B = 10, meta.info = meta.info, group.name = "group")
bootstrapSamples.limRots <- function(B, meta.info, group.name)
{
    labels <- as.numeric(meta.info[, group.name])
    samples <- matrix(nrow = B, ncol = length(labels))
    for (i in seq_len(B)) {
        for (label in unique(labels)) {
            pos <- which(labels == label)
            meta.info.pos <- meta.info[meta.info[, group.name] == label, ]
            meta.info.factors <- c()
            for (j in seq_len(ncol(meta.info))) {
                if (is.factor(meta.info.pos[, j])) {
                    meta.info.factors <- c(meta.info.factors, colnames(meta.info.pos)[j])
                }
            }
            if (is.null(meta.info.factors)) {
                samples <- bootstrapS(
                    B = B,
                    meta.info = meta.info,
                    group.name = group.name,
                    paired = FALSE
                )
                return(samples)
            }
            meta.info.factors <- meta.info.factors[meta.info.factors != group.name]
            meta.info.pos$stratum <- interaction(meta.info.pos[, meta.info.factors])
            stratum_sizes <- table(meta.info.pos$stratum)
            stratum_samples <- round(length(pos) * prop.table(stratum_sizes))
            sampled_indices <- unlist(lapply(names(stratum_samples), function(stratum) {
                stratum_indices <- row.names(meta.info.pos)[which(meta.info.pos$stratum == stratum)]
                sample(stratum_indices, stratum_samples[stratum], replace = TRUE)
            }))
            samples[i, pos] <- sampled_indices
        }
    }
    return(samples)
}
