# /*
#  * Adapted from ROTS Bioconductor (3.19, 2024) file bootstrapSamples.R
#  * Copyright on original version by: F. Seyednasrollah, T. Suomi,
#  * L.L. Elo (2024)
#  * Copyright on modifications (C) 2024-2025 A. Youssef
#  * Modifications: LimROTS::bootstrapS simplifies ROTS::bootstrapS 
#  * to unpaired group-wise resampling and returning 
#  * row names instead of indices.
#  */
#
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
#'
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
            meta.info.factors <-
                    meta.info.factors[meta.info.factors != group.name]
            if (is.null(meta.info.factors) | 
                length(meta.info.factors) == 0) {
                samples <- bootstrapS(
                    niter = niter,
                    meta.info = meta.info,
                    group.name = group.name
                )
                return(samples)
            }
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

#' Generate Stratified Bootstrap Samples with Correlation Blocks
#'
#' This function generates stratified bootstrap samples identical to
#' \code{bootstrapSamples_limRots}, but additionally supports correlation
#' blocks. When \code{correlation_block} is specified, all samples sharing
#' the same block ID are always selected together during resampling.
#' When \code{correlation_block} is NULL, the function delegates entirely
#' to \code{bootstrapSamples_limRots}.
#'
#' @param niter Integer. The number of bootstrap samples to generate.
#' @param meta.info Data frame. Metadata containing sample information,
#' where each row corresponds to a sample. Factor columns in `meta.info`
#' are used to define strata for sampling.
#' @param group.name Character. The name of the column in `meta.info` that
#' defines the grouping variable for the samples.
#' @param correlation_block Character or NULL. The name of a column in
#' `meta.info` that defines correlation blocks. Samples sharing the same
#' value in this column are always resampled together as a unit. If NULL,
#' the function behaves identically to \code{bootstrapSamples_limRots}.
#'
#' @details
#' The function follows the same logic as \code{bootstrapSamples_limRots}:
#' within each group defined by \code{group.name}, it identifies factor
#' columns to create strata, then samples proportionally within each stratum.
#' When \code{correlation_block} is not NULL, entire blocks (e.g., repeated
#' measures from the same subject) are resampled together as a unit instead
#' of individual samples.
#'
#' @return A matrix of dimension \code{niter} x \code{n}, where \code{n} is the
#' number of samples. Each row corresponds to a bootstrap sample, and each
#' entry is a resampled row name from the metadata, stratified by group and
#' additional factors.
#'
#'
#'
bootstrapSamples_limRots_block <- function(niter, meta.info, group.name,
                                           correlation_block = NULL) {
    if (is.null(correlation_block)) {
        return(bootstrapSamples_limRots(
            niter = niter,
            meta.info = meta.info,
            group.name = group.name
        ))
    }
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
            meta.info.factors <-
                    meta.info.factors[meta.info.factors != group.name]
            if (is.null(meta.info.factors) |
                length(meta.info.factors) == 0) {
                block_ids <- meta.info.pos[, correlation_block]
                unique_blocks <- unique(block_ids)
                collected <- c()
                while (length(collected) < length(pos)) {
                    b <- sample(unique_blocks, 1)
                    members <-
                        row.names(meta.info.pos)[which(block_ids == b)]
                    collected <- c(collected, members)
                }
                samples[i, pos] <- collected[seq_len(length(pos))]
            } else {
                meta.info.pos$stratum <-
                    interaction(meta.info.pos[, meta.info.factors])
                stratum_sizes <- table(meta.info.pos$stratum)
                stratum_samples <-
                    round(length(pos) * prop.table(stratum_sizes))
                sampled_indices <-
                    unlist(lapply(names(stratum_samples), function(stratum) {
                        stratum_indices_mask <-
                            which(meta.info.pos$stratum == stratum)
                        stratum_meta <- meta.info.pos[stratum_indices_mask, ]
                        block_ids <- stratum_meta[, correlation_block]
                        unique_blocks <- unique(block_ids)
                        collected <- c()
                        while (length(collected) <
                               stratum_samples[stratum]) {
                            b <- sample(unique_blocks, 1)
                            members <-
                                row.names(stratum_meta)[which(block_ids == b)]
                            collected <- c(collected, members)
                        }
                        collected[seq_len(stratum_samples[stratum])]
                    }))
                samples[i, pos] <- sampled_indices
            }
        }
    }
    return(samples)
}

#' Generate Stratified Bootstrap Samples for Cox limRots with Correlation Blocks
#'
#' This function generates stratified bootstrap samples similar to
#' \code{bootstrapSamples_limRots}, but additionally supports correlation
#' blocks. When \code{correlation_block} is specified, all samples sharing
#' the same block ID are always selected together during resampling.
#'
#' @param niter Integer. The number of bootstrap samples to generate.
#' @param meta.info Data frame. Metadata containing sample information,
#' where each row corresponds to a sample. Factor columns in `meta.info`
#' are used to define strata for sampling.
#' @param correlation_block Character or NULL. The name of a column in
#' `meta.info` that defines correlation blocks. Samples sharing the same
#' value in this column are always resampled together as a unit.
#'
#' @details
#' When \code{correlation_block} is not NULL, the function groups samples by
#' their block ID within each group/stratum and resamples entire blocks with
#' replacement, so that correlated samples (e.g., repeated measures from the
#' same subject) are always kept together.
#'
#' @return A matrix of dimension \code{niter} x \code{n}, where \code{n} is the
#' number of samples. Each row corresponds to a bootstrap sample, and each
#' entry is a resampled row name from the metadata.
#'
#'
#'
bootstrapSamples_limRots_cox <- function(niter, meta.info,
                                         correlation_block = NULL) {
  samples <- matrix(nrow = niter, ncol = nrow(meta.info))
  
  meta.info.factors <- c()
  for (j in seq_len(ncol(meta.info))) {
    if (is.factor(meta.info[, j])) {
      meta.info.factors <-
        c(meta.info.factors, colnames(meta.info)[j])
    }
  }
  
  # Always add "event" to stratification factors
  meta.info.factors <- c(meta.info.factors, "event")
  meta.info[,"event"] <- as.factor(meta.info[,"event"])
  
  # Stratified bootstrap
  meta.info$stratum <-
    interaction(meta.info[, meta.info.factors])
  stratum_sizes <- table(meta.info$stratum)
  stratum_samples <-
    round(nrow(meta.info) * prop.table(stratum_sizes))
  
  for (i in seq_len(niter)) {
    
    
    if (is.null(correlation_block)) {
      # Sample individuals directly
      sampled_indices <-
        unlist(lapply(names(stratum_samples), function(stratum) {
          stratum_indices <-
            row.names(meta.info)[which(meta.info$stratum == stratum)]
          sample(stratum_indices, stratum_samples[stratum],
                 replace = TRUE
          )
        }))
    } else {
      # Sample correlation blocks
      sampled_indices <-
        unlist(lapply(names(stratum_samples), function(stratum) {
          stratum_mask <-
            which(meta.info$stratum == stratum)
          stratum_meta <- meta.info[stratum_mask, ]
          block_ids <- stratum_meta[, correlation_block]
          unique_blocks <- unique(block_ids)
          # Resample blocks until we reach the target count
          collected <- c()
          while (length(collected) <
                 stratum_samples[stratum]) {
            b <- sample(unique_blocks, 1)
            members <-
              row.names(stratum_meta)[which(block_ids == b)]
            collected <- c(collected, members)
          }
          # Trim to exact target size
          collected[seq_len(stratum_samples[stratum])]
        }))
    }
    samples[i, ] <- sampled_indices
  }
  return(samples)
}
