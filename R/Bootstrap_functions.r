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
#' @param group Character. The name of the column in `meta.info` that
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
bootstrapS <- function(niter, meta.info, group) {
    groups <- meta.info[, group]
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
#' @param group Character. The name of the column in `meta.info` that
#' defines the grouping variable for the samples.
#'
#' @details
#' The function works by first identifying the factors in the `meta.info` data
#' frame that are used to create strata for sampling. Within each group defined
#' by `group`, the function samples according to the strata proportions,
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
bootstrapSamples_limRots <- function(niter, meta.info, group) {
    labels <- as.numeric(meta.info[, group])
    samples <- matrix(nrow = niter, ncol = length(labels))
    for (i in seq_len(niter)) {
        for (label in unique(labels)) {
            pos <- which(labels == label)
            meta.info.pos <- meta.info[meta.info[, group] == label, ]
            meta.info.factors <- c()
            for (j in seq_len(ncol(meta.info))) {
                if (is.factor(meta.info.pos[, j])) {
                    meta.info.factors <-
                        c(meta.info.factors, colnames(meta.info.pos)[j])
                }
            }
            meta.info.factors <-
                    meta.info.factors[meta.info.factors != group]
            if (is.null(meta.info.factors) | 
                length(meta.info.factors) == 0) {
                samples <- bootstrapS(
                    niter = niter,
                    meta.info = meta.info,
                    group = group
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
#' This function generates bootstrap samples while preserving the grouping
#' structure and, when specified, the correlation between paired samples.
#' When \code{correlation_block} is provided and block IDs occur in both
#' groups, complete blocks are resampled so that the corresponding samples
#' from both groups are always selected together.
#'
#' @param niter Integer. The number of bootstrap samples to generate.
#' @param meta.info Data frame. Metadata containing sample information,
#' where each row corresponds to one sample. Row names must contain the
#' sample names. Factor columns may be used to define strata when
#' \code{stratum = TRUE}.
#' @param group Character. The name of the column in \code{meta.info} that
#' defines the grouping variable. For repeated-measures block resampling,
#' the grouping variable is expected to contain two levels.
#' @param correlation_block Character or NULL. The name of the column in
#' \code{meta.info} that identifies correlated samples. Samples sharing the
#' same block ID across the two groups are resampled together. Each block
#' must contain exactly one sample from each group. If NULL, the function
#' delegates entirely to \code{bootstrapSamples_limRots}.
#' @param stratum Logical. Indicates whether stratified block resampling
#' should be performed. If TRUE, factor columns in \code{meta.info} are
#' automatically used to define strata, excluding the columns specified by
#' \code{group} and \code{correlation_block}. If FALSE, all correlation
#' blocks are sampled from a single stratum.
#'
#' @details
#' When \code{correlation_block} is NULL, the function behaves identically
#' to \code{bootstrapSamples_limRots}.
#'
#' When correlation blocks span the levels of \code{group}, the function
#' performs block-level bootstrap resampling. A block ID is sampled once,
#' and the corresponding sample from each group is included in the same
#' bootstrap replicate. This preserves the paired or repeated-measures
#' structure of the data.
#'
#' Each correlation block must have exactly one sample in every group.
#' The function stops with an error if a block is missing from one group or
#' has more than one sample in the same group.
#'
#' If \code{stratum = TRUE}, all factor columns in \code{meta.info}, except
#' \code{group} and \code{correlation_block}, are used to define the strata.
#' Each stratification variable must be constant across all samples belonging
#' to the same correlation block. Complete blocks are then sampled with
#' replacement separately within each stratum. The original number of blocks
#' in every stratum is preserved in each bootstrap replicate.
#'
#' If no additional factor columns are available, or if
#' \code{stratum = FALSE}, all correlation blocks are sampled together
#' without additional stratification.
#'
#' @return A matrix of dimension \code{niter} x \code{n}, where \code{n} is
#' the number of samples in \code{meta.info}. Each row represents one
#' bootstrap replicate, and each entry contains a resampled row name from
#' the metadata. Samples belonging to the same correlation block are selected
#' together across the two groups.
#'

bootstrapSamples_limRots_block <- function(
    niter,
    meta.info,
    group,
    correlation_block = NULL
) {
    
    # If no correlation block is provided, use the original
    # limRots bootstrap procedure.
    if (is.null(correlation_block)) {
        return(
            bootstrapSamples_limRots(
                niter = niter,
                meta.info = meta.info,
                group = group
            )
        )
    }

    message(
        "Using correlation blocks for bootstrap sampling. ",
        "All samples sharing the same block ID will be sampled together. ",
        "No stratification will be applied."
    )
    
    groups <- as.character(meta.info[, group])
    blocks <- as.character(meta.info[, correlation_block])
    
    if (anyNA(blocks)) {
        stop("The correlation-block column contains missing values.")
    }
    
    # Every block must contain exactly one sample from each group.
    if (any(table(blocks, groups) != 1L)) {
        stop(
            "Every correlation block must contain exactly one sample ",
            "from each of the two groups."
        )
    }
    
    unique_groups <- unique(groups)    
    unique_blocks <- unique(blocks)
    
    # Connect each block-group combination to its sample.
    key <- paste(blocks, groups, sep = "\r")
    
    # Each row represents one bootstrap iteration.
    # Each column represents one sample position.
    samples <- matrix(nrow = niter, ncol = nrow(meta.info))    
    for (i in seq_len(niter)) {
        
        # Sample complete pairs with replacement.
        drawn_blocks <- sample(
            unique_blocks,
            size = length(unique_blocks),
            replace = TRUE
        )
        
        # Retrieve the corresponding sample from each group using
        # the same drawn-block order, preserving the pairing.
        for (g in unique_groups) {
            
            group_positions <- which(groups == g)
            
            matched_positions <- match(
                paste(drawn_blocks, g, sep = "\r"),
                key
            )
            
            # Defensive check: all block-group combinations
            # should have corresponding samples.
            if (anyNA(matched_positions)) {
                stop(
                    "A sample could not be found for one or more ",
                    "block-group combinations."
                )
            }
            
            samples[i, group_positions] <-
                row.names(meta.info)[matched_positions]
        }
    }
    
    return(samples)
}


# bootstrapSamples_limRots_block <- function(niter, meta.info, group,
#                                            stratum = TRUE,
#                                            correlation_block = NULL) {
#   
#   if (is.null(correlation_block)) {
#     return(bootstrapSamples_limRots(
#       niter = niter,
#       meta.info = meta.info,
#       group = group
#     ))
#   }
#   
#   groups <- as.character(meta.info[, group])
#   blocks <- as.character(meta.info[, correlation_block])
#   
#   if (any(tapply(groups, blocks, function(g) length(unique(g)) > 1))) {
#     # Repeated measures: blocks span levels of `group`, so each replicate
#     # draws whole blocks and takes every drawn block's sample from each
#     # level, keeping the pairing intact.
#     if (any(table(blocks, groups) != 1)) {
#       stop("correlation_block spans levels of '", group, "', but not ",
#            "every block has exactly one sample in each level.")
#     }
#     
#     unique_blocks <- unique(blocks)
#     
#     if(stratum == TRUE){
#       
#       # Detect factor columns defining the strata.
#       meta.info.factors <- c()
#       
#       for (j in seq_len(ncol(meta.info))) {
#         if (is.factor(meta.info[, j])) {
#           meta.info.factors <- c(
#             meta.info.factors,
#             colnames(meta.info)[j]
#           )
#         }
#       }
#       
#       # Group and correlation block must not define strata.
#       meta.info.factors <- meta.info.factors[
#         !meta.info.factors %in% c(group, correlation_block)
#       ]
#     }else{
#       meta.info.factors = character(0)
#     }
#     
#     # Assign each block to a stratum.
#     if (length(meta.info.factors) == 0) {
#       blocks_by_stratum <- list(all = unique_blocks)
#     } else {
#       block_strata <- sapply(unique_blocks, function(b) {
#         block_meta <- meta.info[
#           blocks == b,
#           meta.info.factors,
#           drop = FALSE
#         ]
#         
#         # Strata must be the same for all samples in one block.
#         if (any(vapply(
#           block_meta,
#           function(x) length(unique(x)) != 1,
#           logical(1)
#         ))) {
#           stop(
#             "Strata variables must be constant within ",
#             "correlation block '", b, "'."
#           )
#         }
#         
#         paste(
#           vapply(block_meta[1, , drop = FALSE], as.character, character(1)),
#           collapse = "."
#         )
#       })
#       
#       blocks_by_stratum <- split(
#         unique_blocks,
#         block_strata
#       )
#     }
#     
#     key <- paste(blocks, groups, sep = "\r")
#     samples <- matrix(nrow = niter, ncol = nrow(meta.info))
#     
#     for (i in seq_len(niter)) {
#       
#       # Sample the same number of complete blocks within each stratum.
#       drawn <- unlist(
#         lapply(blocks_by_stratum, function(x) {
#           sample(x, length(x), replace = TRUE)
#         }),
#         use.names = FALSE
#       )
#       
#       for (g in unique(groups)) {
#         pos <- which(groups == g)
#         
#         samples[i, pos] <- row.names(meta.info)[
#           match(paste(drawn, g, sep = "\r"), key)
#         ]
#       }
#     }
#     
#     return(samples)
#   }
# }

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
