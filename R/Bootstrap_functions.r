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
#' Generate Group-Wise Bootstrap Samples
#'
#' Generate bootstrap replicates by sampling samples with replacement
#' independently within each group.
#'
#' @param niter Integer. Number of bootstrap or permutation replicates to
#' generate.
#' @param meta.info Data frame containing sample-level metadata, with one row
#' per sample and sample names as row names.
#' @param group Character string naming the column in `meta.info` that
#' identifies each sample's group.
#'
#' @details
#' Within each group, sample names are drawn with replacement. Each replicate
#' contains the same number of sample positions from each group as the input.
#'
#' @return A character matrix with `niter` rows and one column per sample.
#' Each row is a bootstrap replicate containing resampled sample names.
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


#' Generate Stratified Group-Wise Bootstrap Samples
#'
#' Generate bootstrap replicates by sampling with replacement within groups
#' and requested metadata strata.
#'
#' @param niter Integer. Number of bootstrap or permutation replicates to
#' generate.
#' @param meta.info Data frame containing sample-level metadata, with one row
#' per sample and sample names as row names.
#' @param group Character string naming the column in `meta.info` that
#' identifies each sample's group.
#' @param strata Character vector naming columns used to define strata.
#' Only character or factor columns in `meta.info` are used; the `group`
#' column is excluded.
#'
#' @details
#' Within each group, samples are stratified by the interaction of the
#' selected columns. The number sampled from each stratum is proportional
#' to its size in that group, rounded to the nearest integer.
#'
#' @return A character matrix with `niter` rows and one column per sample.
#' Each row is a bootstrap replicate containing resampled sample names.
bootstrapSamples_limRots <- function(niter, meta.info, group, strata) {
    labels <- meta.info[, group]
    samples <- matrix(nrow = niter, ncol = length(labels))

    meta.info.factors <- names(meta.info)[sapply(meta.info, 
                                                function(x) is.character(x) || is.factor(x))]
    meta.info.factors <-
        meta.info.factors[meta.info.factors != group]

    meta.info.factors <- meta.info.factors[meta.info.factors %in% strata]

    message("Strata used: ", paste(meta.info.factors, collapse = ", "))

    for (i in seq_len(niter)) {
        for (label_i in unique(labels)) {
            pos <- which(labels == label_i)
            meta.info.pos <- meta.info[meta.info[, group] == label_i, ]
            

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

#' Generate Bootstrap Samples by Correlation Block
#'
#' Generate bootstrap replicates by resampling complete correlation blocks,
#' preserving the group composition of each block.
#'
#' @param niter Integer. Number of bootstrap or permutation replicates to
#' generate.
#' @param meta.info Data frame containing sample-level metadata, with one row
#' per sample and sample names as row names.
#' @param group Character string naming the column in `meta.info` that
#' identifies each sample's group.
#' @param correlation_block Character string naming the column in `meta.info`
#' that identifies correlated samples. Every block must contain exactly one
#' sample from each group.
#'
#' @details
#' In each replicate, correlation-block IDs are sampled with replacement.
#' For each selected block, the corresponding sample from every group is
#' included in that replicate, keeping the within-block pairing intact. The
#' function stops if any block does not contain exactly one sample per group
#' or a block-group sample cannot be matched.
#'
#' @return A character matrix with `niter` rows and one column per sample.
#' Each row is a bootstrap replicate containing resampled sample names, with
#' samples from each selected block kept together across groups.
#'

bootstrapSamples_limRots_block <- function(
    niter,
    meta.info,
    group,
    correlation_block
) {
    
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

#' Generate Stratified Permutation Samples
#'
#' Generate permutation replicates by shuffling group labels within the
#' specified metadata strata while leaving the other metadata unchanged.
#'
#' @param niter Integer. Number of bootstrap or permutation replicates to
#' generate.
#' @param meta.info Data frame containing sample-level metadata, with one row
#' per sample and sample names as row names.
#' @param group Character string naming the column in `meta.info` that
#' identifies each sample's group.
#' @param strata Character vector naming columns used to define strata.
#' Only character or factor columns in `meta.info` are used; the `group`
#' column is excluded.
#'
#' @details
#' Group labels are shuffled without replacement independently within each
#' stratum. Only character or factor columns named in `strata` are used to
#' define strata; all other columns are left unchanged.
#'
#' @return A list of length `niter`. Each element is a copy of `meta.info`
#' with the `group` labels permuted within strata.
permutationSamples_limRots <- function(niter, meta.info, group, strata) {

    meta.info.factors <- names(meta.info)[sapply(meta.info, 
                                                function(x) is.character(x) || is.factor(x))]
    meta.info.factors <-
        meta.info.factors[meta.info.factors != group]

    meta.info.factors <- meta.info.factors[meta.info.factors %in% strata]

    message("Strata used: ", paste(meta.info.factors, collapse = ", "))
    ## Create one combined stratum when multiple variables are supplied
    stratum <- interaction(
        meta.info[, meta.info.factors, drop = FALSE],
        drop = TRUE,
        lex.order = TRUE
    )
    ## Row positions belonging to each stratum
    stratum_indices <- split(seq_len(nrow(meta.info)), stratum)
    pSamples <- vector("list", niter)
    for (i in seq_len(niter)) {
        shuffle_df <- meta.info
        permuted_group <- meta.info[[group]]
        ## Shuffle group labels independently inside every stratum
        for (indices in stratum_indices) {
            permuted_group[indices] <- sample(
                meta.info[[group]][indices],
                size = length(indices),
                replace = FALSE
            )
        }
        shuffle_df[[group]] <- permuted_group
        pSamples[[i]] <- shuffle_df
    }
    return(pSamples)
}