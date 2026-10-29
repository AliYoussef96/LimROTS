# /*
#  * Adapted from ROTS Bioconductor (3.19, 2024) file ROTS.R
#  * Copyright on original version by: F. Seyednasrollah, T. Suomi,
#  * L.L. Elo (2024)
#  * Copyright on modifications (C) 2024-2025 A. Youssef
#  * Modifications: ROTS.R lines 143–260, ~35 lines were separated
#  * into its own optimization step function.
#  */
#
#' Optimize Parameters Based on Overlap Calculations
#'
#' This function optimizes parameters by calculating overlaps between observed
#' and permuted data for multiple values of a smoothing constant (`ssq`) and a
#' single-label replicate (SLR) comparison.
#'
#' @param niter Integer. Number of bootstrap samples or resampling iterations.
#' @param ssq Numeric vector. Smoothing constants to be evaluated.
#' @param N Integer vector. Number of top values to consider for overlap
#' calculation.
#' @param D Numeric matrix. Observed data values.
#' @param S Numeric matrix. Standard errors or related values for observed data.
#' @param pD Numeric matrix. Permuted data values.
#' @param pS Numeric matrix. Standard errors or related values for
#' permuted data.
#'
#' @details
#' The function calculates overlaps for a range of smoothing constants and
#' identifies the optimal set of parameters by maximizing a z-score-based
#' metric, which compares the overlap of observed data to permuted data.
#' It computes overlap matrices for both observed (`D` and `S`) and permuted
#' (`pD` and `pS`) data and returns the optimal parameters based on the
#' highest z-score.
#'
#' @return A list containing the optimal parameters:
#' \itemize{
#'   \item \code{a1}: Optimal smoothing constant or 1 for SLR.
#'   \item \code{a2}: SLR flag (1 if smoothing constant is optimal,
#'   0 if SLR is optimal).
#'   \item \code{k}: Optimal number of top values to consider for overlap.
#'   \item \code{R}: Optimal overlap value.
#'   \item \code{Z}: Optimal z-score.
#'   \item \code{ztable}: Matrix of z-scores for all evaluated parameters.
#' }
#'
#'
#'
Optimizing <- function(niter, k, N, D, S, pD, pS) {
    ssq <- c(
        seq(0, 20) / 100, 
        seq(11, 50) / 50, 
        seq(6, 25) / 5
        )
    N <- c(
        seq(1, 20) * 5,
        seq(11, 50) * 10,
        seq(21, 40) * 25,
        seq(11, 1000) * 100
    )
    N <- N[N < K]

    reprotable <- matrix(nrow = length(ssq) + 1, ncol = length(N))
    colnames(reprotable) <- N
    row.names(reprotable) <- c(ssq, "slr")
    reprotable.P <- matrix(nrow = length(ssq) + 1, ncol = length(N))
    colnames(reprotable.P) <- N
    row.names(reprotable.P) <- c(ssq, "slr")
    reprotable.sd <- matrix(nrow = length(ssq) + 1, ncol = length(N))
    colnames(reprotable.sd) <- N
    row.names(reprotable.sd) <- c(ssq, "slr")
    for (i in seq_len(length(ssq))) {
        overlaps <- matrix(0, nrow = niter, ncol = length(N))
        overlaps.P <- matrix(0, nrow = niter, ncol = length(N))
        cResults <- .calOverlaps(
            D, S, pD, pS, nrow(D), as.integer(N), length(N),
            ssq[i], as.integer(niter), overlaps, overlaps.P
        )
        reprotable[i, ] <- colMeans(cResults[["overlaps"]])
        reprotable.P[i, ] <- colMeans(cResults[["overlaps_P"]])
        reprotable.sd[i, ] <- sqrt(rowSums((t(cResults[["overlaps"]]) -
            reprotable[i, ])^2) /
            (nrow(cResults[["overlaps"]]) - 1))
    }
    i <- length(ssq) + 1
    overlaps <- matrix(0, nrow = niter, ncol = length(N))
    overlaps.P <- matrix(0, nrow = niter, ncol = length(N))
    cResults <- .calOverlaps_slr(
        D, pD, nrow(D), as.integer(N), length(N),
        as.integer(niter), overlaps, overlaps.P
    )
    reprotable[i, ] <- colMeans(cResults[["overlaps"]])
    reprotable.P[i, ] <- colMeans(cResults[["overlaps_P"]])
    reprotable.sd[i, ] <- sqrt(rowSums((t(cResults[["overlaps"]]) -
        reprotable[i, ])^2) /
        (nrow(cResults[["overlaps"]]) - 1))
    ztable <- (reprotable - reprotable.P) / reprotable.sd
    sel <- which(ztable == max(ztable[is.finite(ztable)]), arr.ind = TRUE)
    if (length(sel) > 2) {
        sel <- sel[1, ]
    }
    if (sel[1] < nrow(reprotable)) {
        a1 <- as.numeric(row.names(reprotable)[sel[1]])
        a2 <- 1
    }
    if (sel[1] == nrow(reprotable)) {
        a1 <- 1
        a2 <- 0
    }
    k <- as.numeric(colnames(reprotable)[sel[2]])
    R <- reprotable[sel[1], sel[2]]
    Z <- ztable[sel[1], sel[2]]
    rm(reprotable, D, S, overlaps, overlaps.P, cResults, reprotable.P, 
                                                            reprotable.sd)
    gc()
    return(list(a1 = a1, a2 = a2, k = k, R = R, Z = Z, ztable = ztable))
}



# /*
#  * Adapted from ROTS Bioconductor (3.19, 2024) files NeedForSpeed1.cpp 
#  * and NeedForSpeed2.cpp
#  * Copyright on original version by: F. Seyednasrollah, T. Suomi,
#  * L.L. Elo (2024)
#  * Copyright on modifications (C) 2024-2025 A. Youssef
#  * Modifications: ROTS::NeedForSpeed1 and ROTS::NeedForSpeed2 were
#  * translated from C++ to R (LimROTS::.calOverlaps and 
#  * LimROTS::.calOverlaps_slr).
#  */
#
#' Calculate Overlaps Between Observed and Permuted Data
#'
#' This function calculates the overlap between observed and permuted data for
#' two sets of comparisons. It computes the ratio of overlap between pairs of
#' vectors (res1/res2 and pres1/pres2) after sorting the values.
#'
#' @param D Numeric vector. Observed data values (e.g., differences).
#' @param S Numeric vector. Standard errors or related values associated with
#' the observed data.
#' @param pD Numeric vector. Permuted data values (e.g., differences).
#' @param pS Numeric vector. Standard errors or related values associated with
#' the permuted data.
#' @param nrow Integer. Number of rows in each block of data.
#' @param N Integer vector. Number of top values to consider for overlap
#' calculation.
#' @param N_len Integer. Length of the `N` vector.
#' @param ssq Numeric. A small constant added to standard errors for stability.
#' @param niter Integer. Number of bootstrap samples or resampling iterations.
#' @param overlaps Numeric matrix. Matrix to store overlap results for observed
#' data.
#' @param overlaps_P Numeric matrix. Matrix to store overlap results for
#' permuted data.
#'
#' @details
#' The function calculates overlaps for two sets of comparisons: one for
#' observed data (res1/res2) and one for permuted data (pres1/pres2).For each
#' bootstrap sample, the function orders the two vectors being compared, then
#' calculates the proportion of overlap for the top `N` values.
#'
#' @return A list containing two matrices: \code{overlaps} for observed data and
#'  \code{overlaps_P} for permuted data.
#'
#'

.calOverlaps <- function(D, S, pD, pS, nrow, N, N_len, ssq, niter, overlaps,
    overlaps_P) {
    sort2_1R <- function(a, b) {
        order_a <- order(a, b, decreasing = TRUE)
        a <- a[order_a]
        b <- b[order_a]
        list(a = a, b = b)
    }
    idx_b <- seq_len(niter)
    idx_offset <- niter
    D <- abs(D)
    for (b in idx_b) {
        res1 <-
            abs(D[((b - 1) * nrow + 1):(b * nrow)] /
                (S[((b - 1) * nrow + 1):(b * nrow)] + ssq))
        res2 <-
            abs(D[((b + idx_offset - 1) * nrow + 1):((b + idx_offset) * nrow)] /
                (S[((b + idx_offset - 1) *
                    nrow + 1):((b + idx_offset) * nrow)] + ssq))
        pres1 <-
            abs(pD[((b - 1) * nrow + 1):(b * nrow)] /
                (pS[((b - 1) * nrow + 1):(b * nrow)] + ssq))
        pres2 <-
            abs(pD[((b + idx_offset - 1) * nrow + 1):((b + idx_offset) * nrow)]
            / (pS[((b + idx_offset - 1) * nrow + 1):((b + idx_offset) *
                    nrow)] + ssq))

        sorted_res <- sort2_1R(res1, res2)
        res1 <- sorted_res$a
        res2 <- sorted_res$b
        r3_res <- sort(res2, decreasing = TRUE)
        for (i in seq_len(N_len)) {
            N_i <- N[i]
            sum_overlap <- sum(res2[seq_len(N_i)] >= r3_res[N_i])
            overlaps[b, i] <- sum_overlap / N_i
        }
        sorted_pres <- sort2_1R(pres1, pres2)
        pres1 <- sorted_pres$a
        pres2 <- sorted_pres$b
        r3_pres <- sort(pres2, decreasing = TRUE)
        for (i in seq_len(N_len)) {
            N_i <- N[i]
            sum_overlap <- sum(pres2[seq_len(N_i)] >= r3_pres[N_i])
            overlaps_P[b, i] <- sum_overlap / N_i
        }
    }
    return(list(overlaps = overlaps, overlaps_P = overlaps_P))
}


# /*
#  * Adapted from ROTS Bioconductor (3.19, 2024) files NeedForSpeed1.cpp 
#  * and NeedForSpeed2.cpp
#  * Copyright on original version by: F. Seyednasrollah, T. Suomi,
#  * L.L. Elo (2024)
#  * Copyright on modifications (C) 2024-2025 A. Youssef
#  * Modifications: ROTS::NeedForSpeed1 and ROTS::NeedForSpeed2 were
#  * translated from C++ to R (LimROTS::.calOverlaps and 
#  * LimROTS::.calOverlaps_slr).
#  */
#
#' Calculate Overlaps for Single-Label Replicates (SLR)
#'
#' This function computes the overlap between two sets of observed and permuted
#' 'values for single-label replicates (SLR).
#' It calculates the proportion of overlap between pairs of vectors (res1/res2
#' and pres1/pres2) after sorting them.
#'
#' @param D Numeric vector. Observed data values (e.g., differences).
#' @param pD Numeric vector. Permuted data values.
#' @param nrow Integer. Number of rows in each block of data.
#' @param N Integer vector. Number of top values to consider for overlap
#' calculation.
#' @param N_len Integer. Length of the `N` vector.
#' @param niter Integer. Number of bootstrap samples or resampling iterations.
#' @param overlaps Numeric matrix. Matrix to store overlap results for
#' observed data.
#' @param overlaps_P Numeric matrix. Matrix to store overlap results for
#' permuted data.
#'
#' @details
#' The function calculates the overlap for two sets of comparisons: one for
#'  observed data (`res1`/`res2`) and one for permuted data (`pres1`/`pres2`).
#' For each bootstrap sample, the function orders the two vectors being
#' compared, then computes the proportion of overlap for the top `N` values.
#'
#' @return A list containing two matrices: \code{overlaps} for observed data and
#' \code{overlaps_P} for permuted data.
#'
#'

.calOverlaps_slr <- function(D, pD, nrow, N, N_len, niter, 
                                            overlaps, overlaps_P) {
    sort2_1R <- function(a, b) {
        order_a <- order(a, b, decreasing = TRUE)
        a <- a[order_a]
        b <- b[order_a]
        list(a = a, b = b)
    }
    idx_b <- seq_len(niter)
    idx_offset <- niter
    D <- abs(D)
    for (b in idx_b) {
        res1 <- abs(D[((b - 1) * nrow + 1):(b * nrow)])
        res2 <-
            abs(D[((b + idx_offset - 1) * nrow + 1):((b + idx_offset) * nrow)])
        pres1 <- abs(pD[((b - 1) * nrow + 1):(b * nrow)])
        pres2 <-
            abs(pD[((b + idx_offset - 1) * nrow + 1):((b + idx_offset) * nrow)])
        sorted_res <- sort2_1R(res1, res2)
        res1 <- sorted_res$a
        res2 <- sorted_res$b
        r3_res <- sort(res2, decreasing = TRUE)
        for (i in seq_len(N_len)) {
            N_i <- N[i]
            sum_overlap <- sum(res2[seq_len(N_i)] >= r3_res[N_i])
            overlaps[b, i] <- sum_overlap / N_i
        }
        sorted_pres <- sort2_1R(pres1, pres2)
        pres1 <- sorted_pres$a
        pres2 <- sorted_pres$b
        r3_pres <- sort(pres2, decreasing = TRUE)
        for (i in seq_len(N_len)) {
            N_i <- N[i]
            sum_overlap <- sum(pres2[seq_len(N_i)] >= r3_pres[N_i])
            overlaps_P[b, i] <- sum_overlap / N_i
        }
    }
    return(list(overlaps = overlaps, overlaps_P = overlaps_P))
}
