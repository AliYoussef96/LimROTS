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

calOverlaps <- function(D, S, pD, pS, nrow, N, N_len, ssq, niter, overlaps,
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
