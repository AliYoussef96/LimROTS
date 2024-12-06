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
#' @param verbose Logical. If `TRUE`, progress messages will be displayed.
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
Optimizing <- function(niter, ssq, N, D, S, pD, pS, verbose) {
    if (verbose) {
        message("Optimizing a1 and a2")
    }
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
        cResults <- calOverlaps(
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
    cResults <- calOverlaps_slr(
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
