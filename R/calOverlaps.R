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
                        overlaps_P, cluster = NULL) {
    # Helper function to sort paired vectors
    sort2_1R <- function(a, b) {
        order_a <- order(a, b, decreasing = TRUE)
        list(a = a[order_a], b = b[order_a])
    }
    
    # Create a parallel backend if not provided
    if (is.null(cluster)) {
        cluster <- if (isWindows()) {
            SnowParam(workers = 2)
        } else {
            MulticoreParam(workers = 2)
        }
        message("Using ", class(cluster)[1], " with two workers.")
    }
    
    # Prepare data for parallel processing
    idx_b <- seq_len(niter)
    idx_offset <- niter
    D <- abs(D)
    
    # Parallelized loop
    results <- bplapply(idx_b, function(b) {
        # Compute indices for the current iteration
        idx1 <- ((b - 1) * nrow + 1):(b * nrow)
        idx2 <- ((b + idx_offset - 1) * nrow + 1):((b + idx_offset) * nrow)
        
        # Compute observed ratios
        res1 <- abs(D[idx1] / (S[idx1] + ssq))
        res2 <- abs(D[idx2] / (S[idx2] + ssq))
        
        # Compute permuted ratios
        pres1 <- abs(pD[idx1] / (pS[idx1] + ssq))
        pres2 <- abs(pD[idx2] / (pS[idx2] + ssq))
        
        # Sort results
        sorted_res <- sort2_1R(res1, res2)
        r3_res <- sort(sorted_res$b, decreasing = TRUE)
        
        sorted_pres <- sort2_1R(pres1, pres2)
        r3_pres <- sort(sorted_pres$b, decreasing = TRUE)
        
        # Calculate overlaps for observed and permuted data
        overlaps_b <- numeric(N_len)
        overlaps_P_b <- numeric(N_len)
        for (i in seq_len(N_len)) {
            N_i <- N[i]
            overlaps_b[i] <- sum(sorted_res$b[seq_len(N_i)] >= r3_res[N_i]) / N_i
            overlaps_P_b[i] <- sum(sorted_pres$b[seq_len(N_i)] >= r3_pres[N_i]) / N_i
        }
        
        # Return results for this iteration
        list(overlaps = overlaps_b, overlaps_P = overlaps_P_b)
    }, BPPARAM = cluster)
    
    # Combine results
    for (b in seq_len(niter)) {
        overlaps[b, ] <- results[[b]]$overlaps
        overlaps_P[b, ] <- results[[b]]$overlaps_P
    }
    
    return(list(overlaps = overlaps, overlaps_P = overlaps_P))
}
