#' Compute overlaps for single-label replicate statistics (SLR).
#'
#' @param D Numeric matrix of observed statistics/features.
#' @param pD Numeric matrix of permuted statistics/features.
#' @param nrow Integer, number of features per resample block.
#' @param N Integer vector of "top list sizes" for overlap.
#' @param N_len Integer, length of N.
#' @param niter Numeric, total number of bootstrap/permutation rounds.
#' @param mat_obs Matrix (niter x length(N)) to be filled with 
#' overlaps (observed).
#' @param mat_perm Matrix (niter x length(N)) to be filled with 
#' overlaps (permuted).
#' @return List with mat_obs=overlaps (observed) and 
#' mat_perm=overlaps (permuted).
calOverlaps_slr <- function(
    D, pD, nrow, N, N_len, niter,
    mat_obs, mat_perm
) {
  get_overlap <- function(x, y, kvec) {
    x_ord <- order(x, y, decreasing=TRUE)
    x_sorted <- x[x_ord]
    y_sorted <- y[x_ord]
    res <- numeric(length(kvec))
    threshold <- sort(y_sorted, decreasing=TRUE)
    for (i in seq_along(kvec)) {
      k <- kvec[i]
      res[i] <- mean(y_sorted[seq_len(k)] >= threshold[k])
    }
    res
  }
  for (b in seq_len(niter)) {
    startA <- ((b-1)*nrow + 1)
    endA   <- (b*nrow)
    startB <- ((b+niter-1)*nrow + 1)
    endB   <- ((b+niter)*nrow)
    # Observed overlap
    obsA <- abs(D[startA:endA])
    obsB <- abs(D[startB:endB])
    mat_obs[b, ] <- get_overlap(obsA, obsB, N)
    # Permuted overlap
    permA <- abs(pD[startA:endA])
    permB <- abs(pD[startB:endB])
    mat_perm[b, ] <- get_overlap(permA, permB, N)
  }
  list(overlaps = mat_obs, overlaps_P = mat_perm)
}