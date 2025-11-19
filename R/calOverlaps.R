#' Compute overlaps between bootstrap and permuted features.
#'
#' @param D Numeric matrix of observed statistics/features.
#' @param S Numeric matrix of observed standard errors or relevant values.
#' @param pD Numeric matrix of permuted statistics/features.
#' @param pS Numeric matrix of permuted standard errors or relevant values.
#' @param nrow Integer, number of features per resample block.
#' @param N Integer vector of "top list sizes" for overlap.
#' @param N_len Integer, length of N.
#' @param ssq Numeric, regularization constant.
#' @param niter Numeric, total number of bootstrap/permutation rounds.
#' @param mat_obs Matrix (niter x length(N)) to be filled 
#' with overlaps(observed).
#' @param mat_perm Matrix (niter x length(N)) to be filled with 
#' overlaps (permuted).
#' @return List with mat_obs=overlaps (observed) and 
#' mat_perm=overlaps (permuted).
calOverlaps <- function(
    D, S, pD, pS, nrow, N, N_len, ssq, niter,
    mat_obs, mat_perm
) {
  get_top_overlap <- function(x, y, kvec) {
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
    norm_obs1 <- abs(D[startA:endA] / (S[startA:endA] + ssq))
    norm_obs2 <- abs(D[startB:endB] / (S[startB:endB] + ssq))
    mat_obs[b, ] <- get_top_overlap(norm_obs1, norm_obs2, N)
    norm_perm1 <- abs(pD[startA:endA] / (pS[startA:endA] + ssq))
    norm_perm2 <- abs(pD[startB:endB] / (pS[startB:endB] + ssq))
    mat_perm[b, ] <- get_top_overlap(norm_perm1, norm_perm2, N)
  }
  list(overlaps = mat_obs, overlaps_P = mat_perm)
}