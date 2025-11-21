#' Calculate False Discovery Rate (FDR) Using Permuted Values
#'
#' This function calculates the false discovery rate (FDR) by comparing
#' observed values to permuted values.The function sorts observed values,
#' compares them against permuted data, and computes FDR using the median of
#' permutation results.
#'
#' @param observedValues Numeric vector. The observed test statistics or values
#' to be evaluated for significance.
#' @param permutedValues Numeric matrix. The permuted test statistics or values,
#'  with rows corresponding to the same values as in `observedValues` and
#'  columns representing different permutations.
#'
#' @return A numeric vector of the same length as `observedValues`, containing
#' the estimated FDR for each observed value.
#' @importFrom stats median
#'
#'
#'
#'

calculateFalseDiscoveryRate <- function(observedValues, permutedValues) {
  obs_abs <- abs(observedValues)
  perm_abs <- abs(permutedValues)
  ord <- order(obs_abs, decreasing = TRUE, na.last = TRUE)
  obs_sorted <- obs_abs[ord]
  numPermutations <- ncol(perm_abs)
  
  perm_sorted <- apply(perm_abs, 2, function(col) {
    sort(col, decreasing = TRUE, na.last = TRUE)
  })
  
  perm_counts <- vapply(
    seq_len(numPermutations),
    function(i) {
      countLargerThan(obs_sorted, perm_sorted[, i])
    },
    numeric(length(obs_sorted))
  )
  
  ranks <- seq_along(obs_sorted)
  fdr_mat <- sweep(perm_counts, 1, ranks, FUN = "/")
  
  fdr <- apply(fdr_mat, 1, median)
  fdr[fdr > 1] <- 1
  
  fdr_rev <- rev(cummin(rev(fdr)))
  
  out <- numeric(length(fdr_rev))
  out[ord] <- fdr_rev
  return(out)
}


#' Count Larger Permuted Values
#'
#' This helper function compares observed values against permuted values and
#' counts the number of permuted values that are greater than or equal to each
#' observed value.
#'
#' @param observedVec Numeric vector. The observed values.
#' @param permutedVec Numeric vector. The permuted values to compare against
#' the observed values.
#'
#' @return A numeric vector containing the counts of permuted values greater
#' than or equal to the corresponding observed values.
#'
#'

countLargerThan <- function(observedVec, permutedVec) {
  obs <- sort(observedVec, decreasing = TRUE, na.last = TRUE)
  perm <- sort(permutedVec, decreasing = TRUE, na.last = TRUE)
  combined <- sort(c(obs, perm), decreasing = TRUE, na.last = TRUE)
  obs_pos_combined <- match(obs, combined)
  obs_pos_internal <- seq_along(obs)
  in_perm <- obs %in% perm
  return(obs_pos_combined - obs_pos_internal + in_perm)
}
