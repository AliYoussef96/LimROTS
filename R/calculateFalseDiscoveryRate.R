
#' Calculate False Discovery Rate (FDR) Using Permuted Values
#'
#' This function calculates the false discovery rate (FDR) by comparing observed values to permuted values.
#' The function sorts observed values, compares them against permuted data, and computes FDR using the median of permutation results.
#'
#' @param observedValues Numeric vector. The observed test statistics or values to be evaluated for significance.
#' @param permutedValues Numeric matrix. The permuted test statistics or values, with rows corresponding to the same values as in `observedValues` and columns representing different permutations.
#' @param showProgress Logical. If `TRUE`, a progress bar will be shown during the computation.
#'
#' @details
#' This function computes the FDR by comparing the magnitude of observed values against permuted values in a sorted order.
#' The FDR is calculated as the median ratio of the number of permuted values greater than or equal to the observed values.
#' The results are adjusted to ensure that FDR values are non-increasing, as required for multiple testing procedures.
#'
#' @return A numeric vector of the same length as `observedValues`, containing the estimated FDR for each observed value.
#'
#' @export
#'

calculateFalseDiscoveryRate <- function(observedValues, permutedValues, showProgress) {
  observedAbs <- abs(observedValues)
  permutedAbs <- abs(permutedValues)

  ord <- order(observedAbs, decreasing = TRUE, na.last = TRUE)
  a <- observedAbs[ord]

  numPermutations <- ncol(permutedValues)
  FDRmatrix <- matrix(NA, nrow = length(a), ncol = numPermutations)

  if (showProgress) {
    progressBar <- txtProgressBar(min = 0, max = numPermutations, style = 3)
  }

  for (i in seq_len(numPermutations)) {
    a.rand <- sort(permutedAbs[, i], decreasing = TRUE, na.last = TRUE)
    n.bigger <- countLargerThan(a, a.rand)  # Use the new function
    FDRmatrix[ord, i] <- n.bigger / seq_along(a)

    if (showProgress) {
      setTxtProgressBar(progressBar, i)
    }
  }

  if (showProgress) {
    close(progressBar)
  }

  falseDiscoveryRate <- apply(FDRmatrix, 1, median)
  falseDiscoveryRate[falseDiscoveryRate > 1] <- 1
  falseDiscoveryRate[ord] <- rev(sapply(length(falseDiscoveryRate):1, function(x) {
    return(min(falseDiscoveryRate[ord][x:length(falseDiscoveryRate)]))
  }))

  return(falseDiscoveryRate)
}



#' Count Larger Permuted Values
#'
#' This helper function compares observed values against permuted values and counts the number of permuted values that are greater than or equal to each observed value.
#'
#' @param observed Numeric vector. The observed values.
#' @param permuted Numeric vector. The permuted values to compare against the observed values.
#'
#' @return A numeric vector containing the counts of permuted values greater than or equal to the corresponding observed values.
#'
#' @export


countLargerThan <- function(x, y) {
  # Sort both vectors in decreasing order
  sortedX <- sort(x, decreasing = TRUE, na.last = TRUE)
  sortedY <- sort(y, decreasing = TRUE, na.last = TRUE)

  # Create a matrix to compare each element of sortedX against sortedY
  comparisonMatrix <- outer(sortedX, sortedY, ">")

  # Count how many elements in sortedY are greater than each element in sortedX
  counts <- colSums(comparisonMatrix)

  return(counts)
}
