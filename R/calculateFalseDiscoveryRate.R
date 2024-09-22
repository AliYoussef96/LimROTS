#' Calculate False Discovery Rate (FDR) Using Permuted Values (Optimized)
#'
#' This function calculates the false discovery rate (FDR) by comparing observed values to permuted values.
#' The function sorts observed values, compares them against permuted data, and computes FDR using the median of permutation results.
#'
#' @param observedValues Numeric vector. The observed test statistics or values to be evaluated for significance.
#' @param permutedValues Numeric matrix. The permuted test statistics or values, with rows corresponding to the same values as in `observedValues` and columns representing different permutations.
#' @param showProgress Logical. If `TRUE`, a progress bar will be shown during the computation.
#'
#' @return A numeric vector of the same length as `observedValues`, containing the estimated FDR for each observed value.
#' @importFrom stats median
#'
#' @examples
#' # example code
#'
#' observedValues <- c(2.5, 1.8, 3.1, 0.7, 2.9)
#' set.seed(123)
#' permutedValues <- matrix(rnorm(5 * 5, mean = 2, sd = 1), nrow = 5)
#' fdr <- calculateFalseDiscoveryRate(observedValues, permutedValues, showProgress = FALSE)
#' print(fdr)
#'
#' @export
#'
#'

calculateFalseDiscoveryRate <- function(observedValues,
                                        permutedValues,
                                        showProgress = FALSE) {
  observedAbs <- abs(observedValues)
  permutedAbs <- abs(permutedValues)
  ord <- order(observedAbs, decreasing = TRUE, na.last = TRUE)
  a <- observedAbs[ord]
  numPermutations <- ncol(permutedValues)
  FDRmatrix <- matrix(NA, nrow = length(a), ncol = numPermutations)
  if (showProgress) {
    progressBar <- txtProgressBar(min = 0,
                                  max = numPermutations,
                                  style = 3)
  }
  for (i in seq_len(numPermutations)) {
    a.rand <- sort(permutedAbs[, i], decreasing = TRUE, na.last = TRUE)
    n.bigger <- countLargerThan(a, a.rand)
    FDRmatrix[, i] <- n.bigger / seq_along(a)
    if (showProgress) {
      setTxtProgressBar(progressBar, i)
    }
  }
  if (showProgress) {
    close(progressBar)
  }
  falseDiscoveryRate <- apply(FDRmatrix, 1, median)
  falseDiscoveryRate[falseDiscoveryRate > 1] <- 1
  for (i in length(falseDiscoveryRate):1) {
    falseDiscoveryRate[i] <- min(falseDiscoveryRate[i:length(falseDiscoveryRate)])
  }
  return(falseDiscoveryRate)
}

#' Count Larger Permuted Values (Optimized)
#'
#' This helper function compares observed values against permuted values and counts the number of permuted values that are greater than or equal to each observed value.
#'
#' @param x Numeric vector. The observed values.
#' @param y Numeric vector. The permuted values to compare against the observed values.
#'
#' @return A numeric vector containing the counts of permuted values greater than or equal to the corresponding observed values.
#'
#'

countLargerThan <- function(x, y) {
  n <- length(x)
  counts <- numeric(n)
  j <- 1  # Index for y
  for (i in seq_len(n)) {
    while (j <= length(y) && y[j] >= x[i]) {
      j <- j + 1
    }
    counts[i] <- j - 1
  }
  return(counts)
}
