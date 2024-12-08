#' Calculate False Discovery Rate (FDR) Using Permuted Values (Adjusted)
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
#' @examples
#' observedValues <- c(2.5, 1.8, 3.1, 0.7, 2.9)
#' set.seed(123)
#' permutedValues <- matrix(rnorm(5 * 5, mean = 2, sd = 1), nrow = 5)
#' fdr <- calculateFalseDiscoveryRate(observedValues, permutedValues)
#' print(fdr)
#'

calculateFalseDiscoveryRate <- function(observedValues, permutedValues) {
    observedAbs <- abs(observedValues)
    permutedAbs <- abs(permutedValues)
    ord.obs <- order(observedAbs, decreasing = TRUE, na.last = TRUE)
    abs.obs <- observedAbs[ord.obs]
    numPermutations <- ncol(permutedValues)
    FDRmatrix <- matrix(NA, nrow = length(abs.obs), ncol = numPermutations)
    for (i in seq_len(numPermutations)) {
        a.rand <- sort(permutedAbs[, i], decreasing = TRUE, na.last = TRUE)
        bigger <- countLargerThan(abs.obs, a.rand)
        FDRmatrix[ord.obs, i] <- bigger / seq_along(abs.obs)
    }
    falseDiscoveryRate <- apply(FDRmatrix, 1, median)
    falseDiscoveryRate[falseDiscoveryRate > 1] <- 1
    falseDiscoveryRate[ord.obs] <-
        rev(vapply(
            length(falseDiscoveryRate):1,
            function(x) {
                min(falseDiscoveryRate[ord.obs][x:length(falseDiscoveryRate)])
            },
            numeric(1)
        ))
    return(falseDiscoveryRate)
}

#' Count Larger Permuted Values (Modified)
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
    # Sort observed and permuted vectors in decreasing order
    observedVec <- sort(observedVec, decreasing = TRUE, na.last = TRUE)
    permutedVec <- sort(permutedVec, decreasing = TRUE, na.last = TRUE)

    # Get the positions of elements in observed vector
    observedPos <- match(observedVec, observedVec)

    # Check if each observed element exists in permuted vector
    observedInPermuted <- observedVec %in% permutedVec

    # Combine observed and permuted into a single sorted vector
    combinedSorted <- sort(c(observedVec, permutedVec),
        decreasing = TRUE,
        na.last = TRUE
    )

    # Match observed elements to the combined vector positions
    combinedPos <- match(observedVec, combinedSorted)

    # Return the count of larger values for each observed value
    return(combinedPos - observedPos + observedInPermuted)
}
