#' @export
#'
calculateFalseDiscoveryRate <- function(observedValues, permutedValues, showProgress) {
  observedAbs <- abs(observedValues)
  permutedAbs <- abs(permutedValues)

  orderedIndices <- order(observedAbs, decreasing = TRUE, na.last = TRUE)
  sortedObserved <- observedAbs[orderedIndices]

  numPermutations <- ncol(permutedValues)
  FDRmatrix <- matrix(NA, nrow = length(sortedObserved), ncol = numPermutations)

  if (showProgress) {
    progressBar <- txtProgressBar(min = 0, max = numPermutations, style = 3)
  }

  for (permIndex in seq_len(numPermutations)) {
    sortedPermuted <- sort(permutedAbs[, permIndex], decreasing = TRUE, na.last = TRUE)
    numLarger <- biggerN(sortedObserved, sortedPermuted)
    FDRmatrix[orderedIndices, permIndex] <- numLarger / seq_along(sortedObserved)

    if (showProgress) {
      setTxtProgressBar(progressBar, permIndex)
    }
  }

  if (showProgress) {
    close(progressBar)
  }

  falseDiscoveryRate <- apply(FDRmatrix, 1, median)
  falseDiscoveryRate[falseDiscoveryRate > 1] <- 1

  falseDiscoveryRate[orderedIndices] <- rev(sapply(length(falseDiscoveryRate):1, function(x) {
    return(min(falseDiscoveryRate[orderedIndices][x:length(falseDiscoveryRate)]))
  }))

  return(falseDiscoveryRate)
}
