#' Optimized Test Statistic Calculation
#'
#' This function calculates the test statistic (mean differences and standard deviations) for comparing two or more groups of samples, with support for paired or unpaired samples.
#'
#' @param isPaired Logical. If TRUE, calculates test statistics for paired samples. If FALSE, for unpaired samples.
#' @param x List of matrices or data frames. Each element represents a group of samples (columns) with the same set of features (rows).
#'
#' @details
#' The function supports comparison between two groups or multiple groups of samples. For two groups, it computes the mean differences and pooled standard deviations for unpaired samples, or paired standard deviations for paired samples.
#' When comparing more than two groups, it calculates the mean differences and standard deviations across all groups. For unpaired samples, a scaling factor based on sample size is used. Paired comparisons are only supported for two groups.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{d}: Mean differences between groups.
#'   \item \code{s}: Standard deviations (pooled or paired).
#' }
#'
#'




testStatOptimized <- function(isPaired, x) {
    sampleGroups <- x
    if (length(sampleGroups) == 2) {
        groupA <- sampleGroups[[1]]
        groupB <- sampleGroups[[2]]
        meanA <- rowMeans(groupA, na.rm = TRUE)
        meanB <- rowMeans(groupB, na.rm = TRUE)
        sumSqA <- rowSums((groupA - meanA) ^ 2, na.rm = TRUE)
        sumSqB <- rowSums((groupB - meanB) ^ 2, na.rm = TRUE)
        if (!isPaired) {
            nonNAcountA <- rowSums(!is.na(groupA))
            nonNAcountB <- rowSums(!is.na(groupB))
            meanDiff <- meanB - meanA
            pooledSD <- sqrt(((sumSqA + sumSqB) / (nonNAcountA + nonNAcountB - 2)) * (1 / nonNAcountA + 1 / nonNAcountB))
            insufficientSamples <- which(nonNAcountA < 2 | nonNAcountB < 2)
            meanDiff[insufficientSamples] <- 0
            pooledSD[insufficientSamples] <- 1
            return(list(d = meanDiff, s = pooledSD))
        } else {
            covarianceAB <- rowSums((groupA - meanA) * (groupB - meanB), na.rm = TRUE)
            pairedCount <- rowSums(!is.na(groupA * groupB))
            meanDiff <- meanB - meanA
            pairedSD <- sqrt(((sumSqA + sumSqB) / (2 * pairedCount - 2)) * (2 / pairedCount) - 2 / (pairedCount * (pairedCount - 1)) * covarianceAB)
            insufficientSamples <- which(pairedCount < 2)
            meanDiff[insufficientSamples] <- 0
            pairedSD[insufficientSamples] <- 1
        }
        return(list(d = meanDiff, s = pairedSD))
    } else if (length(sampleGroups) > 2) {
        allSamples <- do.call("cbind", sampleGroups)
        if (!isPaired) {
            factorScaling <- sum(sapply(sampleGroups, ncol)) / prod(sapply(sampleGroups, ncol))
            rowVariance <- rowSums(sapply(sampleGroups, function(group)
                (
                    rowMeans(group, na.rm = TRUE) - rowMeans(allSamples, na.rm = TRUE)
                ) ^ 2))
            meanDiff <- sqrt(factorScaling * rowVariance)
            scalingFactor <- 1 / sum(sapply(sampleGroups, ncol) - 1) *
                sum(1 / sapply(sampleGroups, ncol))
            totalVariance <- rowSums(sapply(sampleGroups, function(group)
                rowSums(( group - rowMeans(group, na.rm = TRUE) ) ^ 2, na.rm = TRUE)))
            standardDev <- sqrt(scalingFactor * totalVariance)
        } else {
            stop("Multiple paired groups are not supported!")
        }
        return(list(d = meanDiff, s = standardDev))
    }
}
