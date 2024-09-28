#' @title Test Statistics for Survival Data
#' @description This function calculates mean differences and standard deviations for survival data
#' across different sample groups, considering at-risk samples at each unique event time.
#'
#' @param x A list of matrices or data frames, where each element represents a group of samples
#' (columns) with the same set of features (rows).
#' @param survivalTime A numeric vector containing the survival times for each sample.
#' @param survivalEvent A binary vector indicating the occurrence of an event (1 for event, 0 for censored)
#' for each sample.
#'
#' @return A list containing:
#' \item{d}{A numeric vector of mean differences for each feature across the groups.}
#' \item{s}{A numeric vector of standard deviations for the mean differences.}
#'
#' @details The function identifies unique times from the survival data and computes the mean difference
#' for each time point by comparing samples that experienced the event against those at risk.
#' It calculates the standard deviation based on the variation among at-risk samples.
#'
#'




testStatSurvivalOptimized <- function(x,
                                      survivalTime,
                                      survivalEvent) {
    sampleGroups <- x
    allSamples <- do.call("cbind", sampleGroups)
    uniqueTimes <- unique(survivalTime[survivalEvent == 1])

    meanDifference <- vector(mode = "numeric", length = nrow(allSamples))

    for (currentTime in uniqueTimes) {
        indicesAtRisk <- which(survivalTime >= currentTime)
        indicesEvent <- which(survivalTime == currentTime)
        eventIndices <- indicesEvent[which(survivalEvent[indicesEvent] == 1)]

        if (length(indicesAtRisk) > 1) {
            meanDifference <- meanDifference + (
                rowSums(as.data.frame(allSamples[, eventIndices]), na.rm = TRUE) -
                    length(eventIndices) * rowMeans(allSamples[, indicesAtRisk], na.rm = TRUE)
            )
        }
    }

    survivalSD <- vector(mode = "numeric", length = nrow(allSamples))

    for (currentTime in uniqueTimes) {
        indicesAtRisk <- which(survivalTime >= currentTime)
        indicesEvent <- which(survivalTime == currentTime)
        eventIndices <- indicesEvent[which(survivalEvent[indicesEvent] == 1)]

        if (length(indicesAtRisk) > 1) {
            survivalSD <- survivalSD + ((length(eventIndices) / length(indicesAtRisk)) *
                                            rowSums((
                                                allSamples[, indicesAtRisk] -
                                                    rowMeans(allSamples[, indicesAtRisk], na.rm = TRUE)
                                            ) ^ 2, na.rm = TRUE))
        }
    }

    survivalSD <- sqrt(survivalSD)

    return(list(d = meanDifference, s = survivalSD))
}
