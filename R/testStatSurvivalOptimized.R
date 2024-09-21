testStatSurvivalOptimized <- function(sampleGroups, survivalTime, survivalEvent) {
  allSamples <- do.call("cbind", sampleGroups)
  uniqueTimes <- unique(survivalTime[survivalEvent == 1])

  meanDifference <- vector(mode = "numeric", length = nrow(allSamples))

  for (currentTime in uniqueTimes) {
    indicesAtRisk <- which(survivalTime >= currentTime)
    indicesEvent <- which(survivalTime == currentTime)
    eventIndices <- indicesEvent[which(survivalEvent[indicesEvent] == 1)]

    if (length(indicesAtRisk) > 1) {
      meanDifference <- meanDifference + (rowSums(as.data.frame(allSamples[, eventIndices]), na.rm = TRUE) -
                                            length(eventIndices) * rowMeans(allSamples[, indicesAtRisk], na.rm = TRUE))
    }
  }

  survivalSD <- vector(mode = "numeric", length = nrow(allSamples))

  for (currentTime in uniqueTimes) {
    indicesAtRisk <- which(survivalTime >= currentTime)
    indicesEvent <- which(survivalTime == currentTime)
    eventIndices <- indicesEvent[which(survivalEvent[indicesEvent] == 1)]

    if (length(indicesAtRisk) > 1) {
      survivalSD <- survivalSD + ((length(eventIndices) / length(indicesAtRisk)) *
                                    rowSums((allSamples[, indicesAtRisk] -
                                               rowMeans(allSamples[, indicesAtRisk], na.rm = TRUE))^2, na.rm = TRUE))
    }
  }

  survivalSD <- sqrt(survivalSD)

  return(list(d = meanDifference, s = survivalSD))
}
