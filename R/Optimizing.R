Optimizing <- function(ssq, N, D, S, pD, pS,
                       verbose, progress){


  if (verbose)
    message("Optimizing parameters")
  reprotable <- matrix(nrow = length(ssq) + 1, ncol = length(N))
  colnames(reprotable) <- N
  row.names(reprotable) <- c(ssq, "slr")
  reprotable.P <- matrix(nrow = length(ssq) + 1, ncol = length(N))
  colnames(reprotable.P) <- N
  row.names(reprotable.P) <- c(ssq, "slr")
  reprotable.sd <- matrix(nrow = length(ssq) + 1, ncol = length(N))
  colnames(reprotable.sd) <- N
  row.names(reprotable.sd) <- c(ssq, "slr")
  if (progress)
    pb <- txtProgressBar(min = 0, max = length(ssq),
                         style = 3)
  for (i in 1:length(ssq)) {
    overlaps <- matrix(0, nrow = B, ncol = length(N))
    overlaps.P <- matrix(0, nrow = B, ncol = length(N))
    cResults = calOverlaps(D, S, pD, pS, nrow(D),
                           as.integer(N), length(N), ssq[i], as.integer(B),
                           overlaps, overlaps.P)
    reprotable[i, ] <- colMeans(cResults[["overlaps"]])
    reprotable.P[i, ] <- colMeans(cResults[["overlaps_P"]])
    reprotable.sd[i, ] <- sqrt(rowSums((t(cResults[["overlaps"]]) -
                                          reprotable[i, ])^2)/(nrow(cResults[["overlaps"]]) -
                                                                 1))
    if (progress)
      setTxtProgressBar(pb, i)
  }
  if (progress)
    close(pb)
  i <- length(ssq) + 1
  overlaps <- matrix(0, nrow = B, ncol = length(N))
  overlaps.P <- matrix(0, nrow = B, ncol = length(N))
  cResults = calOverlaps.slr(D, pD, nrow(D), as.integer(N),
                             length(N), as.integer(B), overlaps, overlaps.P)
  rm(D, S)
  gc()
  reprotable[i, ] <- colMeans(cResults[["overlaps"]])
  reprotable.P[i, ] <- colMeans(cResults[["overlaps_P"]])
  reprotable.sd[i, ] <- sqrt(rowSums((t(cResults[["overlaps"]]) -
                                        reprotable[i, ])^2)/(nrow(cResults[["overlaps"]]) -
                                                               1))
  rm(overlaps, overlaps.P, cResults)
  gc()
  ztable <- (reprotable - reprotable.P)/reprotable.sd
  rm(reprotable.P, reprotable.sd)
  gc()
  sel <- which(ztable == max(ztable[is.finite(ztable)]),
               arr.ind = TRUE)
  if (length(sel) > 2)
    sel <- sel[1, ]
  if (sel[1] < nrow(reprotable)) {
    a1 <- as.numeric(row.names(reprotable)[sel[1]])
    a2 <- 1
  }
  if (sel[1] == nrow(reprotable)) {
    a1 <- 1
    a2 <- 0
  }
  k <- as.numeric(colnames(reprotable)[sel[2]])
  R <- reprotable[sel[1], sel[2]]
  Z <- ztable[sel[1], sel[2]]
  rm(reprotable)
  gc()

  return(list(a1 = a1, a2 = a2,
              k = k, R = R, Z = Z))

}
