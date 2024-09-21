
bootstrapS <- function (B, meta.info, paired)
{
  groups <- meta.info[,group.name]
  bootsamples <- matrix(nrow = B, ncol = length(groups))
  for (i in seq_len(B)) {
    for (g in unique(groups)) {
      g.names <- row.names(meta.info)[which(groups == g)]
      bootsamples[i, which(groups == g)] <- sample(g.names, length(g.names), replace = TRUE)
    }
  }
  if (paired) {
    for (i in 1:B) {
        g.names1 <-  bootsamples[i, which(groups == unique(groups)[1])]

        g.names2 <- match(g.names1 , row.names(meta.info)) + length(g.names1)

        bootsamples[i, which(groups == unique(groups)[2])] <- row.names(meta.info)[g.names2]

    }
  }
  return(bootsamples)
}



permutatedS <- function (meta.info, B)
{
  persamples <- matrix(nrow = B, ncol = nrow(meta.info))
  for (i in seq_len(B)) {
    persamples[i, ] <- sample(row.names(meta.info))
  }
  return(persamples)
}





#' calculateFDR
#'
#' Calculates the False Discovery Rate (FDR) based on observed and permuted values.
#'
#' @param observed A numeric vector of observed values.
#' @param permuted A matrix of permuted values.
#' @param progress Logical, whether to display a progress bar (default = FALSE).
#'
#' @return A numeric vector of FDR values.
#'
#' @importFrom stats median
#'
#' @export


calculateFDR <- function (observed, permuted, progress)
{
  observed <- abs(observed)
  permuted <- abs(permuted)
  ord <- order(observed, decreasing = TRUE, na.last = TRUE)
  a <- observed[ord]
  A <- matrix(NA, nrow = length(a), ncol = ncol(permuted))
  if (progress)
    pb <- txtProgressBar(min = 0, max = ncol(A), style = 3)
  for (i in seq_len(ncol(A))) {
    a.rand <- sort(permuted[, i], decreasing = TRUE, na.last = TRUE)
    n.bigger <- biggerN(a, a.rand)
    A[ord, i] <- n.bigger/seq_along(a)
    if (progress)
      setTxtProgressBar(pb, i)
  }
  if (progress)
    close(pb)
  FDR <- apply(A, 1, median)
  FDR[FDR > 1] <- 1
  FDR[ord] <- rev(sapply(length(FDR):1, function(x) return(min(FDR[ord][x:length(FDR)]))))
  return(FDR)
}


#' biggerN
#'
#' Computes the number of elements in x that are bigger than elements in y.
#'
#' @param x A numeric vector.
#' @param y A numeric vector.
#'
#' @return A numeric vector of counts.
#' @export


biggerN <- function (x, y)
{
  x <- sort(x, decreasing = TRUE, na.last = TRUE)
  y <- sort(y, decreasing = TRUE, na.last = TRUE)
  a <- match(x, x)
  b <- x %in% y
  z <- sort(c(x, y), decreasing = TRUE, na.last = TRUE)
  d <- match(x, z)
  return(d - a + b)
}


#' testStatistic.surv
#'
#' Computes test statistics for survival data.
#'
#' @param samples A list of matrices containing expression data.
#' @param time A numeric vector of survival times.
#' @param event A binary vector indicating the event occurrence (1 for event, 0 for censored).
#'
#' @return A list with two components:
#'   \item{d}{A numeric vector of test statistics.}
#'   \item{s}{A numeric vector of standard errors.}
#' @export




testStatistic.surv <- function (samples, time, event)
{
  samples.all <- do.call("cbind", samples)
  t <- unique(time[event == 1])
  r <- vector(mode = "numeric", length = nrow(samples.all))
  for (k in t) {
    i <- which(time >= k)
    z <- which(time == k)
    d <- z[which(event[which(time == k)] == 1)]
    if (length(i) > 1) {
      r <- r + (rowSums(as.data.frame(samples.all[, d]),
                        na.rm = TRUE) - length(d) * rowMeans(samples.all[,
                                                                         i], na.rm = TRUE))
    }
  }
  s <- vector(mode = "numeric", length = nrow(samples.all))
  for (k in t) {
    i <- which(time >= k)
    z <- which(time == k)
    d <- z[which(event[which(time == k)] == 1)]
    if (length(i) > 1) {
      s <- s + ((length(d)/length(i)) * rowSums((samples.all[,
                                                             i] - rowMeans(samples.all[, i], na.rm = TRUE))^2,
                                                na.rm = TRUE))
    }
  }
  s <- s^0.5
  return(list(d = r, s = s))
}




#' testStatistic
#'
#' Computes test statistics for given samples, optionally accounting for paired design.
#'
#' @param paired Logical, indicating whether the data is paired (default = FALSE).
#' @param samples A list of matrices containing expression data.
#'
#' @return A list with two components:
#'   \item{d}{A numeric vector of test statistics.}
#'   \item{s}{A numeric vector of standard errors.}
#' @export



testStatistic <- function (paired, samples)
{
  if (length(samples) == 2) {
    X <- samples[[1]]
    Y <- samples[[2]]
    mX <- rowMeans(X, na.rm = TRUE)
    mY <- rowMeans(Y, na.rm = TRUE)
    sX <- rowSums((X - mX)^2, na.rm = TRUE)
    sY <- rowSums((Y - mY)^2, na.rm = TRUE)
    if (!paired) {
      nX <- rowSums(!is.na(X))
      nY <- rowSums(!is.na(Y))
      d <- mY - mX
      s <- sqrt(((sX + sY)/(nX + nY - 2)) * (1/nX + 1/nY))
      ind <- which(nY < 2 | nX < 2)
      d[ind] <- 0
      s[ind] <- 1
    }
    if (paired) {
      sXY <- rowSums((X - mX) * (Y - mY), na.rm = TRUE)
      n <- rowSums(!is.na(X * Y))
      d <- mY - mX
      s <- sqrt(((sX + sY)/(n + n - 2)) * (2/n) - 2/(n *
                                                       n - n) * sXY)
      ind <- which(n < 2)
      d[ind] <- 0
      s[ind] <- 1
    }
    return(list(d = d, s = s))
  }
  if (length(samples) > 2) {
    samples.all <- do.call("cbind", samples)
    if (!paired) {
      f <- sum(sapply(samples, ncol))/prod(sapply(samples,
                                                  ncol))
      r <- vector(mode = "numeric", length = nrow(samples.all))
      for (k in 1:length(samples)) {
        r <- r + (rowMeans(samples[[k]], na.rm = TRUE) -
                    rowMeans(samples.all, na.rm = TRUE))^2
      }
      d <- (f * r)^0.5
      f <- 1/sum(sapply(samples, ncol) - 1) * sum(1/sapply(samples,
                                                           ncol))
      s <- vector(mode = "numeric", length = nrow(samples.all))
      for (k in 1:length(samples)) {
        s <- s + colSums(apply(samples[[k]], 1, function(x) (x -
                                                               mean(x, na.rm = TRUE))^2), na.rm = TRUE)
      }
      s <- (f * s)^0.5
    }
    if (paired) {
      stop("Multiple paired groups not supported!")
    }
    return(list(d = d, s = s))
  }
}



#' bootstrapSamples.limRots
#'
#' Generates bootstrap samples for limRots with covariate adjustment.
#'
#' @param data A matrix of the data to be bootstrapped.
#' @param B Number of bootstrap samples to generate.
#' @param meta.info Data frame containing meta.info (Confounding Variables).
#' @param group.name Column name in meta.info indicating group labels (optional).
#'
#' @return A matrix of bootstrap samples.
#' @export




bootstrapSamples.limRots <- function (data, B, meta.info ,group.name)
{
  labels <- as.numeric( meta.info[,group.name] )
  samples <- matrix(nrow = B, ncol = length(labels))
  for (i in 1:B) {
    for (label in unique(labels)) {
      pos <- which(labels == label)
      meta.info.pos <- meta.info[meta.info[,group.name] == label,]
      meta.info.factors <- c()
      for (j in 1:ncol(meta.info)){
        if(is.factor(meta.info.pos[,j])){
          meta.info.factors <- c(meta.info.factors, colnames(meta.info.pos)[j])
        }
      }
      meta.info.factors <- meta.info.factors[meta.info.factors != group.name]
      meta.info.pos$stratum <- interaction(meta.info.pos[,meta.info.factors])
      stratum_sizes <- table(meta.info.pos$stratum)
      stratum_samples <- round(length(pos) * prop.table(stratum_sizes))
      sampled_indices <- unlist(lapply(names(stratum_samples), function(stratum) {
        stratum_indices <- row.names(meta.info.pos)[which(meta.info.pos$stratum ==  stratum )]
        sample(stratum_indices, stratum_samples[stratum], replace = TRUE)
      }))
      samples[i, pos] <- sampled_indices
    }
  }
  return(samples)
}


