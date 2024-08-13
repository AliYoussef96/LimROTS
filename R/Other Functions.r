#' bootstrapSamples
#' 
#' Generates bootstrap samples for given data, optionally accounting for paired design.
#'
#' @param data A matrix or data frame of the data to be bootstrapped.
#' @param B Number of bootstrap samples to generate.
#' @param labels A vector of labels indicating group membership for each sample.
#' @param paired Logical, indicating whether the data is paired (default = FALSE).
#'
#' @return A matrix of bootstrap samples.
#' @export
#' @importFrom stats sample


bootstrapSamples <- function (data, B, labels, paired) 
{
  samples <- matrix(nrow = B, ncol = length(labels))
  for (i in 1:B) {
    for (label in unique(labels)) {
      pos <- which(labels == label)
      samples[i, pos] <- sample(pos, length(pos), replace = TRUE)
    }
  }
  if (paired) {
    for (i in 1:B) {
      for (label in unique(labels)[-1]) {
        pos <- which(labels == label)
        samples[i, pos] <- samples[i, which(labels == 
                                              1)] + pos[1] - 1
      }
    }
  }
  return(samples)
}

#' permutatedSamples
#' 
#' Generates permuted samples for given data.
#'
#' @param data A matrix or data frame of the data to be permuted.
#' @param B Number of permuted samples to generate.
#' @param cl Cluster object for parallel processing (optional).
#'
#' @return A matrix of permuted samples.
#' @export
#' @importFrom stats sample



permutatedSamples <- function (data, B, cl) 
{
  samples <- matrix(nrow = B, ncol = ncol(data))
  for (i in seq_len(B)) {
    samples[i, ] <- sample(seq_len(ncol(data)))
  }
  return(samples)
}

#' calculateOverlaps1
#' 
#' Calculates overlaps using the NeedForSpeed1 function.
#'
#' @param D A numeric vector of observed values.
#' @param S A numeric vector of observed values.
#' @param pD A numeric vector of permuted values.
#' @param pS A numeric vector of permuted values.
#' @param nrow Number of rows in the data.
#' @param N Number of permutations.
#' @param N_len Length of permutations.
#' @param ssq Numeric, a variance-related parameter.
#' @param B Number of bootstrap samples.
#' @param overlaps A numeric vector of overlaps.
#' @param overlaps.P A numeric vector of permuted overlaps.
#'
#' @return A numeric vector of calculated overlaps.
#' @export
#' @importFrom stats .Call


calculateOverlaps1 <- function (D, S, pD, pS, nrow, N, N_len, ssq, B, overlaps, overlaps.P) 
{
  overlap <- NeedForSpeed1(D, S, pD, pS, nrow, N, N_len, ssq, 
                           B, overlaps, overlaps.P)
  overlap
}



#' NeedForSpeed1
#' 
#' Calls the C++ function for speed-optimized overlap calculations.
#'
#' @param D A numeric vector of observed values.
#' @param S A numeric vector of observed values.
#' @param pD A numeric vector of permuted values.
#' @param pS A numeric vector of permuted values.
#' @param nrow Number of rows in the data.
#' @param N Number of permutations.
#' @param N_len Length of permutations.
#' @param ssq Numeric, a variance-related parameter.
#' @param B Number of bootstrap samples.
#' @param overlaps A numeric vector of overlaps.
#' @param overlaps_P A numeric vector of permuted overlaps.
#'
#' @return A numeric vector of calculated overlaps.
#' @export
#' @importFrom stats .Call


NeedForSpeed1 <- function (D, S, pD, pS, nrow, N, N_len, ssq, B, overlaps, overlaps_P) 
{
  .Call("ROTS_NeedForSpeed1", PACKAGE = "ROTS", D, S, pD, 
        pS, nrow, N, N_len, ssq, B, overlaps, overlaps_P)
}


#' calculateOverlaps2
#' 
#' Calculates overlaps using the NeedForSpeed2 function.
#'
#' @param D A numeric vector of observed values.
#' @param pD A numeric vector of permuted values.
#' @param nrow Number of rows in the data.
#' @param N Number of permutations.
#' @param N_len Length of permutations.
#' @param B Number of bootstrap samples.
#' @param overlaps A numeric vector of overlaps.
#' @param overlaps.P A numeric vector of permuted overlaps.
#'
#' @return A numeric vector of calculated overlaps.
#' @export
#' @importFrom stats .Call



calculateOverlaps2 <- function (D, pD, nrow, N, N_len, B, overlaps, overlaps.P) 
{
  overlap <- NeedForSpeed2(D, pD, nrow, N, N_len, B, overlaps, 
                           overlaps.P)
  overlap
}


#' NeedForSpeed2
#' 
#' Calls the C++ function for speed-optimized overlap calculations.
#'
#' @param D A numeric vector of observed values.
#' @param pD A numeric vector of permuted values.
#' @param nrow Number of rows in the data.
#' @param N Number of permutations.
#' @param N_len Length of permutations.
#' @param B Number of bootstrap samples.
#' @param overlaps A numeric vector of overlaps.
#' @param overlaps_P A numeric vector of permuted overlaps.
#'
#' @return A numeric vector of calculated overlaps.
#' @export
#' @importFrom stats .Call




NeedForSpeed2 <- function (D, pD, nrow, N, N_len, B, overlaps, overlaps_P) 
{
  .Call("ROTS_NeedForSpeed2", PACKAGE = "ROTS", D, pD, nrow, 
        N, N_len, B, overlaps, overlaps_P)
}


#' calculateP
#' 
#' Calculates p-values based on observed and permuted values.
#'
#' @param observed A numeric vector of observed values.
#' @param permuted A matrix of permuted values.
#'
#' @return A numeric vector of p-values.
#' @export
#' @importFrom stats order sort


calculateP <- function (observed, permuted) 
{
  observed_order <- order(abs(observed), decreasing = TRUE)
  observed <- sort(abs(observed), decreasing = TRUE)
  permuted <- sort(abs(as.vector(permuted)), decreasing = TRUE)
  p <- pvalue(observed, permuted)
  results <- vector(mode = "numeric", length = length(p))
  for (i in seq_along(results)) {
    results[observed_order[i]] <- p[i]
  }
  return(results)
}



#' pvalue
#' 
#' Calls the C++ function to compute p-values.
#'
#' @param a A numeric vector of observed values.
#' @param b A numeric vector of permuted values.
#'
#' @return A numeric vector of p-values.
#' @export
#' @importFrom stats .Call


pvalue <- function (a, b) 
{
  .Call("ROTS_pvalue", PACKAGE = "ROTS", a, b)
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
#' @export
#' @importFrom stats order apply


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
#' @importFrom stats sort match


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
#' @importFrom stats unique rowMeans rowSums




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
#' @importFrom stats rowMeans rowSums



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
#' @param labels A vector of labels indicating group membership for each sample.
#' @param paired Logical, indicating whether the data is paired (default = FALSE).
#' @param covariates Data frame containing covariates (Confounding Variables).
#' @param group.name Column name in covariates indicating group labels (optional).
#'
#' @return A matrix of bootstrap samples.
#' @export
#' @importFrom stats sample interaction prop.table




bootstrapSamples.limRots <- function (data, B, labels, paired, covariates, group.name) 
{
  samples <- matrix(nrow = B, ncol = length(labels))
  
  
  for (i in 1:B) {
    
    for (label in unique(labels)) {
      
      pos <- which(labels == label)
      
      covariates.pos <- covariates[covariates$sample.id %in% colnames(data)[pos],]
      
      ### Get Factor covariates
      covariates.factors <- c()
      for (j in 1:ncol(covariates.pos)){
        
        if(is.factor(covariates.pos[,j])){
          covariates.factors <- c(covariates.factors, colnames(covariates.pos)[j])
        }
        
      }
      
      covariates.factors <- covariates.factors[covariates.factors != group.name]
      
      # Combine gender and batch into a single factor
      
      covariates.pos$stratum <- interaction(covariates.pos[,covariates.factors])
      stratum_sizes <- table(covariates.pos$stratum)
      
      # Calculate the number of samples needed from each stratum
      stratum_samples <- round(length(pos) * prop.table(stratum_sizes))
      
      # Sample from each stratum
      sampled_indices <- unlist(lapply(names(stratum_samples), function(stratum) {
        stratum_indices <- covariates.pos[which(covariates.pos$stratum ==  stratum ),]$sample.id
        sample(stratum_indices, stratum_samples[stratum], replace = TRUE)
      }))
      
      samples[i, pos] <- sampled_indices 
    }
  }
  if (paired) {
    for (i in 1:B) {
      for (label in unique(labels)[-1]) {
        pos <- which(labels == label)
        samples[i, pos] <- samples[i, which(labels == 
                                              1)] + pos[1] - 1
      }
    }
  }
  return(samples)
}


