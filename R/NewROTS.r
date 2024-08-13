#' ROTS_with_covariates2
#' 
#' Perform reproducibility-optimized test statistic (ROTS) analysis considering covariates.
#'
#' @param data A numeric data matrix or an ExpressionSet instance, in which rows correspond to genes and columns correspond to samples.
#' @param groups A vector indicating group labels for each sample.
#' @param B An integer specifying the number of bootstrap and permutation resamplings (default 1000).
#' @param K An integer indicating the largest top list size considered. If no value is given, 1/4 of the features are used.
#' @param paired Logical, indicating whether the data is paired (default = FALSE).
#' @param seed Random seed for reproducibility (optional).
#' @param a1 ROTS parameters (Non-negative parameters.) (optional).
#' @param a2 ROTS parameters (Non-negative parameters.) (optional).
#' @param log A logical (deafult TRUE) indicating whether input data is log2 scaled. This information is only used to calculate log fold change.
#' @param progress Logical, display progress bar (default = FALSE).
#' @param verbose Logical, display additional messages (default = TRUE).
#' @param time Survival time data (optional).
#' @param event Event data for survival analysis (optional).
#' @param covariates Data frame containing covariates (Confounding Variables) (optional).
#' @param cluster Cluster object for parallel processing (optional).
#' @param group.name Column name in covariates indicating group labels (optional).
#' @param formula.str Formula string for covariate adjustment (optional).
#' @param trend Logical, adjust for trend in covariate analysis (default = TRUE).
#' @param robust Logical, use robust fitting in covariate analysis (default = TRUE).
#'
#' @return A list containing the ROTS analysis results.
#' @export
#' @import Biobase
#' @import ROTS
#' @import methods
#' @import stats
#' @import utils
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel stopCluster clusterExport
#' @importFrom dplyr %>%
#' @importFrom stringr str_detect
#' @importFrom doParallel registerDoParallel




ROTS_with_covariates2 <- function (data, groups, B = 1000, K = NULL, paired = FALSE, 
                                   seed = NULL, a1 = NULL, a2 = NULL, log = TRUE, progress = FALSE, 
                                   verbose = TRUE, time = NULL, event = NULL, 
                                   covariates = NULL, cluster = NULL ,
                                   group.name = NULL , formula.str = NULL, trend=TRUE, robust=TRUE) 
{
  if (is(data, "ExpressionSet")) 
    data <- Biobase::exprs(data)
  if (!is.null(seed)) 
    set.seed(seed, kind = "default")
  if (!is.null(time)) {
    groups <- time
    if (length(time) != length(event)) {
      stop(paste("Number of survival times and events do not match."))
    }
  }
  if (length(groups) != ncol(data)) {
    stop(paste("Number of samples in the data does not match the groups."))
  }
  if (is.null(rownames(data))) 
    rownames(data) <- 1:nrow(data)
  ssq <- c((0:20)/100, (11:50)/50, (6:25)/5)
  N <- c((1:20) * 5, (11:50) * 10, (21:40) * 25, (11:1000) * 
           100)
  if (is.null(K)) {
    K <- floor(nrow(data)/4)
    if (verbose) 
      message(paste("No top list size K given, using", 
                    K))
  }
  K <- min(K, nrow(data))
  N <- N[N < K]
  if (inherits(groups, "character")) {
    groups <- factor(groups)
    groups.levels <- levels(groups)
    groups <- as.numeric(groups)
  }
  else {
    groups.levels <- NULL
  }
  
  if (!is.null(time)) {
    event <- event[order(groups)]
  }
  
  #### Added part for sorting in normal rots
  
  sort.df <- data.frame(sample.id = colnames(data), groups = groups)
  
  sort.df <- sort.df[ order( sort.df$groups ), ]
  data <- data[,sort.df$sample.id]
  groups <- sort.df$groups 
  
  
  if(!is.null(covariates)){
  covariates <- covariates[order(covariates[,group.name]),]
  data <- data[, covariates$sample.id]
  covariates[, group.name] <- factor(covariates[, group.name])
  }
  
  if (is.null(time)) {
    for (i in unique(groups)) {
      if (any(rowSums(is.na(data[, which(groups == i)])) >= 
              length(which(groups == i)) - 1)) {
        if (is.null(groups.levels)) {
          target <- i
        }
        else {
          target <- groups.levels[i]
        }
        stop(paste("The data matrix of group", target, 
                   "contains rows with less than two non-missing values, please remove these rows."))
      }
    }
  }
  cl <- groups + (1 - min(groups))

  
  if (paired) {
    for (i in unique(cl)[-1]) {
      if (length(which(cl == 1)) != length(which(cl == 
                                                 i))) 
        stop("Uneven number of samples for paired test.")
    }
  }
  if (length(unique(cl)) == 2 && is.null(time)) {
    if (log) {
      logfc <- rowMeans(data[, which(cl == 1)], na.rm = TRUE) - 
        rowMeans(data[, which(cl == 2)], na.rm = TRUE)
    }
    else {
      logfc <- rowMeans(log2(data[, which(cl == 1)] + 
                               1), na.rm = TRUE) - rowMeans(log2(data[, which(cl == 
                                                                                2)] + 1), na.rm = TRUE)
    }
  }
  else {
    logfc <- rep(NA, nrow(data))
  }
  if (verbose) 
    message("Bootstrapping samples")
  
  
  
  
  if(!is.null(covariates)){
    samples <- bootstrapSamples.limRots(data, 2 * B, cl, paired, covariates, group.name )
  }else{
    samples <- bootstrapSamples(data, 2 * B, cl, paired)
  }
  
  pSamples <- permutatedSamples(data, nrow(samples), cl)
  
  
  
  D <- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))
  S <- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))
  pD <- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))
  pS <- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))
  
  pb <- txtProgressBar(min = 0, max = 100, style = 3)
  
  registerDoParallel(cluster)
  clusterExport(cluster, varlist = c( "pb", "samples" , "pSamples" , "D", "data",
                                     "S" , "pD" , "pS", "time", "formula.str", "group.name" ,
                                     "cl", "event", "covariates", "testStatistic",
                                     "testStatistic_with_covariates" , "testStatistic.surv", 
                                     "testStatistic_with_covariates_permutating",
                                     "a1" , "a2", "trend", "robust" )  ,envir = environment())
  
  
  
  if (progress) {
    setTxtProgressBar(pb, 50)
  }
  results_list <- foreach(i = seq_len(nrow(samples)), .combine = "c",
                          .packages = c("utils", "dplyr" , "stringr")) %dopar% {
                            
                            samples.R <- split(samples[i, ], cl)
                            pSamples.R <- split(pSamples[i, ], cl)
                            
                            # Initialize placeholders for results
                            d_result <- s_result <- pd_result <- ps_result <- NULL
                            
                            # Compute D and S if conditions are met
                            if (is.null(a1) | is.null(a2)) {
                              if (!is.null(time)) {
                                fit <- testStatistic.surv(lapply(samples.R, function(x) data[, x]), cl, event)
                                
                              }else if(!is.null(covariates)){
                                
                                fit <- testStatistic_with_covariates(paired = paired, 
                                                                     data = lapply(samples.R, function(x) data[, x]),
                                                                     group.name = group.name, covariates = covariates, 
                                                                     formula.str = formula.str,
                                                                     trend=TRUE, robust=TRUE)
                              }else{
                                
                                fit <- testStatistic(paired, lapply(samples.R, 
                                                                              function(x) data[, x]))
                              }
                              d_result <- fit$d
                              s_result <- fit$s
                              
                              df1 <- data.frame(d_result = d_result , s_result = s_result)
                            }
                            
                            # Compute pD and pS
                            if (!is.null(time)) {
                              pFit <- testStatistic.surv(lapply(pSamples.R, function(x) data[, x]), cl, event)
                            }else if(!is.null(covariates)){
                              
                              pFit <- testStatistic_with_covariates_permutating(paired = paired, 
                                                                    data = lapply(pSamples.R, function(x) data[, x]),
                                                                    group.name = group.name, 
                                                                    covariates = covariates, 
                                                                    formula.str = formula.str,
                                                                    trend=TRUE, robust=TRUE)
                            }else{
                              
                              pFit <- testStatistic(paired, lapply(pSamples.R, 
                                                                   function(x) data[, x]))
                            }
                            pd_result <- pFit$d
                            ps_result <- pFit$s
                            
                            df2 <- data.frame(pd_result = pd_result , ps_result = ps_result)
                            
                            
                            # Return results for this iteration as a data frame
                            list(ds = df1, pdps = df2)
                          }

  stopCluster(cluster)
  
    if (progress) {
    setTxtProgressBar(pb, 80)
  }
  
  
  j <-  0
  q <-  0
  # Populate matrices D, S, pD, pS from results
  for (i in seq_along(results_list)) {
  
    
    if (names(results_list)[i] == "ds"){
      j <- j + 1
      D[, j] <- results_list[[i]]$d_result
      S[, j] <- results_list[[i]]$s_result
    }else{
      q <- q + 1
      pD[, q] <- results_list[[i]]$pd_result
      pS[, q] <- results_list[[i]]$ps_result
    }

  }

  
  if (progress) {
    setTxtProgressBar(pb, 100)
  }
  

  if (progress) 
    close(pb)
  rm(samples, pSamples)
  gc()
  if (is.null(a1) | is.null(a2)) {
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
      cResults = calculateOverlaps1(D, S, pD, pS, nrow(D), 
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
    cResults = calculateOverlaps2(D, pD, nrow(D), as.integer(N), 
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
    if (!is.null(time)) {
      fit <- testStatistic.surv(lapply(split(1:length(cl), 
                                             cl), function(x) data[, x]), cl, event)
    }else if(!is.null(covariates)){
      fit <- testStatistic_with_covariates_permutating(paired = paired, data = lapply(split(1:length(cl), 
                                                                                cl), function(x) data[, x]),
                                           group.name = group.name , covariates = covariates , 
                                           formula.str = formula.str,
                                           trend=TRUE, robust=TRUE)
    }else{
      fit <- testStatistic(paired, lapply(split(1:length(cl), 
                                                cl), function(x) data[, x]))
    }
    d <- fit$d/(a1 + a2 * fit$s)
    pD <- pD/(a1 + a2 * pS)
    rm(pS)
    gc()
    if (verbose) 
      message("Calculating p-values")
    p <- calculateP(d, pD)
    if (verbose) 
      message("Calculating FDR")
    FDR <- calculateFDR(d, pD, progress)
    corrected.logfc <- fit$corrected.logfc
    rm(pD)
    gc()
    ROTS.output <- list(data = data, B = B, d = d, logfc = logfc,
                        pvalue = p, FDR = FDR, a1 = a1, a2 = a2, k = k, 
                        R = R, Z = Z, ztable = ztable, cl = cl , corrected.logfc = corrected.logfc)
  }
  else {
    if (!is.null(time)) {
      fit <- testStatistic.surv(lapply(split(1:length(cl), 
                                             cl), function(x) data[, x]), cl, event)
    }
    else {
      fit <- testStatistic_with_covariates_permutating(paired = paired, data = lapply(split(1:length(cl), 
                                                                                cl), function(x) data[, x]),
                                           group.name = group.name , covariates = covariates, 
                                           formula.str = formula.str,
                                           trend=TRUE, robust=TRUE)
    }
    d <- fit$d/(a1 + a2 * fit$s)
    if (verbose) 
      message("Calculating p-values")
    p <- calculateP(d, pD/(a1 + a2 * pS))
    if (verbose) 
      message("Calculating FDR")
    FDR <- calculateFDR(d, pD/(a1 + a2 * pS), progress)
    corrected.logfc <- fit$corrected.logfc
    ROTS.output <- list(data = data, B = B, d = d, logfc = logfc , 
                        pvalue = p, FDR = FDR, a1 = a1, a2 = a2, k = NULL, 
                        R = NULL, Z = NULL, cl = cl , corrected.logfc = corrected.logfc)
  }
  class(ROTS.output) <- "ROTS"
  return(ROTS.output)
}
