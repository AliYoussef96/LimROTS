#' LimROTS: An Extension of the ROTS Method with Limma Integration
#'
#' The `LimROTS` function performs robust ranking of differential expression statistics for omics data,
#' incorporating covariates from metadata and optionally integrating survival analysis, paired data, and more.
#'
#' @param data.exp A matrix where rows represent features (e.g., genes, proteins), and columns represent samples.
#'             The values should be log-transformed `log`, or a SummarizedExperiment object.
#' @param B An integer specifying the number of bootstrap iterations. Default is 1000.
#' @param K An optional integer representing the top list size for ranking. If not specified, it is set to one-fourth of the number of features.
#' @param a1 Optional numeric value used in the optimization process.
#' @param a2 Optional numeric value used in the optimization process.
#' @param log Logical, indicating whether the data is already log-transformed. Default is TRUE.
#' @param progress Logical, indicating whether to display a progress bar during bootstrap sampling. Default is FALSE.
#' @param verbose Logical, indicating whether to display messages during the function's execution. Default is TRUE.
#' @param meta.info A data frame containing sample-level metadata, where each row corresponds to a sample. It should include the grouping variable specified in `group.name`.
#' @param cluster A parallel cluster object for distributed computation, e.g., created by `makeCluster()`. Default is NULL.
#' @param group.name A string specifying the column in `meta.info` that represents the groups or conditions for comparison.
#' @param formula.str An optional formula string used when covariates are present in `meta.info` for modeling.
#' @param trend Logical, indicating whether to include trend fitting in the differential expression analysis. Default is TRUE.
#' @param robust Logical, indicating whether robust fitting should be used. Default is TRUE.
#' @param paired Logical, indicating whether the data represent paired samples. Default is FALSE.
#' @param n.ROTS Default is FALSE. If TRUE, all parameters related to LimROTS will be ignored, and the normal ROTS analysis will run.
#' @param survival To enable survival analysis, If TRUE then `meta.info` should contains time and event
#' @param seed.cl A seed should be set for randomization; if not, the default is 1234
#'
#' @return A list of class `"list"` with the following elements:
#' \item{data}{The original data matrix.}
#' \item{B}{The number of bootstrap samples used.}
#' \item{d}{Differential expression statistics for each feature.}
#' \item{logfc}{Log-fold change values between groups.}
#' \item{pvalue}{P-values computed based on the permutation samples.}
#' \item{FDR}{False discovery rate estimates.}
#' \item{a1}{Optimized parameter used in differential expression ranking.}
#' \item{a2}{Optimized parameter used in differential expression ranking.}
#' \item{k}{Top list size used for ranking.}
#' \item{corrected.logfc}{Corrected log-fold change values, if applicable.}
#' \item{q_values}{Estimated q-values using the `qvalue` package.}
#' \item{BH.pvalue}{Benjamini-Hochberg adjusted p-values.}
#'
#' @details
#' This function extends the ROTS (Reproducibility Optimized Test Statistic) method by integrating functionality from `limma` for differential expression analysis with optional covariates, survival data, and paired samples.
#' It allows the user to specify additional covariates and models complex experimental designs.
#'
#' @examples
#' # Example usage:
#'
#' data <- data.frame(matrix(rnorm(1000), nrow = 100, ncol = 10)) # Simulated data
#' meta.info <- data.frame(group = factor(rep(1:2, each = 5)), row.names = colnames(data))
#' formula.str <- "~ 0 + group"
#' result <- LimROTS(data, meta.info = meta.info, group.name = "group", formula.str = formula.str)
#'
#' @importFrom limma voom lmFit eBayes
#' @importFrom stats model.matrix formula p.adjust
#' @importFrom dplyr bind_cols
#' @importFrom parallel makeCluster clusterSetRNGStream clusterExport
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @import doRNG
#' @importFrom qvalue empPvals qvalue
#' @import utils
#' @import SummarizedExperiment
#' @export

LimROTS <- function (data.exp, B = 1000, K = NULL, a1 = NULL, a2 = NULL,
                     log = TRUE, progress = FALSE,
                                   verbose = TRUE, meta.info = NULL, cluster = NULL ,
                                  group.name = NULL , formula.str = NULL,
                     survival = FALSE, paired = FALSE,
                     n.ROTS = FALSE, seed.cl = 1234)
{


  SanityChecK.list <- SanityChecK(data.exp, B = B, K = K, a1 = a1, a2 = a2,
                                  meta.info = meta.info,
                          group.name = group.name , formula.str = formula.str,
                          survival = survival, paired = paired,
                          n.ROTS = n.ROTS, seed.cl = 1234)


  meta.info <-  SanityChecK.list$meta.info
  data <-  SanityChecK.list$data
  groups <-  SanityChecK.list$groups
  event <-  SanityChecK.list$event
  K <-  SanityChecK.list$K

  #### FC

  if (length(unique(groups)) == 2 && survival == FALSE) {
    group1_data <- data[, groups == 1]
    group2_data <- data[, groups == 2]

    if (log) {
      logfc <- rowMeans(group1_data, na.rm = TRUE) - rowMeans(group2_data, na.rm = TRUE)
    } else {
      logfc <- rowMeans(log2(group1_data + 1), na.rm = TRUE) - rowMeans(log2(group2_data + 1), na.rm = TRUE)
    }
  } else {
    logfc <- rep(NA, nrow(data))
  }


  if (verbose)
    message("Bootstrapping samples")


  if(n.ROTS == FALSE){
    if (ncol(meta.info) > 1){
    samples <- bootstrapSamples.limRots(data = data, B = 2 * B, meta.info = meta.info, group.name =  group.name )
    pSamples <- NULL
    }else{
      paired <- FALSE
      samples <- bootstrapS(2 * B, meta.info ,group.name, paired)
      pSamples <- NULL
    }
  }else{
    samples <- bootstrapS(2 * B, meta.info , group.name, paired)
    pSamples <- permutatedS(meta.info, 2 * B)
  }

  D <- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))
  S <- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))
  pD <- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))
  pS <- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))

  pb <- txtProgressBar(min = 0, max = 100, style = 3)

  if(is.null(cluster)){
    cluster <- makeCluster(2)
    registerDoParallel(cluster)
    message("No cluster found; only two cores will be used!")
  }else{
    registerDoParallel(cluster)
  }

  clusterSetRNGStream(cluster, iseed = seed.cl)
  clusterExport(cluster, varlist = c( "pb", "samples" , "pSamples" , "D", "data",
                                     "S" , "pD" , "pS", "time", "formula.str", "group.name" ,
                                     "groups", "event", "meta.info",
                                     "a1" , "a2", "trend", "robust" , "n.ROTS", "survival" )  ,envir = environment())



  if (progress) {
    setTxtProgressBar(pb, 50)
  }
  results_list <- foreach(i = seq_len(nrow(samples)), .combine = "c", .options.RNG = seed.cl,
                          .packages = c("utils", "dplyr" , "stringr", "stats" ,"LimROTS")) %dorng% {

                            samples.R <- split(samples[i, ], groups)

                            # Initialize placeholders for results
                            d_result <- s_result <- pd_result <- ps_result <- NULL

                            # Compute D and S if conditions are met
                            if (is.null(a1) | is.null(a2)) {
                              if (survival == TRUE) {

                                fit <- testStatSurvivalOptimized(lapply(samples.R, function(x) data[, x]), groups, event)

                              }else if(n.ROTS == FALSE){

                                fit <- testStatistic_with_covariates(data = lapply(samples.R, function(x) data[, x]),
                                                                     group.name = group.name, meta.info = meta.info,
                                                                     formula.str = formula.str,
                                                                     trend=trend, robust=robust)

                              }else{

                                fit <- testStatOptimized(paired, lapply(samples.R,
                                                                              function(x) data[, x]))
                              }
                              d_result <- fit$d
                              s_result <- fit$s

                              df1 <- data.frame(d_result = d_result , s_result = s_result)
                            }

                            # Compute pD and pS
                            if (survival == TRUE) {
                              pSamples.R <- split(pSamples[i, ], groups)

                              pFit <- testStatSurvivalOptimized(lapply(pSamples.R, function(x) data[, x]), groups, event)
                            }else if(n.ROTS == FALSE){

                              pFit <- testStatistic_with_covariates_permutating(data = lapply(split(1:length(groups),groups), function(x) data[, x]),
                                                                    group.name = group.name,
                                                                    meta.info = meta.info,
                                                                    formula.str = formula.str,
                                                                    trend=trend, robust=robust)
                            }else{
                              pSamples.R <- split(pSamples[i, ], groups)

                              pFit <- testStatOptimized(paired, lapply(pSamples.R,
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

  names(results_list) <- paste0(names(results_list) , seq(1,length(names(results_list))))

  j <-  0
  q <-  0
  # Populate matrices D, S, pD, pS from results
  for (i in seq_along(results_list)) {


    if (grepl("ds", names(results_list)[i], fixed = TRUE)){
      j <- j + 1
      D[, j] <- results_list[[names(results_list)[i]]]$d_result
      S[, j] <- results_list[[names(results_list)[i]]]$s_result
    }else{
      q <- q + 1
      pD[, q] <- results_list[[names(results_list)[i]]]$pd_result
      pS[, q] <- results_list[[names(results_list)[i]]]$ps_result
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

    ssq <- c((0:20)/100, (11:50)/50, (6:25)/5)
    N <- c((1:20) * 5, (11:50) * 10, (21:40) * 25, (11:1000) *
             100)
    K <- min(K, nrow(data))
    N <- N[N < K]

    optimized.parameters <- Optimizing(B, ssq, N, D, S, pD, pS,
                           verbose, progress)


    a1 <- optimized.parameters$a1
    a2 <- optimized.parameters$a2
    k <- optimized.parameters$k
    R <- optimized.parameters$R
    Z <- optimized.parameters$Z
    ztable <- optimized.parameters$ztable

    if (survival == TRUE) {
      fit <- testStatSurvivalOptimized(lapply(split(1:length(groups),
                                             groups), function(x) data[, x]), groups, event)
    }else if(n.ROTS == FALSE){
      fit <- testStatistic_with_covariates_Fit(data = lapply(split(1:length(groups),groups), function(x) data[, x]),
                                           group.name = group.name , meta.info = meta.info,
                                           formula.str = formula.str,
                                           trend=trend, robust=robust)
    }else{
      fit <- testStatOptimized(paired, lapply(split(1:length(groups),
                                                groups), function(x) data[, x]))
    }
    d <- fit$d/(a1 + a2 * fit$s)
    pD <- pD/(a1 + a2 * pS)
    rm(pS)
    gc()
    if (verbose)
      message("Calculating p-values")
    p <- empPvals(stat = d, stat0 = pD,
                  pool = TRUE)
    if (verbose)
      message("Calculating FDR")

    FDR <- calculateFalseDiscoveryRate(d, pD, progress)
    corrected.logfc <- fit$corrected.logfc
    rm(pD)
    gc()

    q_values <- qvalue(p, pi0.method = "bootstrap", lambda = seq(0.01,0.95, 0.01))
    BH.pvalue <- p.adjust(p, method = "BH")

    LimROTS.output <- list(data = data, B = B, d = d, logfc = logfc,
                        pvalue = p, FDR = FDR, a1 = a1, a2 = a2, k = k,
                        R = R, Z = Z, ztable = ztable, groups = groups , corrected.logfc = corrected.logfc,
                        q_values = q_values , BH.pvalue = BH.pvalue)
  }
  else {
    if (survival == TRUE) {
      fit <- testStatSurvivalOptimized(lapply(split(1:length(groups),
                                             groups), function(x) data[, x]), groups, event)
    }else if(n.ROTS == FALSE){
      fit <- testStatistic_with_covariates_Fit(data = lapply(split(1:length(groups),
                                               groups), function(x) data[, x]),
                                               group.name = group.name , meta.info = meta.info ,
                                               formula.str = formula.str,
                                               trend=trend, robust=robust)
    }else{
      fit <- testStatOptimized(paired,  lapply(split(1:length(groups),
                                                groups), function(x) data[, x]))
    }

    d <- fit$d/(a1 + a2 * fit$s)
    pD <- pD/(a1 + a2 * pS)
    if (verbose)
      message("Calculating p-values")
    p <- empPvals(stat = d, stat0 = pD,
                  pool = TRUE)
    if (verbose)
      message("Calculating FDR")
    FDR <- calculateFalseDiscoveryRate(d, pD, progress)
    corrected.logfc <- fit$corrected.logfc
    q_values <-  qvalue(p, pi0.method = "bootstrap", lambda = seq(0.01,0.95, 0.01))
    BH.pvalue <- p.adjust(p, method = "BH")
    LimROTS.output <- list(data = data, B = B, d = d, logfc = logfc ,
                        pvalue = p, FDR = FDR, a1 = a1, a2 = a2, k = NULL,
                        R = NULL, Z = NULL, groups = groups , corrected.logfc = corrected.logfc,
                        q_values = q_values , BH.pvalue = BH.pvalue)
  }
  return(LimROTS.output)
}
