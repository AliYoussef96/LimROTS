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
#' @param time A numeric vector representing survival times, if survival analysis is used.
#' @param event A numeric vector indicating the survival event status (1 = event occurred, 0 = censored) corresponding to `time`.
#' @param paired Logical, indicating whether the data represent paired samples. Default is FALSE.
#' @param n.ROTS Default is FALSE. If TRUE, all parameters related to LimROTS will be ignored, and the normal ROTS analysis will run.
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
#' @import parallel
#' @import foreach
#' @import doRNG
#' @importFrom qvalue empPvals qvalue
#' @import utils
#' @import SummarizedExperiment
#' @importFrom doParallel registerDoParallel
#' @export

LimROTS <- function (data.exp, B = 1000, K = NULL, a1 = NULL, a2 = NULL,
                     log = TRUE, progress = FALSE,
                                   verbose = TRUE, meta.info = NULL, cluster = NULL ,
                                  group.name = NULL , formula.str = NULL, trend = TRUE, robust = TRUE,
                     survival = FALSE, paired = FALSE,
                     n.ROTS = FALSE, seed.cl = 1234)
{


  ### SummarizedExperiment

  if(inherits(data.exp, "SummarizedExperiment")){
    message("Data is SummarizedExperiment object")

    if(is.null(meta.info)){
      stop("meta.info should be a vector of colData names to be used")
    }else{
      meta.info.colnames <- meta.info
      meta.info <- data.frame(colData(data.exp)[,meta.info], check.names = FALSE, row.names = row.names(colData(data.exp)))

      if(length(meta.info) != length(meta.info.colnames)){
        stop("meta.info should be a vector of colData names to be used")
      }else{
        colnames(meta.info) <- meta.info.colnames
      }

    }

    if(!group.name %in% colnames(meta.info)){
      stop("group.name should be a string specifying the column in `meta.info` that represents the groups or conditions for comparison.")
    }else{
      groups <- meta.info[,group.name]
    }

    message( sprintf("Assay: %s will be used" , assayNames(data.exp)[1]) )
    data <- assay(data.exp , assayNames(data.exp)[1])

  }else{
    data <- data.exp
    }


  ### meta.info

  if(any(!row.names(meta.info) %in% colnames(data))){
    stop("rownames for meta.info should match the data colnames (samples names)")
  }


  if(any(grepl("." , colnames(data) , fixed = TRUE))){
    stop("Sample names should contains no '.', please remove it if any")
  }


  if(!is.null(meta.info) & n.ROTS == FALSE){
    if(ncol(meta.info) == 1){
      message("A meta.info table is provided with only group infomration >>> LimROTS with no covariates will be used")
      if(is.null(formula.str)){
        stop("formula.str should by provided for the model")
      }
    }else{
      message("A meta.info table is provided with covariates >>> LimROTS with covariates will be used")
      if(is.null(formula.str)){
        stop("formula.str should by provided for the model")
      }
    }
  }else{
    message("n.ROTS is TRUE >>> ROTS will be used")

  }

  ### Sort

  if (nrow(meta.info) != ncol(data)) {
    stop("Number of samples in the data does not match the groups.")
  }

  sort.df <- data.frame(sample.id = colnames(data), groups = meta.info[,group.name])
  sort.df <- sort.df[ order( sort.df$groups ), ]
  data <- data[,sort.df$sample.id]


  meta.info$temp <- row.names(meta.info)
  meta.info <- data.frame(meta.info[colnames(data),], check.rows = F, check.names = F)
  meta.info$temp <- NULL


  ### Groups


  if(!inherits(meta.info[,group.name], "character")){
    meta.info[,group.name] <- factor(meta.info[,group.name])
    groups.levels <- levels(groups)
    meta.info[,group.name] <- as.numeric(meta.info[,group.name])
  }

  groups <- meta.info[,group.name]

  if(survival  == TRUE){
    if(all(c("time", "event") %in% colnames(meta.info))){
      stop("meta.info must have two columns time and event. Also, group.name must be time")
    }
    event <- meta.info[,"event"]
    groups <- meta.info[,"time"]
  }




  if (is.null(K)) {
    K <- floor(nrow(data)/4)
    if (verbose)
      message(sprintf("No top list size K given, using %s",
                      K))
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



  if(n.ROTS == FALSE)
    {
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
                                     "cl", "event", "meta.info",
                                     "a1" , "a2", "trend", "robust" )  ,envir = environment())



  if (progress) {
    setTxtProgressBar(pb, 50)
  }
  results_list <- foreach(i = seq_len(nrow(samples)), .combine = "c", .options.RNG = seed.cl,
                          .packages = c("utils", "dplyr" , "stringr", "stats" ,"LimROTS")) %dorng% {

                            samples.R <- split(samples[i, ], cl)

                            # Initialize placeholders for results
                            d_result <- s_result <- pd_result <- ps_result <- NULL

                            # Compute D and S if conditions are met
                            if (is.null(a1) | is.null(a2)) {
                              if (!is.null(time)) {

                                fit <- testStatSurvivalOptimized(lapply(samples.R, function(x) data[, x]), cl, event)

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
                            if (!is.null(time)) {
                              pSamples.R <- split(pSamples[i, ], cl)

                              pFit <- testStatSurvivalOptimized(lapply(pSamples.R, function(x) data[, x]), cl, event)
                            }else if(n.ROTS == FALSE){

                              pFit <- testStatistic_with_covariates_permutating(data = lapply(split(1:length(cl),cl), function(x) data[, x]),
                                                                    group.name = group.name,
                                                                    meta.info = meta.info,
                                                                    formula.str = formula.str,
                                                                    trend=trend, robust=robust)
                            }else{
                              pSamples.R <- split(pSamples[i, ], cl)

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

    if (!is.null(time)) {
      fit <- testStatSurvivalOptimized(lapply(split(1:length(cl),
                                             cl), function(x) data[, x]), cl, event)
    }else if(n.ROTS == FALSE){
      fit <- testStatistic_with_covariates_Fit(data = lapply(split(1:length(cl),cl), function(x) data[, x]),
                                           group.name = group.name , meta.info = meta.info,
                                           formula.str = formula.str,
                                           trend=trend, robust=robust)
    }else{
      fit <- testStatOptimized(paired, lapply(split(1:length(cl),
                                                cl), function(x) data[, x]))
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

    ROTS.output <- list(data = data, B = B, d = d, logfc = logfc,
                        pvalue = p, FDR = FDR, a1 = a1, a2 = a2, k = k,
                        R = R, Z = Z, ztable = ztable, cl = cl , corrected.logfc = corrected.logfc,
                        q_values = q_values , BH.pvalue = BH.pvalue)
  }
  else {
    if (!is.null(time)) {
      fit <- testStatSurvivalOptimized(lapply(split(1:length(cl),
                                             cl), function(x) data[, x]), cl, event)
    }else if(n.ROTS == FALSE){
      fit <- testStatistic_with_covariates_Fit(data = lapply(split(1:length(cl),
                                               cl), function(x) data[, x]),
                                               group.name = group.name , meta.info = meta.info ,
                                               formula.str = formula.str,
                                               trend=trend, robust=robust)
    }else{
      fit <- testStatOptimized(paired,  lapply(split(1:length(cl),
                                                cl), function(x) data[, x]))
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
    ROTS.output <- list(data = data, B = B, d = d, logfc = logfc ,
                        pvalue = p, FDR = FDR, a1 = a1, a2 = a2, k = NULL,
                        R = NULL, Z = NULL, cl = cl , corrected.logfc = corrected.logfc,
                        q_values = q_values , BH.pvalue = BH.pvalue)
  }
  return(ROTS.output)
}
