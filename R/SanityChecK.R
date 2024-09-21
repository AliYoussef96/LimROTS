SanityChecK <- function(data.exp, B = 1000, K = NULL, a1 = NULL, a2 = NULL,
                       meta.info = NULL,
                       group.name = NULL , formula.str = NULL,
                       survival = FALSE, paired = FALSE,
                       n.ROTS = FALSE, seed.cl = 1234){


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
  }else{
    event = NULL
    groups = NULL
  }

  groups <- groups + (1 - min(groups))

  ### paired

  if (paired) {
    for (i in unique(groups)[-1]) {
      if (length(which(groups == 1)) != length(which(groups ==
                                                     i)))
        stop("Uneven number of samples for paired test.")
    }
  }


  if (is.null(K)) {
    K <- floor(nrow(data)/4)
    if (verbose)
      message(sprintf("No top list size K given, using %s",
                      K))
  }

  return(list(meta.info = meta.info, data = data, groups = groups, event = event, K = K))

}
