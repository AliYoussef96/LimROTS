# Copyright (C) 2025 Ali Youssef.
# This program is licensed under GPL (>= 2).
#
#' Sanity Check for Input Data and Parameters
#'
#' Validate analysis inputs and prepare data and metadata for downstream
#' LimROTS analysis.
#'
#' @param x A `SummarizedExperiment` containing the feature-by-sample data.
#' @param assay.type Character string or numeric index selecting the assay in
#' `x` to analyze. If `NULL`, the first assay is used.
#' @param niter Integer. Number of bootstrap or permutation replicates to
#' generate.
#' @param K Integer. Number of top-ranked features to use. If `NULL`, defaults
#' to one quarter of the number of features.
#' @param meta.info Data frame of sample-level metadata, or for a
#' `SummarizedExperiment`, a character vector naming the columns of `colData(x)`
#' to use.
#' @param group Character string naming the column in `meta.info` that
#' identifies each sample's group.
#' @param model Character string specifying the fitting method, such as
#' `"limma"` or `"dream"`.
#' @param formula.str Character string containing the model formula.
#' @param survival Logical. If `TRUE`, apply survival-analysis validation,
#' requiring `time` and `event` metadata and a survival formula. Defaults to
#' `FALSE`.
#' @param pairwise_contrasts Optional pairwise contrast specification retained
#' in the function interface; it is not evaluated by this validation function.
#' @param correlation_block Character string or `NULL`. Name of the metadata
#' column identifying correlated samples or blocks.
#' @param strata Character vector or `NULL`. Names of character or factor
#' metadata columns used to define sampling strata.
#' @param competing_risks Logical. If `TRUE`, validate metadata for competing
#' risks analysis; otherwise validate the standard survival-analysis event
#' groups.
#'
#' @details
#' Extracts the selected assay and metadata, checks their sample alignment and
#' requested analysis parameters, and returns the prepared objects. For
#' survival analysis, it delegates formula and event validation to
#' `.check_survival()`.
#'
#' @return A list containing the prepared `meta.info`, feature-by-sample
#' `data`, grouping column name `group`, and top-list size `K`. In survival
#' mode, the returned list also contains the processed `formula.str`.
SanityChecK <- function(x,
        assay.type,
        niter, K,
        meta.info, group, model,
        formula.str, survival, pairwise_contrasts,
        correlation_block, strata , competing_risks) {

    cli_alert_info("Validating analysis inputs and metadata.")
    ### Check if x is a SummarizedExperiment and extract data and metadata
    ### and extract the assay if specified
    Check_SExp <- .Check_SummarizedExperiment(
        x = x,
        assay.type = assay.type,
        meta.info = meta.info,
        group = group,
        survival = survival,
        K = K
    )
    data <- Check_SExp[["data"]]
    group <- Check_SExp[["group"]]
    meta.info <- Check_SExp[["meta.info"]]
    K <- Check_SExp[["K"]]
    
    .Check_meta_info(
        meta.info = meta.info,
        data = data,
        group = group,
        K = K,
        competing_risks = competing_risks,
        survival = survival
    )

    ### Check if strata in meta.info and it is character or factor
    if(!is.null(strata)){
        if(!all(strata %in% colnames(meta.info))){
            cli_abort("'strata' must be NULL or names of a column in meta.info.")
        }
        if(!all(sapply(meta.info[strata], 
        function(x) is.character(x)))){
            cli_abort("strata' must be character")
        }
        if(group %in% strata){
            cli_abort("group column should not be included in strata")
            }    
        }



    ### Check if survival analysis is requested and perform additional checks
    if(survival){
        cli_alert_info("Survival analysis mode activated. Performing additional checks.")
        survival_check <- .check_survival(meta.info = meta.info, 
        data = data, 
        K = K, 
        formula.str = formula.str,
        competing_risks = competing_risks,
        group  = group )
        return(survival_check)
      }

    if(!(group %in% colnames(meta.info))){
        cli_abort("The specified group column is not present in meta.info.")
    }

    ### Checks for the model type and random effects in the formula
    .check_formula(formula.str = formula.str, model = model, 
    correlation_block = correlation_block , meta.info = meta.info , group = group)

    ### Final checks and sorting of meta.info and data based on the group column
    meta.info <- meta.info[order(meta.info[,group]), , drop = FALSE]
    data <- data[, row.names(meta.info)]
    if(!identical(colnames(data) , row.names(meta.info)) ){
        cli_abort("After sorting, column names of data do not match row names of
             meta info in colData")
    }

    if(!inherits(meta.info[, group], "character")) {
        cli_abort("The group column in meta info from colData must be of type character.")
    }

    .Check_meta_info(
        meta.info = meta.info,
        data = data,
        group = group,
        K = K,
        competing_risks = competing_risks,
        survival = survival
    )

    cli_alert_success("Analysis input and metadata checks passed.")
    return(list(
        meta.info = meta.info,
        data = data,
        group = group,
        K = K
    ))
}


#' Extract and validate data and metadata from a SummarizedExperiment
#'
#' Select the requested assay and metadata columns and supply a default top-list
#' size when one was not provided.
#'
#' @param x A `SummarizedExperiment` containing the feature-by-sample data.
#' @param assay.type Character string or numeric index selecting the assay in
#' `x` to analyze. If `NULL`, the first assay is used.
#' @param meta.info Data frame of sample-level metadata, or for a
#' `SummarizedExperiment`, a character vector naming the columns of `colData(x)`
#' to use.
#' @param group Character string naming the column in `meta.info` that
#' identifies each sample's group. Not required for survival analysis.
#' @param survival Logical. If `TRUE`, skip group-column validation specific to
#' non-survival comparisons. Defaults to `FALSE`.
#' @param K Integer. Number of top-ranked features to use. If `NULL`, defaults
#' to one quarter of the number of features.
#'
#' @import SummarizedExperiment
#'
#' @return A list containing the selected assay as `data`, processed
#' `meta.info`, `group`, and `K`.

.Check_SummarizedExperiment <- function(x,
                                       assay.type, meta.info,
                                       group, survival, K) {

    
    if(!inherits(x, "SummarizedExperiment")){
        cli_abort("Input 'x' must be a SummarizedExperiment object.")
    }

    if (inherits(x, "SummarizedExperiment")) {
        cli_alert_info("Input is a SummarizedExperiment object.")
        if (is.null(meta.info)) {
            cli_abort("meta.info should be a vector of colData names to be used")
        } else {
            meta.info.colnames <- meta.info
            meta.info <- data.frame(
                colData(x)[, meta.info ,  drop = FALSE],
                check.names = FALSE,
                row.names = row.names(colData(x))
            )
            if (!identical(colnames(meta.info) , meta.info.colnames)) {
                cli_abort("meta.info should be a vector of colData names to be used")
            }
        }
        
        if (is.null(assay.type)) {
            assay.type <- assayNames(x)[1]
        }
        cli_alert_info("Using assay: {assay.type}")
        data <- assay(x, assay.type)
    } else {
        cli_abort("Input 'x' must be a SummarizedExperiment object.")
    }

    ### Check if K is NULL and set it to a quarter of the number of rows in the data
    if (is.null(K)) {
        K <- floor(nrow(data) / 4)
        cli_alert_info("Top-list size `K` not provided; using one quarter of the features: {K}")
    }

    if (anyDuplicated(rownames(data)) ||
        anyDuplicated(colnames(data)) ||
        anyDuplicated(rownames(meta.info))) {
        cli_abort("Feature names and sample names must be unique.")
    }

    if (anyNA(meta.info)) {
        cli_abort("meta.info contains missing values.")
    }

    cli_alert_success("Input data and metadata extracted.")
    return(list(
        data = data, group = group,
        meta.info = meta.info,
        K = K
    ))
}


#' Validate sample metadata and analysis groups
#'
#' Check that metadata rows correspond uniquely to data columns and that the
#' requested group and top-list size are valid.
#'
#' @param meta.info Data frame of sample-level metadata, with sample names as
#' row names.
#' @param data Feature-by-sample matrix with sample names as column names.
#' @param group Character string naming the column in `meta.info` that
#' identifies each sample's group.
#' @param K Integer. Number of top-ranked features to use.
#' @param survival Logical. If `TRUE`, validate event groups for survival
#' analysis; otherwise require two groups for a standard comparison.
#' @param competing_risks Logical. If `TRUE` during survival analysis, require
#' three unique event groups; otherwise require two.
#'
#' @return `TRUE` indicating that the checks passed; errors are raised for
#' invalid metadata or parameters.

.Check_meta_info <- function(meta.info, data , group, K, survival, competing_risks) {
    if (any(!row.names(meta.info) %in% colnames(data))) {
        cli_abort("rownames
            for colData should match the data colnames (samples names)")
    }
    if (any(grepl(".", colnames(data), fixed = TRUE))) {
        cli_abort("Sample names should contain no '.', please remove it if any.")
    }
    if (nrow(meta.info) != ncol(data)) {
        cli_abort("Number of samples in the data does not match the group.")
    }
    if(!identical(sort(rownames(meta.info)), sort(colnames(data)))){
        cli_abort("Sample names in the data and meta info in colData do not match.")
    }
    if (is.null(row.names(meta.info)) ||
        anyDuplicated(row.names(meta.info))) {
        cli_abort("meta.info must have unique sample names as row names.")
    }

    if (!group %in% colnames(meta.info)) {
        cli_abort(
            paste0(
                "group should be a string specifying",
                " the column in `meta.info` that",
                " represents the group or conditions",
                " for comparison."
            )
        )
    }

    if(!survival){
        if(length(unique(meta.info[, group])) != 2){
        cli_abort("LimROTS only supports two-group comparisons.")
        }
    }else if(survival){
        
        if(competing_risks == FALSE){ 
        
        if(length(unique(meta.info[, group])) != 2 ){
        cli_abort("Survival analysis requires at two unique groups in the 'event' column.")
        }
        }else if(competing_risks == TRUE){
            
            if(length(unique(meta.info[, group])) != 3){
        cli_abort("Competing risks analysis requires three unique groups in the 'event' column.")
            }
            
        }
    }

    if (anyNA(meta.info[, group])) {
        cli_abort("The group column contains missing values.")
    }
    if(K < 1 || K > nrow(data)){
        cli_abort("K must be between 1 and the number of rows in the data matrix.")
    }
    return(TRUE)
}


#' Validate survival-analysis metadata and formula
#'
#' Verify that required survival columns and formula covariates are present,
#' align metadata with the data columns, and validate the event groups.
#'
#' @param meta.info Data frame of sample-level metadata, with sample names as
#' row names.
#' @param data Feature-by-sample matrix with sample names as column names.
#' @param K Integer. Number of top-ranked features to use. If `NULL`, defaults
#' to one quarter of the number of features.
#' @param formula.str Character string containing the survival model formula.
#' It must include `Surv(time,event)~y`; additional covariates must be columns
#' in `meta.info`.
#' @param competing_risks Logical. If `TRUE`, validate three unique event
#' groups; otherwise validate two.
#' @param group Character string naming the column in `meta.info` that
#' identifies the event groups.
#'
#' @return A list containing the ordered `meta.info`, aligned `data`, top-list
#' size `K`, and `formula.str`.
.check_survival <- function(meta.info, data, K, formula.str, competing_risks, group ) {
  
  if(grepl("[*:]", formula.str)){
    cli_abort("For survival analysis, 'formula.str' should not contain interaction terms or 
             special characters like '*' or ':'.")
  }
  
  if(!str_detect( formula.str , fixed("Surv(time,event)~y"))){
    cli_abort("For survival analysis, 'formula.str' must include a 
             'Surv(time, event)~y'")
  }
  
  if(!all(c("time", "event") %in% colnames(meta.info))){
    cli_abort("For survival analysis, 'meta.info' must contain 'time' and 
             'event' columns.")
  }

 formula.str.temp <- str_remove(formula.str , fixed("Surv(time,event)~"))
 formula.str.temp <- unlist(str_split(formula.str.temp , fixed("+") ))
 formula.str.temp <- formula.str.temp[formula.str.temp != "y"]
 
  if(!all(formula.str.temp %in% colnames(meta.info))){
    cli_abort("All covariates specified in 'formula.str' must be present in 'meta.info'.")
  }

  if (is.null(K)) {
    K <- floor(nrow(data) / 4)
  }

  meta.info <- meta.info[order(meta.info[,"event"]), , drop = FALSE]
  data <- data[, row.names(meta.info)]
  if(!identical(colnames(data) , row.names(meta.info)) ){
    cli_abort("After sorting, column names of data do not match row names of
             meta info in colData")
  }

  .Check_meta_info(
    meta.info = meta.info,
    data = data,
    group = group,
    K = K,
    competing_risks = competing_risks,
    survival = TRUE
  )
  
    cli_alert_success("Survival-analysis checks passed.")

  return(list(
    meta.info = meta.info,
    data = data,
    K = K,
    formula.str = formula.str
  ))
  
}

#' Validate a model formula and correlation-block specification
#'
#' Check that formula variables are present in the metadata, that the selected
#' model supports the formula's random-effects terms, and that any requested
#' correlation-block column exists.
#'
#' @param formula.str Character string containing the model formula.
#' @param model Character string specifying `"limma"` or `"dream"`.
#' @param correlation_block Character string or `NULL`. Name of the metadata
#' column identifying correlated samples or blocks.
#' @param meta.info Data frame of sample-level metadata.
#'
#' @return `NULL` if all checks pass; errors are raised for invalid formulas or
#' correlation-block specifications.
.check_formula <- function(formula.str, model, correlation_block,
meta.info , group) {

        #### Check if correlation_block is provided and is a single character string
        if(!is.null(correlation_block)){
            if(length(correlation_block) != 1){
                cli_abort("correlation_block should be a single character string specifying the 
                column name in meta.info that defines correlation blocks.")
            }
            if (!correlation_block %in% colnames(meta.info)) {
            cli_abort(
                paste0(
                    "The correlation-block column '",
                    correlation_block,
                    "' was not found in meta.info."
                )
            )
        }
        }

        ### check the formula terms and meta.info agree
        vars <- unique(all.vars(formula(formula.str)))
        if(!all(vars %in% colnames(meta.info))){
            cli_abort("All variables in the formula must be present in meta.info.")
        }

        ### check if the first variable after '~' is the group variable
        formula_terms <- attr(terms(formula(formula.str)), "term.labels")
        first_var <- all.vars(as.formula(paste("~", formula_terms[1])))[1]
        if(first_var != group){
            cli_abort("The first variable in the formula must be the group variable: {group}")
        }

        ### Check for random effects in the formula
        random_terms <- lme4::findbars(formula(formula.str))
        if(model == "limma"){
            if (!is.null(random_terms)) { 
                cli_abort("Random terms are not supported with model `limma`. Use `dream` for random effects.")
                Please use model='dream' for random effects analysis.")
            }else if(model == "lmer"){
            if (!is.null(random_terms)) { 
                cli_abort("Random terms with model `lmer` are not supported in LimROTS_peptide.")
            }
            }}else if(model == "dream"){
                if (is.null(random_terms)) { 
                    cli_abort("`dream` requires at least one random term in `formula.str`.")
                }
                if(!is.null(correlation_block)){            
                    if (!any(vapply(random_terms, function(x) {
                    rlang::as_string(rlang::get_expr(x[[3]])) == correlation_block
                }, logical(1)))) {
                    cli_abort("`correlation_block` must match a random-effect grouping variable in the formula.")
                    }
                }
        }
}
