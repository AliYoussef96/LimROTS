# ----------------------------------------------------------------------------
# Peptide-level fits (Julia backend)
#
# Same statistical model as Fit_peptide.R:
#     y ~ <fixed part> + (1 | peptide_id) + (1 | sample_id), REML,
#     group contrast from the pairwise_contrasts string (asymptotic SE).
#
# Division of work
#   R     : all sampling, all data preparation, all checks, the design matrix
#           X and the contrast vector L (so direction/coding are decided here).
#   Julia : fits the models of ONE replicate, all proteins in parallel threads.
#
# Julia is only needed when LimROTS_peptide() is called: nothing here runs
# at package build/check time. The
# JuliaConnectoR package is only listed under Suggests and is only used
# (via `::`) when a Julia function is actually called.
# ----------------------------------------------------------------------------

## Session-level state (the Julia process and the loaded Julia function)
.limrots_julia <- new.env(parent = emptyenv())


#' Start Julia, install the Julia packages if needed, load limrots_peptide.jl
#'
#' Runs once per R session. Julia packages are installed in a private
#' environment under `tools::R_user_dir("LimROTS", "data")`, never in the
#' user's global Julia environment, and only after consent (interactive prompt
#' or `options(LimROTS.julia_install = TRUE)`).
#'
#' @param n_threads Number of Julia threads: a number or `"auto"`. Only effective
#' if Julia has not been started yet in this R session.
#' @return `TRUE` invisibly.
.setup_julia_peptide <- function(n_threads) {

  ## Threads must be requested BEFORE Julia starts.
  Sys.setenv(JULIA_NUM_THREADS = as.character(n_threads))

  if (isTRUE(.limrots_julia$ready)) {
    return(invisible(TRUE))
  }

  if (!requireNamespace("JuliaConnectoR", quietly = TRUE)) {
    cli_abort(c(
      "The Julia backend needs the {.pkg JuliaConnectoR} R package.",
      "i" = "Install it with {.code install.packages(\"JuliaConnectoR\")}."
    ))
  }

  if (!JuliaConnectoR::juliaSetupOk()) {
    cli_abort(c(
      "Julia was not found or could not be started.",
      "i" = "Install a recent Julia from https://julialang.org and make sure
             `julia` is on your PATH (or set the JULIA_BINDIR variable)."
    ))
  }

  ## Private Julia environment for LimROTS
  env_dir <- file.path(tools::R_user_dir("LimROTS", which = "data"),
                       "julia_env")
  dir.create(env_dir, recursive = TRUE, showWarnings = FALSE)
  env_dir <- normalizePath(env_dir, winslash = "/")

  missing_pkgs <- JuliaConnectoR::juliaEval(sprintf('begin
    import Pkg
    Pkg.activate("%s"; io = devnull)
    [p for p in ["MixedModels", "DataFrames"]
       if !(p in keys(Pkg.project().dependencies))]
  end', env_dir))

  if (length(missing_pkgs) > 0) {
    missing_pkgs <- as.character(unlist(missing_pkgs))
    allowed <- isTRUE(getOption("LimROTS.julia_install", FALSE)) ||
      (interactive() && isTRUE(utils::askYesNo(
        paste0("LimROTS needs the Julia packages ",
               paste(missing_pkgs, collapse = ", "),
               ". Install them now (one time only)?"))))
    if (!allowed) {
      cli_abort(c(
        "Missing Julia packages: {.val {missing_pkgs}}.",
        "i" = "Run interactively to be asked, or set
               {.code options(LimROTS.julia_install = TRUE)}."
      ))
    }
    cli_alert_info("Installing Julia packages (one time only)...")
    JuliaConnectoR::juliaCall("Pkg.add", missing_pkgs)
  }

  ## Load the Julia file that ships with the package
  julia_file <- system.file("julia", "limrots_peptide.jl", package = "LimROTS")
  if (!nzchar(julia_file)) {
    cli_abort("Could not find {.file inst/julia/limrots_peptide.jl}.")
  }
  JuliaConnectoR::juliaCall("include", julia_file)
  .limrots_julia$fit <- JuliaConnectoR::juliaFun("limrots_fit_proteins")

  n_active <- JuliaConnectoR::juliaEval("Threads.nthreads()")
  cli_alert_info("Julia backend ready, using {n_active} thread{?s}.")
  if (is.numeric(n_threads) && n_active != n_threads) {
    cli_alert_warning(
      "Requested {n_threads} thread{?s} but Julia has {n_active}
       (Julia was probably already running in this session)."
    )
  }

  .limrots_julia$ready <- TRUE
  invisible(TRUE)
}


#' Abort if any protein has fewer than two peptides
#'
#' With a single peptide, the number of `sample_id` levels equals the number
#' of observations and the model cannot be fitted.
#'
#' @param peptide_metadata Data frame with `peptide_id` and `protein_id`.
#' @return `TRUE` invisibly; errors otherwise.
.check_peptides_per_protein <- function(peptide_metadata) {
  n_peptides <- table(peptide_metadata[["protein_id"]])
  single <- names(n_peptides)[n_peptides < 2]
  if (length(single) > 0) {
    cli_abort(c(
      "{length(single)} protein{?s} ha{?s/ve} only one peptide.",
      "x" = "First: {.val {utils::head(single, 10)}}",
      "i" = "Remove these proteins (or add peptides) before running
             LimROTS_peptide()."
    ))
  }
  invisible(TRUE)
}


#' Build the fixed-effect design matrix and the contrast vector in R
#'
#' This reproduces what `lmer()` + `emmeans(specs = group)` + `contrast()` use:
#' treatment coding, numeric covariates at their mean, factor covariates
#' averaged with equal weights, contrast weights from `limma::makeContrasts`.
#'
#' @param covariates Sample-level data frame (rows in the same order as the
#' columns of the intensity matrix; for bootstrap it contains the duplicated
#' resampled rows).
#' @param group Name of the group column.
#' @param model_formula Full model formula (random terms are removed here).
#' @param pairwise_contrasts Single contrast string such as `"groupA - groupB"`.
#'
#' @return A list with `X` (samples x k, intercept column removed, unnamed) and
#' `L` (length k + 1, intercept first).
.peptide_design <- function(covariates, group, model_formula,
                            pairwise_contrasts) {

  if (length(pairwise_contrasts) != 1L) {
    cli_abort("Exactly one contrast string is supported in `pairwise_contrasts`.")
  }

  ## Fixed part of the formula, one-sided: ~ group + covariates
  fixed <- stats::formula(
    stats::delete.response(stats::terms(lme4::nobars(model_formula)))
  )

  ## --- Formula limitations (only simple terms are supported) ---------------
  vars <- all.vars(fixed)
  allowed_symbols <- c("~", "+", "*", ":", vars)
  unsupported <- setdiff(all.names(fixed), allowed_symbols)
  if (length(unsupported) > 0) {
    cli_abort(c(
      "The formula contains unsupported elements: {.val {unsupported}}.",
      "i" = "The Julia backend supports variables combined with
             {.code +}, {.code *} and {.code :} only (no poly(), I(), log(),
             {.code -}, {.code ^}, {.code /}, {.code %in%}, ...)."
    ))
  }
  reserved <- intersect(vars, c("y", "peptide_id", "sample_id", "protein_id"))
  if (length(reserved) > 0) {
    cli_abort("Reserved names used as covariates: {.val {reserved}}.")
  }
  if (!group %in% vars) {
    cli_abort("The group variable {.val {group}} must be in the formula.")
  }

  ## --- Covariates: characters -> factors, reject unsupported types ---------
  cov <- as.data.frame(covariates)[, vars, drop = FALSE]
  for (v in vars) {
    x <- cov[[v]]
    if (is.ordered(x) || is.logical(x)) {
      cli_abort("Covariate {.val {v}} is ordered/logical, which is not supported.")
    } else if (is.character(x)) {
      cov[[v]] <- factor(x)
    } else if (is.factor(x)) {
      cov[[v]] <- droplevels(x)
    } else if (!is.numeric(x)) {
      cli_abort("Covariate {.val {v}} must be numeric, character or factor.")
    }
  }
  if (!is.factor(cov[[group]])) {
    cli_abort("The group column {.val {group}} must be character or factor.")
  }
  group_levels <- levels(cov[[group]])

  ## --- Design matrix of the samples ---------------------------------------
  X_full <- stats::model.matrix(fixed, cov)
  if (colnames(X_full)[1] != "(Intercept)") {
    cli_abort("The model must contain an intercept.")
  }
  if (qr(X_full)$rank < ncol(X_full)) {
    cli_abort(c(
      "The fixed-effect design matrix is rank deficient.",
      "i" = "A covariate is constant, collinear, or a factor level is missing
             in this replicate."
    ))
  }

  ## --- Reference grid (same idea as emmeans) -------------------------------
  grid_values <- lapply(cov, function(v) {
    if (is.factor(v)) levels(v) else mean(v)
  })
  ref_grid <- expand.grid(grid_values, KEEP.OUT.ATTRS = FALSE,
                          stringsAsFactors = FALSE)
  for (v in vars) {
    if (is.factor(cov[[v]])) {
      ref_grid[[v]] <- factor(ref_grid[[v]], levels = levels(cov[[v]]))
    }
  }
  X_ref <- stats::model.matrix(fixed, ref_grid)

  ## Mean design row of each group level (r x n_levels)
  group_rows <- vapply(
    group_levels,
    function(g) colMeans(X_ref[ref_grid[[group]] == g, , drop = FALSE]),
    numeric(ncol(X_ref))
  )

  ## --- Contrast weights: identical call as in Fit_peptide.R ----------------
  weights <- limma::makeContrasts(
    contrasts = stringr::str_remove_all(pairwise_contrasts,
                                        stringr::fixed(group)),
    levels = group_levels
  )[, 1]

  L <- drop(group_rows %*% weights)

  list(X = unname(X_full[, -1, drop = FALSE]), L = unname(L))
}


#' Fit all proteins of one replicate in Julia (shared by the three fit types)
#'
#' @param combined_data Peptide x sample numeric matrix.
#' @param covariates.p Sample-level data frame, row names equal to the columns
#' of `combined_data`.
#' @param group Name of the group column.
#' @param model_formula Full model formula.
#' @param pairwise_contrasts Single contrast string, e.g. `"groupA - groupB"`.
#' @param peptide_metadata Data frame mapping `peptide_id` to `protein_id`.
#'
#' @return A list with the SIGNED contrast estimates `d` and standard errors
#' `s`, named by protein in order of first appearance in `peptide_metadata`.
.fit_peptide_julia_core <- function(combined_data, covariates.p, group,
                                    model_formula, pairwise_contrasts,
                                    peptide_metadata) {

  ## Same alignment checks as the R implementation
  stopifnot(identical(colnames(combined_data), row.names(covariates.p)))
  stopifnot(identical(row.names(combined_data), peptide_metadata[["peptide_id"]]))
  if (anyDuplicated(peptide_metadata[["peptide_id"]])) {
    cli_abort("peptide_metadata contains duplicated peptide_id values.")
  }
  if (anyNA(combined_data)) {
    cli_abort("The Julia backend does not support missing values yet.")
  }

  design <- .peptide_design(covariates.p, group, model_formula,
                            pairwise_contrasts)

  ## Proteins in order of first appearance (same as unique(long_data$protein_id))
  protein_names <- unique(as.character(peptide_metadata[["protein_id"]]))
  protein_index <- match(as.character(peptide_metadata[["protein_id"]]),
                         protein_names)

  Y <- unname(as.matrix(combined_data))
  storage.mode(Y) <- "double"

  res <- .limrots_julia$fit(
    Y,
    as.integer(protein_index),
    design$X,
    as.numeric(design$L),
    protein_names
  )
  if (is.null(names(res))) {
    names(res) <- c("protein", "estimate", "se")
  }

  ## Julia must return exactly the proteins we sent, in the same order
  if (!identical(as.character(unlist(res$protein)), protein_names)) {
    cli_abort("Protein names/order returned by Julia do not match the input.")
  }
  estimate <- as.numeric(res$estimate)
  se <- as.numeric(res$se)
  if (length(estimate) != length(protein_names) ||
      length(se) != length(protein_names)) {
    cli_abort("Julia returned results of the wrong length.")
  }

  list(d = setNames(estimate, protein_names),
       s = setNames(se, protein_names))
}


#' Fit the original data in Julia (drop-in for fit_original_peptide)
#'
#' @param x Peptide x sample matrix (sample names as column names).
#' @param group Name of the group column.
#' @param meta.info Sample metadata, sample names as row names.
#' @param model_formula Full model formula.
#' @param pairwise_contrasts Single contrast string.
#' @param peptide_metadata Data frame mapping `peptide_id` to `protein_id`.
#' @return List with signed `d` and `s`.
fit_original_peptide <- function(x, group, meta.info, model_formula,
                                       pairwise_contrasts, peptide_metadata) {
  stopifnot(identical(colnames(x), row.names(meta.info)))
  .setup_julia_peptide()
  .fit_peptide_julia_core(
    combined_data = x, covariates.p = meta.info, group = group,
    model_formula = model_formula, pairwise_contrasts = pairwise_contrasts,
    peptide_metadata = peptide_metadata
  )
}


#' Fit one bootstrap sample in Julia (drop-in for fit_bootstrap_peptide)
#'
#' Uses `.generate_data_for_fit_peptide()` unchanged.
#'
#' @param x List of peptide x sample matrices, one per group.
#' @param group Name of the group column.
#' @param meta.info Sample metadata, sample names as row names.
#' @param model_formula Full model formula.
#' @param pairwise_contrasts Single contrast string.
#' @param peptide_metadata Data frame mapping `peptide_id` to `protein_id`.
#' @return List with absolute `d` and `s`.
fit_bootstrap_peptide <- function(x, group, meta.info, model_formula,
                                        pairwise_contrasts, peptide_metadata) {
  data_for_fit <- .generate_data_for_fit_peptide(x = x, meta.info = meta.info)
  fit <- .fit_peptide_julia_core(
    combined_data = data_for_fit$combined_data,
    covariates.p = data_for_fit$covariates.p,
    group = group, model_formula = model_formula,
    pairwise_contrasts = pairwise_contrasts,
    peptide_metadata = peptide_metadata
  )
  list(d = abs(fit$d), s = fit$s)
}


#' Fit one permutation in Julia (drop-in for fit_permutation_peptide)
#'
#' @param x Peptide x sample matrix (sample names as column names).
#' @param group Name of the group column.
#' @param meta.info Permuted sample metadata, sample names as row names.
#' @param model_formula Full model formula.
#' @param pairwise_contrasts Single contrast string.
#' @param peptide_metadata Data frame mapping `peptide_id` to `protein_id`.
#' @return List with absolute `d` and `s`.
fit_permutation_peptide <- function(x, group, meta.info, model_formula,
                                          pairwise_contrasts, peptide_metadata) {
  stopifnot(identical(colnames(x), row.names(meta.info)))
  fit <- .fit_peptide_julia_core(
    combined_data = x, covariates.p = meta.info, group = group,
    model_formula = model_formula, pairwise_contrasts = pairwise_contrasts,
    peptide_metadata = peptide_metadata
  )
  list(d = abs(fit$d), s = fit$s)
}


#' Julia version of parallel_peptide()
#'
#' Loops over the replicates in R (sequentially) and fits every replicate in
#' Julia, where the proteins are fitted in parallel threads. Returns exactly
#' the same structure as `parallel_peptide()`.
#'
#' @param samples Matrix of bootstrap replicates (one row per replicate).
#' @param data Peptide x sample matrix.
#' @param model_formula Full model formula.
#' @param group Name of the group column.
#' @param meta.info Sample metadata, sample names as row names.
#' @param peptide_metadata Data frame mapping `peptide_id` to `protein_id`.
#' @param pSamples List of permuted metadata data frames.
#' @param pairwise_contrasts Single contrast string.
#' @param n_threads Number of Julia threads or `"auto"`.
#' @param BPPARAM Ignored (kept for a compatible signature): Julia threads
#' replace BiocParallel workers (forking R workers would also fork Julia).
#'
#' @return List with one element per replicate, each with `ds` and `pdps`.
parallel_peptide <- function(BPPARAM = NULL,
                                   samples,
                                   data,
                                   model_formula,
                                   group,
                                   meta.info,
                                   peptide_metadata,
                                   pSamples,
                                   pairwise_contrasts,
                                   n_threads) {

  ## Cheap checks first, before any fitting
  .check_peptides_per_protein(peptide_metadata)
  if (!is.null(BPPARAM)) {
    cli_alert_warning(
      "BPPARAM is ignored by the Julia backend; Julia threads are used instead."
    )
  }
  .setup_julia_peptide(n_threads)

  n_rep <- nrow(samples)
  results_list <- vector("list", n_rep)
  pb <- progress::progress_bar$new(format = "Fitting [:bar] :percent | ETA: :eta", 
            total = n_rep + 1, clear = FALSE)

  for (i in seq_len(n_rep)) {

    samples.R <- split(samples[i, ], meta.info[samples[i, ], group])
    pSamples_i <- pSamples[[i]]

    fit <- tryCatch(
      fit_bootstrap_peptide(
        x = lapply(samples.R, function(x) data[, x]),
        group = group, meta.info = meta.info,
        model_formula = model_formula,
        pairwise_contrasts = pairwise_contrasts,
        peptide_metadata = peptide_metadata
      ),
      error = function(e) {
        cli_abort("Bootstrap replicate {i} failed.", parent = e)
      }
    )
    df1 <- data.frame(d_result = fit$d, s_result = fit$s)

    pFit <- tryCatch(
      fit_permutation_peptide(
        x = data, group = group, meta.info = pSamples_i,
        model_formula = model_formula,
        pairwise_contrasts = pairwise_contrasts,
        peptide_metadata = peptide_metadata
      ),
      error = function(e) {
        cli_abort("Permutation replicate {i} failed.", parent = e)
      }
    )
    df2 <- data.frame(pd_result = pFit$d, ps_result = pFit$s)

    results_list[[i]] <- list(ds = df1, pdps = df2)
    pb$tick()
  }
  pb$tick()

  return(results_list)
}


#' Assemble Bootstrap Data and Matching Metadata
#'
#' Combine group-specific matrices and construct a metadata row for each
#' resampled sample while preserving their alignment.
#'
#' @param x List of feature-by-sample matrices, one per group in the bootstrap
#' replicate. Matrices must have matching peptide IDs as row names.
#' @param meta.info Data frame of sample-level metadata, with sample names as
#' row names.
#'
#' @return A list containing the combined feature-by-sample matrix
#' `combined_data` and corresponding sample-level metadata `covariates.p`.
#'
.generate_data_for_fit_peptide <- function(x, meta.info) {
  ### Check if x is a list of matrices with identical rownames
  stopifnot(
    all(vapply(
      x,
      function(z) identical(rownames(z), rownames(x[[1L]])),
      logical(1)
    ))
  )
  combined_data <- do.call(cbind, unname(x))
  
  ### create a unique colnames for the combined data to avoid duplicates
  colnames(combined_data) <-
    paste0(colnames(combined_data), ".", seq(1, ncol(combined_data)))
  
  ### create a covariates data frame for the combined data
  covariates.p <- data.frame()
  meta.info.temp <- meta.info
  meta.info.temp$sample.id <- row.names(meta.info.temp)
  for (i in colnames(combined_data)) {
    real_SampleNames <- str_split_fixed(i, fixed("."), 2)[, 1]
    df.temp <- meta.info.temp[real_SampleNames, , drop = FALSE]
    df.temp[,"sample.id"] <- i
    covariates.p <- rbind(covariates.p, df.temp)
  }
  stopifnot(identical(colnames(combined_data), covariates.p[,"sample.id"]))
  row.names(covariates.p) <- covariates.p[,"sample.id"]
  covariates.p[,"sample.id"] <- NULL
  
  stopifnot(identical(colnames(combined_data), row.names(covariates.p)))
  return(list(combined_data = combined_data, covariates.p = covariates.p))
}