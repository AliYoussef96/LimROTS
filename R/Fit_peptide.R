# ----------------------------------------------------------------------------
# Peptide-level fits with Julia (MixedModels.jl): BiocParallel only
#
# Same statistical model as the original code:
#     y ~ <fixed part> + (1 | peptide_id) + (1 | sample_id), REML,
#     group contrast from the pairwise_contrasts string (asymptotic SE).
#
# Division of work
#   R     : sampling, data preparation, all checks, the design matrix X and the
#           contrast vector L (built with emmeans, see .peptide_design).
#   Julia : fits the models of ONE replicate (one model per protein, one
#           thread). There are NO thread settings: all parallelism comes from
#           BiocParallel (BPPARAM handled as in the original code).
#
# Where Julia runs
#   Julia is NEVER started in the main R session before the workers. Every
#   BiocParallel worker (SnowParam or MulticoreParam) starts its own,
#   independent Julia process the first time it needs it. In the main session
#   Julia is only started afterwards, by fit_original_peptide(), and stopped
#   when that fit is done.
#
# Julia is only needed when LimROTS_peptide() is called: nothing here runs
# at package build/check time. The JuliaConnectoR package is only listed under
# Suggests and is only used (via `::`) when a Julia function is actually called.
# ----------------------------------------------------------------------------

## Session-level state (the Julia process and the loaded Julia function)
.limrots_julia <- new.env(parent = emptyenv())


#' Folder of the private Julia environment of LimROTS
#'
#' @return Path (created lazily by Julia's package manager).
.julia_env_dir <- function() {
  normalizePath(
    file.path(tools::R_user_dir("LimROTS", which = "data"), "julia_env"),
    winslash = "/", mustWork = FALSE
  )
}


#' Are the Julia packages already installed? (does NOT start Julia)
#'
#' Reads `Project.toml` of the private LimROTS Julia environment, so the main R
#' session can check this without starting a Julia process.
#'
#' @return `TRUE` or `FALSE`.
.julia_packages_installed <- function() {
  project_file <- file.path(.julia_env_dir(), "Project.toml")
  if (!file.exists(project_file)) {
    return(FALSE)
  }
  lines <- readLines(project_file, warn = FALSE)
  all(vapply(c("MixedModels", "DataFrames"), function(p) {
    any(grepl(paste0("^", p, "[[:space:]]*="), lines))
  }, logical(1)))
}


#' Install the Julia packages (run ONCE, interactively)
#'
#' Installs MixedModels.jl and DataFrames.jl in a private environment under
#' `tools::R_user_dir("LimROTS", "data")`, never in the user's global Julia
#' environment, and only after consent (interactive prompt or
#' `options(LimROTS.julia_install = TRUE)`). It starts a temporary Julia
#' process and stops it again when done.
#'
#' @return `TRUE` invisibly.
.install_julia_peptide <- function() {

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

  env_dir <- .julia_env_dir()
  dir.create(env_dir, recursive = TRUE, showWarnings = FALSE)

  Sys.setenv(JULIA_NUM_THREADS = "1")
  missing_pkgs <- JuliaConnectoR::juliaEval(sprintf('begin
    import Pkg
    Pkg.activate("%s"; io = devnull)
    [p for p in ["MixedModels", "DataFrames"]
       if !(p in keys(Pkg.project().dependencies))]
  end', env_dir))
  missing_pkgs <- as.character(unlist(missing_pkgs))

  if (length(missing_pkgs) == 0) {
    cli_alert_success("The Julia packages are already installed.")
  } else {
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
    cli_alert_success("Julia packages installed.")
  }

  ## Leave no Julia process behind
  tryCatch(JuliaConnectoR::stopJulia(), error = function(e) NULL)
  .limrots_julia$ready <- FALSE
  invisible(TRUE)
}


#' Start Julia in THIS R process and load limrots_peptide_final.jl
#'
#' Runs once per R process (a no-op afterwards). It is called by every
#' BiocParallel worker, so each worker gets its own independent Julia process,
#' and by `fit_original_peptide()` in the main session. It never installs
#' anything: if the Julia packages are missing it aborts and points to
#' `.install_julia_peptide()`.
#'
#' @return `TRUE` invisibly.
.setup_julia_peptide <- function() {

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

  ## One Julia thread per process (the code is single-threaded; this also
  ## overrides a global JULIA_NUM_THREADS so idle threads are not created).
  Sys.setenv(JULIA_NUM_THREADS = "1")

  ## Private Julia environment of LimROTS
  env_dir <- .julia_env_dir()
  missing_pkgs <- JuliaConnectoR::juliaEval(sprintf('begin
    import Pkg
    Pkg.activate("%s"; io = devnull)
    [p for p in ["MixedModels", "DataFrames"]
       if !(p in keys(Pkg.project().dependencies))]
  end', env_dir))
  if (length(missing_pkgs) > 0) {
    cli_abort(c(
      "The Julia packages (MixedModels, DataFrames) are not installed.",
      "i" = "Run {.code LimROTS:::.install_julia_peptide()} once in an
             interactive session."
    ))
  }

  ## Load the Julia file that ships with the package
  julia_file <- system.file("julia", "limrots_peptide_final.jl",
                            package = "LimROTS")
  if (!nzchar(julia_file)) {
    cli_abort("Could not find {.file inst/julia/limrots_peptide_final.jl}.")
  }
  JuliaConnectoR::juliaCall("include", julia_file)
  .limrots_julia$fit <- JuliaConnectoR::juliaFun("limrots_fit_proteins")

  .limrots_julia$ready <- TRUE
  invisible(TRUE)
}


#' Stop the Julia process of this R process (if LimROTS started one)
#'
#' Frees the memory and the threads. The next call that needs Julia simply
#' starts a new Julia process (which pays the start-up and compile time again).
#' Set `options(LimROTS.julia_keep_alive = TRUE)` to keep Julia running after
#' `fit_original_peptide()`, e.g. when testing repeatedly.
#'
#' @return `TRUE` invisibly if Julia was stopped, `FALSE` otherwise.
.stop_julia_peptide <- function() {
  if (!isTRUE(.limrots_julia$ready)) {
    return(invisible(FALSE))
  }
  stopped <- tryCatch({
    JuliaConnectoR::stopJulia()
    TRUE
  }, error = function(e) {
    cli_alert_warning("Could not stop Julia: {conditionMessage(e)}")
    FALSE
  })
  if (stopped) {
    ## Forget everything tied to the old process
    .limrots_julia$ready <- FALSE
    .limrots_julia$faster_ready <- FALSE
    .limrots_julia$fit <- NULL
    .limrots_julia$fit_all <- NULL
  }
  invisible(stopped)
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


#' Build the fixed-effect design matrix and the contrast vector
#'
#' The design matrix comes from `model.matrix()`. The contrast vector `L` comes
#' from the emmeans package: a fixed-effects-only `lm()` with the same fixed
#' formula is used as a carrier, and `emmeans::emmeans()` followed by
#' `emmeans::contrast()` give the linear function of the coefficients (the
#' response of the carrier `lm()` is irrelevant: `L` depends only on the
#' design). The contrast weights come from `limma::makeContrasts`, as in the
#' original code. Julia then only applies `L` to the coefficients it estimates.
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
  reserved <- intersect(vars, c("y", "peptide_id", "sample_id", "protein_id",
                                ".y_dummy"))
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

  ## --- Contrast weights: identical call as in the original code ------------
  weights <- limma::makeContrasts(
    contrasts = stringr::str_remove_all(pairwise_contrasts,
                                        stringr::fixed(group)),
    levels = group_levels
  )[, 1]

  ## --- Contrast vector L from emmeans --------------------------------------
  ## Carrier model: same fixed formula, dummy response (L does not depend on it)
  cov[[".y_dummy"]] <- seq_len(nrow(cov))
  carrier_formula <- stats::reformulate(
    attr(stats::terms(fixed), "term.labels"), response = ".y_dummy"
  )
  carrier <- stats::lm(carrier_formula, data = cov)

  emm <- suppressMessages(
    emmeans::emmeans(carrier, specs = group, data = cov)
  )
  ## The weights are defined per group level, so the grid order must match
  if (!identical(as.character(emm@grid[[group]]), group_levels)) {
    cli_abort("emmeans group levels do not match the factor levels of {.val {group}}.")
  }
  con <- suppressMessages(
    emmeans::contrast(
      emm,
      method = stats::setNames(list(unname(weights)), pairwise_contrasts),
      adjust = "none"
    )
  )

  ## Linear function of the lm coefficients (same order as model.matrix columns)
  L <- as.numeric(con@linfct[1, ])
  if (length(L) != ncol(X_full) || anyNA(L)) {
    cli_abort("The group contrast is not estimable (or has the wrong length).")
  }
  if (!is.null(colnames(con@linfct)) &&
      !identical(colnames(con@linfct), colnames(X_full))) {
    cli_abort("emmeans coefficient names do not match the design matrix columns.")
  }

  list(X = unname(X_full[, -1, drop = FALSE]), L = L)
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
#' Runs in the main R session, AFTER the workers are done. Starts a Julia
#' process here and stops it again at the end (the last Julia use of a
#' `LimROTS_peptide()` run). Set `options(LimROTS.julia_keep_alive = TRUE)` to
#' keep Julia running, e.g. when testing repeatedly.
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
  if (!isTRUE(getOption("LimROTS.julia_keep_alive", FALSE))) {
    on.exit(.stop_julia_peptide(), add = TRUE)
  }
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


#' Julia version of parallel_peptide(): BPPARAM as in the original code
#'
#' The replicates are distributed by `bplapply()` exactly as before, with the
#' same default backend and the same `BPPARAM` handling. There are no thread
#' settings: Julia runs single-threaded in every process.
#'
#' Julia is NOT started in the main session. Every worker (`SnowParam` or
#' `MulticoreParam`) starts its own independent Julia process on its first
#' replicate. Therefore:
#' * the Julia packages must be installed beforehand (a worker cannot ask for
#'   consent): run `LimROTS:::.install_julia_peptide()` once, interactively;
#' * with `MulticoreParam` (forked workers), Julia must not be running in the
#'   main session, otherwise the forks would inherit its connection.
#'   (`fit_original_peptide()` stops Julia when it is done, so a second run in
#'   the same session is fine.)
#'
#' With `SerialParam` the "worker" is the main session, so Julia starts there
#' and is stopped by `fit_original_peptide()` at the end.
#'
#' @param BPPARAM A `BiocParallelParam` object. If `NULL`, a two-worker
#' backend appropriate to the operating system is used.
#' @param samples Matrix of bootstrap replicates (one row per replicate).
#' @param data Peptide x sample matrix.
#' @param model_formula Full model formula.
#' @param group Name of the group column.
#' @param meta.info Sample metadata, sample names as row names.
#' @param peptide_metadata Data frame mapping `peptide_id` to `protein_id`.
#' @param pSamples List of permuted metadata data frames.
#' @param pairwise_contrasts Single contrast string.
#'
#' @return List with one element per replicate, each with `ds`
#' (`d_result`, `s_result`) and `pdps` (`pd_result`, `ps_result`), row names =
#' proteins.
parallel_peptide <- function(BPPARAM = NULL,
                             samples,
                             data,
                             model_formula,
                             group,
                             meta.info,
                             peptide_metadata,
                             pSamples,
                             pairwise_contrasts) {

  .check_peptides_per_protein(peptide_metadata)

  if (is.null(BPPARAM)) {
    if (.Platform$OS.type == "windows") {
      BPPARAM <- SnowParam(workers = 2)
      cli_alert_info("Using SnowParam (Windows) with two workers.")
    } else {
      BPPARAM <- MulticoreParam(workers = 2)
      cli_alert_info("Using MulticoreParam (Unix-like OS) with two workers.")
    }
  } else {
    cli_alert_info("Using the provided parallel backend.")
  }
  if (inherits(BPPARAM, "SnowParam")) {
    BPPARAM$exportglobals <- FALSE
  }

  ## ---- Julia is NOT started here; only checks that do not start it ----------
  if (!requireNamespace("JuliaConnectoR", quietly = TRUE)) {
    cli_abort("The Julia backend needs the {.pkg JuliaConnectoR} R package.")
  }
  if (!.julia_packages_installed()) {
    cli_abort(c(
      "The Julia packages (MixedModels, DataFrames) are not installed yet.",
      "i" = "Run {.code LimROTS:::.install_julia_peptide()} once in an
             interactive session."
    ))
  }
  if (inherits(BPPARAM, "MulticoreParam") && isTRUE(.limrots_julia$ready)) {
    cli_abort(c(
      "Julia is already running in this R session, so forked workers would inherit it.",
      "i" = "Restart R, or use {.fn BiocParallel::SnowParam}."
    ))
  }

  export_vars <- list(
    samples = samples,
    data = data,
    model_formula = model_formula,
    group = group,
    meta.info = meta.info,
    peptide_metadata = peptide_metadata,
    pSamples = pSamples,
    pairwise_contrasts = pairwise_contrasts
  )
  for (name in names(export_vars)) {
    assign(name, export_vars[[name]])
  }


  results_list <- bplapply(seq_len(nrow(samples)), function(i) {

    ## Once per R process (a no-op afterwards): starts THIS worker's own Julia
    .setup_julia_peptide()

    samples.R <- split(samples[i, ], meta.info[samples[i, ], group])
    pSamples_i <- pSamples[[i]]
    d_result <- NULL
    s_result <- NULL
    pd_result <- NULL
    ps_result <- NULL
    fit <- tryCatch(
      fit_bootstrap_peptide(
        x = lapply(samples.R, function(x)
          data[, x]
        ),
        group = group,
        meta.info = meta.info,
        model_formula = model_formula,
        pairwise_contrasts = pairwise_contrasts,
        peptide_metadata = peptide_metadata
      ),
      error = function(e) {
        cli_abort("Bootstrap replicate {i} failed.", parent = e)
      }
    )
    d_result <- fit$d
    s_result <- fit$s
    df1 <- data.frame(d_result = d_result, s_result = s_result)
    pFit <- tryCatch(
      fit_permutation_peptide(
        x = data,
        group = group,
        meta.info = pSamples_i,
        model_formula = model_formula,
        pairwise_contrasts = pairwise_contrasts,
        peptide_metadata = peptide_metadata
      ),
      error = function(e) {
        cli_abort("Permutation replicate {i} failed.", parent = e)
      }
    )
    pd_result <- pFit$d
    ps_result <- pFit$s
    df2 <- data.frame(pd_result = pd_result, ps_result = ps_result)
    list(ds = df1, pdps = df2)
  }, BPPARAM = BPPARAM,
  BPOPTIONS = bpoptions(packages = c("utils",
                                     "stringr", "stats",
                                     "data.table", "lme4", "emmeans",
                                     "limma")))
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