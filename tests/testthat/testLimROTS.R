library(testthat)
library(SummarizedExperiment)

# Helper: build a small SummarizedExperiment for two-group comparison
make_se_two_group <- function(nfeatures = 50, nsamples = 10, seed = 42) {
    set.seed(seed)
    mat <- matrix(rnorm(nfeatures * nsamples), nrow = nfeatures,
                  ncol = nsamples)
    rownames(mat) <- paste0("feature", seq_len(nfeatures))
    colnames(mat) <- paste0("sample", seq_len(nsamples))
    col_data <- data.frame(
        group = factor(rep(1:2, each = nsamples / 2)),
        row.names = colnames(mat)
    )
    SummarizedExperiment(assays = list(expr = mat),
                         colData = col_data)
}

# Helper: multi-group SE (3 groups)
make_se_multi_group <- function(nfeatures = 50, nsamples = 15, seed = 7) {
    set.seed(seed)
    mat <- matrix(rnorm(nfeatures * nsamples), nrow = nfeatures,
                  ncol = nsamples)
    rownames(mat) <- paste0("feature", seq_len(nfeatures))
    colnames(mat) <- paste0("sample", seq_len(nsamples))
    col_data <- data.frame(
        group = factor(rep(1:3, each = nsamples / 3)),
        row.names = colnames(mat)
    )
    SummarizedExperiment(assays = list(expr = mat),
                         colData = col_data)
}

# ── LimROTS returns SummarizedExperiment ──────────────────────────────────────
test_that("LimROTS returns a SummarizedExperiment when input is SE", {
    se <- make_se_two_group()
    result <- LimROTS(
        x = se,
        assay.type = "expr",
        meta.info = c("group"),
        group = "group",
        formula.str = "~ 0 + group",
        niter = 10,
        verbose = FALSE,
        BPPARAM = BiocParallel::SerialParam()
    )
    expect_s4_class(result, "SummarizedExperiment")
})

# ── rowData contains expected result columns ───────────────────────────────────
test_that("LimROTS rowData contains expected result columns", {
    se <- make_se_two_group()
    result <- LimROTS(
        x = se,
        assay.type = "expr",
        meta.info = c("group"),
        group = "group",
        formula.str = "~ 0 + group",
        niter = 10,
        verbose = FALSE,
        BPPARAM = BiocParallel::SerialParam()
    )
    expected_cols <- c("statistics", "pvalue", "FDR",
                       "corrected.logfc", "BH.pvalue")
    expect_true(all(expected_cols %in% colnames(rowData(result))))
})

# ── metadata contains optimisation parameters ──────────────────────────────────
test_that("LimROTS metadata contains optimisation parameters", {
    se <- make_se_two_group()
    result <- LimROTS(
        x = se,
        assay.type = "expr",
        meta.info = c("group"),
        group = "group",
        formula.str = "~ 0 + group",
        niter = 10,
        verbose = FALSE,
        BPPARAM = BiocParallel::SerialParam()
    )
    meta <- metadata(result)
    expect_true(all(c("a1", "a2", "k", "Z", "R") %in% names(meta)))
    expect_true(is.numeric(meta$a1))
    expect_true(is.numeric(meta$a2))
})

# ── result has same features as input ─────────────────────────────────────────
test_that("LimROTS output retains all input features", {
    se <- make_se_two_group(nfeatures = 50)
    result <- LimROTS(
        x = se,
        assay.type = "expr",
        meta.info = c("group"),
        group = "group",
        formula.str = "~ 0 + group",
        niter = 10,
        verbose = FALSE,
        BPPARAM = BiocParallel::SerialParam()
    )
    expect_equal(nrow(result), nrow(se))
    expect_equal(rownames(result), rownames(se))
})

# ── pvalues are in [0, 1] ──────────────────────────────────────────────────────
test_that("LimROTS p-values are in [0, 1]", {
    se <- make_se_two_group()
    result <- LimROTS(
        x = se,
        assay.type = "expr",
        meta.info = c("group"),
        group = "group",
        formula.str = "~ 0 + group",
        niter = 10,
        verbose = FALSE,
        BPPARAM = BiocParallel::SerialParam()
    )
    pvals <- rowData(result)$pvalue
    expect_true(all(pvals >= 0 & pvals <= 1, na.rm = TRUE))
})

# ── FDR values are in [0, 1] ──────────────────────────────────────────────────
test_that("LimROTS FDR values are in [0, 1]", {
    se <- make_se_two_group()
    result <- LimROTS(
        x = se,
        assay.type = "expr",
        meta.info = c("group"),
        group = "group",
        formula.str = "~ 0 + group",
        niter = 10,
        verbose = FALSE,
        BPPARAM = BiocParallel::SerialParam()
    )
    fdr <- rowData(result)$FDR
    expect_true(all(fdr >= 0 & fdr <= 1, na.rm = TRUE))
})

# ── multi-group analysis runs without error ────────────────────────────────────
test_that("LimROTS handles multi-group SE input", {
    se <- make_se_multi_group()
    result <- LimROTS(
        x = se,
        assay.type = "expr",
        meta.info = c("group"),
        group = "group",
        formula.str = "~ 0 + group",
        niter = 10,
        verbose = FALSE,
        BPPARAM = BiocParallel::SerialParam()
    )
    expect_s4_class(result, "SummarizedExperiment")
    expect_equal(nrow(result), nrow(se))
})

# ── covariate model runs without error ────────────────────────────────────────
test_that("LimROTS handles SE with covariate in formula", {
    set.seed(1)
    nfeatures <- 50; nsamples <- 12
    mat <- matrix(rnorm(nfeatures * nsamples), nrow = nfeatures,
                  ncol = nsamples)
    rownames(mat) <- paste0("feature", seq_len(nfeatures))
    colnames(mat) <- paste0("sample", seq_len(nsamples))
    col_data <- data.frame(
        group = factor(rep(1:2, each = nsamples / 2)),
        batch = factor(rep(c("A", "B"), times = nsamples / 2)),
        row.names = colnames(mat)
    )
    se <- SummarizedExperiment(assays = list(expr = mat),
                               colData = col_data)
    result <- LimROTS(
        x = se,
        assay.type = "expr",
        meta.info = c("group", "batch"),
        group = "group",
        formula.str = "~ 0 + group + batch",
        niter = 10,
        verbose = FALSE,
        BPPARAM = BiocParallel::SerialParam()
    )
    expect_s4_class(result, "SummarizedExperiment")
})
