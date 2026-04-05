library(testthat)
library(SummarizedExperiment)

# Helper: build a SummarizedExperiment with survival colData
make_se_survival <- function(nfeatures = 50, nsamples = 30, seed = 42) {
    set.seed(seed)
    mat <- matrix(rnorm(nfeatures * nsamples), nrow = nfeatures,
                  ncol = nsamples)
    rownames(mat) <- paste0("feature", seq_len(nfeatures))
    colnames(mat) <- paste0("sample", seq_len(nsamples))
    col_data <- data.frame(
        time  = abs(rnorm(nsamples, mean = 5, sd = 2)) + 0.1,
        event = sample(0:1, nsamples, replace = TRUE),
        group = factor(rep(1:2, each = nsamples / 2)),
        row.names = colnames(mat)
    )
    SummarizedExperiment(assays = list(expr = mat),
                         colData = col_data)
}

# ── LimROTS_survival returns SummarizedExperiment ────────────────────────────
test_that("LimROTS_survival returns a SummarizedExperiment when input is SE", {
    se <- make_se_survival()
    result <- LimROTS_survival(
        x = se,
        meta.info = c("time", "event", "group"),
        formula.str = "~ Surv(time, event) + group",
        niter = 10,
        verbose = FALSE,
        BPPARAM = BiocParallel::SerialParam()
    )
    expect_s4_class(result, "SummarizedExperiment")
})

# ── rowData contains expected result columns ──────────────────────────────────
test_that("LimROTS_survival rowData contains expected result columns", {
    se <- make_se_survival()
    result <- LimROTS_survival(
        x = se,
        meta.info = c("time", "event", "group"),
        formula.str = "~ Surv(time, event) + group",
        niter = 10,
        verbose = FALSE,
        BPPARAM = BiocParallel::SerialParam()
    )
    expected_cols <- c("statistics", "pvalue", "FDR", "exp_coef", "BH.pvalue")
    expect_true(all(expected_cols %in% colnames(rowData(result))))
})

# ── metadata contains optimisation parameters ─────────────────────────────────
test_that("LimROTS_survival metadata contains optimisation parameters", {
    se <- make_se_survival()
    result <- LimROTS_survival(
        x = se,
        meta.info = c("time", "event", "group"),
        formula.str = "~ Surv(time, event) + group",
        niter = 10,
        verbose = FALSE,
        BPPARAM = BiocParallel::SerialParam()
    )
    meta <- metadata(result)
    expect_true(all(c("a1", "a2", "k", "Z", "R") %in% names(meta)))
    expect_true(is.numeric(meta$a1))
    expect_true(is.numeric(meta$a2))
})

# ── output retains all input features ────────────────────────────────────────
test_that("LimROTS_survival output retains all input features", {
    se <- make_se_survival(nfeatures = 50)
    result <- LimROTS_survival(
        x = se,
        meta.info = c("time", "event", "group"),
        formula.str = "~ Surv(time, event) + group",
        niter = 10,
        verbose = FALSE,
        BPPARAM = BiocParallel::SerialParam()
    )
    expect_equal(nrow(result), nrow(se))
    expect_equal(rownames(result), rownames(se))
})

# ── p-values are in [0, 1] ────────────────────────────────────────────────────
test_that("LimROTS_survival p-values are in [0, 1]", {
    se <- make_se_survival()
    result <- LimROTS_survival(
        x = se,
        meta.info = c("time", "event", "group"),
        formula.str = "~ Surv(time, event) + group",
        niter = 10,
        verbose = FALSE,
        BPPARAM = BiocParallel::SerialParam()
    )
    pvals <- rowData(result)$pvalue
    expect_true(all(pvals >= 0 & pvals <= 1, na.rm = TRUE))
})

# ── FDR values are in [0, 1] ─────────────────────────────────────────────────
test_that("LimROTS_survival FDR values are in [0, 1]", {
    se <- make_se_survival()
    result <- LimROTS_survival(
        x = se,
        meta.info = c("time", "event", "group"),
        formula.str = "~ Surv(time, event) + group",
        niter = 10,
        verbose = FALSE,
        BPPARAM = BiocParallel::SerialParam()
    )
    fdr <- rowData(result)$FDR
    expect_true(all(fdr >= 0 & fdr <= 1, na.rm = TRUE))
})

# ── exp_coef contains numeric hazard ratio estimates ─────────────────────────
test_that("LimROTS_survival exp_coef contains numeric values", {
    se <- make_se_survival()
    result <- LimROTS_survival(
        x = se,
        meta.info = c("time", "event", "group"),
        formula.str = "~ Surv(time, event) + group",
        niter = 10,
        verbose = FALSE,
        BPPARAM = BiocParallel::SerialParam()
    )
    exp_coef <- rowData(result)$exp_coef
    expect_true(is.numeric(exp_coef))
    expect_true(all(exp_coef > 0, na.rm = TRUE))
})

# ── user-supplied a1/a2 bypasses optimisation ────────────────────────────────
test_that("LimROTS_survival accepts user-supplied a1 and a2", {
    se <- make_se_survival()
    result <- LimROTS_survival(
        x = se,
        meta.info = c("time", "event", "group"),
        formula.str = "~ Surv(time, event) + group",
        niter = 10,
        a1 = 0.3,
        a2 = 0.2,
        verbose = FALSE,
        BPPARAM = BiocParallel::SerialParam()
    )
    expect_s4_class(result, "SummarizedExperiment")
    expect_equal(metadata(result)$a1, 0.3)
    expect_equal(metadata(result)$a2, 0.2)
})

# ── missing time/event columns raises an error ───────────────────────────────
test_that("LimROTS_survival errors when time or event columns are missing", {
    set.seed(1)
    nfeatures <- 20; nsamples <- 20
    mat <- matrix(rnorm(nfeatures * nsamples), nrow = nfeatures,
                  ncol = nsamples)
    rownames(mat) <- paste0("feature", seq_len(nfeatures))
    colnames(mat) <- paste0("sample", seq_len(nsamples))
    col_data <- data.frame(
        group = factor(rep(1:2, each = nsamples / 2)),
        row.names = colnames(mat)
    )
    se <- SummarizedExperiment(assays = list(expr = mat),
                               colData = col_data)
    expect_error(
        LimROTS_survival(
            x = se,
            meta.info = c("group"),
            formula.str = "~ Surv(time, event) + group",
            niter = 10,
            verbose = FALSE,
            BPPARAM = BiocParallel::SerialParam()
        )
    )
})

# ── competing risks model runs without error ─────────────────────────────────
test_that("LimROTS_survival competing_risks = TRUE runs without error", {
    se <- make_se_survival(nsamples = 30)
    result <- LimROTS_survival(
        x = se,
        meta.info = c("time", "event", "group"),
        formula.str = "~ Surv(time, event) + group",
        niter = 10,
        verbose = FALSE,
        competing_risks = TRUE,
        BPPARAM = BiocParallel::SerialParam()
    )
    expect_s4_class(result, "SummarizedExperiment")
    expect_true("statistics" %in% colnames(rowData(result)))
})
