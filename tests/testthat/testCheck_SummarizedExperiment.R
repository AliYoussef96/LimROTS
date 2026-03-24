library(testthat)
library(SummarizedExperiment)

test_that("Check_SummarizedExperiment correctly selects assay", {
  # Create mock data with two distinct assays
  mat1 <- matrix(runif(100), nrow = 10, ncol = 10)
  mat2 <- matrix(rnorm(100), nrow = 10, ncol = 10)
  rownames(mat1) <- rownames(mat2) <- paste0("Feature", 1:10)
  colnames(mat1) <- colnames(mat2) <- paste0("Sample", 1:10)

  col_data <- DataFrame(Group = rep(c("Control", "Treatment"), each = 5))
  rownames(col_data) <- colnames(mat1)

  se <- SummarizedExperiment(
    assay = list(first_assay = mat1, second_assay = mat2),
    colData = col_data
  )

  # Test default behavior (assay.type = NULL) -> should return the first assay
  res_default <- suppressMessages(Check_SummarizedExperiment(
    data.exp = se,
    assay.type = NULL,
    meta.info = "Group",
    group.name = "Group"
  ))
  expect_equal(res_default$data, mat1)

  # Test selection by character string -> should return the second assay
  res_named <- suppressMessages(Check_SummarizedExperiment(
    data.exp = se,
    assay.type = "second_assay",
    meta.info = "Group",
    group.name = "Group"
  ))
  expect_equal(res_named$data, mat2)

  # Test selection by numeric index -> should return the second assay
  res_index <- suppressMessages(Check_SummarizedExperiment(
    data.exp = se,
    assay.type = 2,
    meta.info = "Group",
    group.name = "Group"
  ))
  expect_equal(res_index$data, mat2)
})