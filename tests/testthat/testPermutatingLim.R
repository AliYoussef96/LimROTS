library(testthat)
library(limma)
library(stringr)

test_that("Limma_permutating returns expected output structure", {
    # Sample data with 10 samples per group
    data <- data.frame(matrix(rnorm(500), nrow = 100, ncol = 10))
    meta.info <- data.frame(
        group = factor(rep(1:2, each = 5)),
        row.names = colnames(data)
    )
    formula.str <- "~ 0 + group"

    # Call the function
    result <- Limma_permutating(data, "group", meta.info, formula.str, trend = FALSE, robust = FALSE, permutating.group = TRUE)

    # Check output structure
    expect_type(result, "list")
    expect_equal(length(result), 2) # Expecting 3 elements in the output list
    expect_true("d" %in% names(result)) # Check for 'd' in result
    expect_true("s" %in% names(result)) # Check for 's' in result
})
