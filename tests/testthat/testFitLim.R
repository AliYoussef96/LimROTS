library(testthat)
library(limma)

test_that("Limma_fit returns expected output structure", {
    # Sample data with 10 samples per group
    data <- list(
        group1 = t(matrix(rnorm(400), nrow = 10)), # 10 samples, 4 features
        group2 = t(matrix(rnorm(400), nrow = 10)) # 10 samples, 4 features
    )

    # Sample meta info with corresponding group labels
    meta.info <- data.frame(
        sample_id = c(paste0("sample", 1:10), paste0("sample", 11:20)),
        group = as.factor(c(rep(1, 10), rep(2, 10)))
    )

    # Formula string for modeling
    formula.str <- "~ 0 + group"

    # Call the function
    result <- Limma_fit(data, "group", meta.info, formula.str, trend = FALSE, robust = FALSE)

    # Check output structure
    expect_type(result, "list")
    expect_equal(length(result), 3) # Expecting 3 elements in the output list
    expect_true("d" %in% names(result)) # Check for 'd' in result
    expect_true("s" %in% names(result)) # Check for 's' in result
    expect_true("corrected.logfc" %in% names(result)) # Check for 'corrected.logfc' in result
})
