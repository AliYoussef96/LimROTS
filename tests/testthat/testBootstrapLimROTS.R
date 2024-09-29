library(testthat)

test_that("bootstrapSamples.limRots returns expected output structure", {
    # Sample meta info with two groups and one factor
    meta.info <- data.frame(
        row.names = paste0("sample", 1:20),
        group = c(rep(1, 10), rep(2, 10)),
        factor1 = factor(rep(c("A", "B"), each = 10))
    )
    row.names(meta.info) <- meta.info$sample_id

    # Parameters
    B <- 10 # Number of bootstrap samples
    group.name <- "group"

    # Call the function
    result <- bootstrapSamples.limRots(data = NULL, B = B, meta.info = meta.info, group.name = group.name)

    # Check output structure
    expect_type(result, "character")
    expect_equal(dim(result), c(B, nrow(meta.info))) # Expecting B rows and nrow(meta.info) columns

})
