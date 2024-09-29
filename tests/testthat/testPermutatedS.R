library(testthat)

test_that("permutatedS returns expected output structure", {
    # Sample meta info with some samples
    meta.info <- data.frame(
        sample_id = paste0("sample", 1:20)
    )
    row.names(meta.info) <- meta.info$sample_id

    # Parameters
    B <- 5 # Number of permuted samples

    # Call the function
    result <- permutatedS(meta.info, B)

    # Check output structure
    expect_type(result, "character")
    expect_equal(dim(result), c(B, nrow(meta.info))) # Expecting B rows and nrow(meta.info) columns

    # Check that each row is a permutation of the original row names
    for (i in seq_len(B)) {
        expect_true(all(result[i, ] %in% row.names(meta.info))) # All sampled names should be in meta.info
        expect_equal(length(unique(result[i, ])), nrow(meta.info)) # Each row should have unique samples
    }
})
