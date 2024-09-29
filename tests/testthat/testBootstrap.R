library(testthat)

test_that("bootstrapS returns expected output structure", {
    # Sample meta info with two groups
    meta.info <- data.frame(
        sample_id = paste0("sample", 1:20),
        group = c(rep("group1", 10), rep("group2", 10))
    )
    row.names(meta.info) <- meta.info$sample_id

    # Parameters
    B <- 5 # Number of bootstrap samples
    group.name <- "group"
    paired <- FALSE # Change to TRUE for paired testing

    # Call the function
    result <- bootstrapS(B, meta.info, group.name, paired)

    # Check output structure
    expect_type(result, "character")
    expect_equal(dim(result), c(B, nrow(meta.info))) # Expecting B rows and nrow(meta.info) columns

    # Check that all samples are from the correct groups
    groups <- meta.info[, group.name]
    for (i in seq_len(B)) {
        for (g in unique(groups)) {
            expect_true(all(result[i, groups == g] %in% row.names(meta.info)[groups == g]))
        }
    }
})
