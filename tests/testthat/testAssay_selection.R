library(testthat)
library(SummarizedExperiment)

test_that("Check_SummarizedExperiment correctly selects assay by name and index", {
    # Create mock data with two distinct assays
    mat1 <- matrix(runif(100), nrow = 10, ncol = 10)
    mat2 <- matrix(rnorm(100), nrow = 10, ncol = 10)
    rownames(mat1) <- rownames(mat2) <- paste0("Feature", 1:10)
    colnames(mat1) <- colnames(mat2) <- paste0("Sample", 1:10)
    
    col_data <- DataFrame(Group = rep(c("Control", "Treatment"), each = 5))
    rownames(col_data) <- colnames(mat1)
    
    se <- SummarizedExperiment(
        assays = list(first_assay = mat1, second_assay = mat2),
        colData = col_data
    )
    
    # Test default behavior (assay.name = NULL) -> should return the first assay
    res_default <- suppressMessages(Check_SummarizedExperiment(
        data.exp = se,
        meta.info = "Group",
        group.name = "Group",
        assay.name = NULL
    ))
    expect_equal(res_default$data, mat1)
    
    # Test selection by character name -> should return the second assay
    res_named <- suppressMessages(Check_SummarizedExperiment(
        data.exp = se,
        meta.info = "Group",
        group.name = "Group",
        assay.name = "second_assay"
    ))
    expect_equal(res_named$data, mat2)
    
    # Test selection by numeric index -> should return the second assay
    res_index <- suppressMessages(Check_SummarizedExperiment(
        data.exp = se,
        meta.info = "Group",
        group.name = "Group",
        assay.name = 2
    ))
    expect_equal(res_index$data, mat2)
})
