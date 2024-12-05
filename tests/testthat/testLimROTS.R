library(testthat)

test_that("LiROTS the main function returns expected output structure", {
  
  data <- data.frame(matrix(rnorm(500), nrow = 100, ncol = 10))
  meta.info <- data.frame(
      group = factor(rep(1:2, each = 5)),
      row.names = colnames(data)
  )
  formula.str <- "~ 0 + group"
  result <- LimROTS(data,
      meta.info = meta.info, group.name = "group",
      formula.str = formula.str, niter = 10, seed.cl = 1234 )
  
  expect_type(result, "list")
  expect_length(result, 16)
  
  expected_names <- c("data", "niter", "d", "logfc", "pvalue", "FDR",
                      "a1", "a2", "k", "R", "Z", "ztable", "groups",
                      "corrected.logfc" , "q_values", "BH.pvalue")
  
  expect_true(all(expected_names %in% names(result)), 
              info = paste("Expected names not found in the list."))
  
})
