#' Optimize Parameters Based on Overlap Calculations
#'
#' This function optimizes parameters by calculating overlaps between observed
#' and permuted data for multiple values of a smoothing constant and a
#' single-label replicate (SLR) comparison.
#'
#' @param niter Integer. Number of bootstrap samples or resampling iterations.
#' @param smoothing_constants Numeric vector. Smoothing constants to be evaluated.
#' @param top_n_values Integer vector. Number of top values to consider for overlap
#' calculation.
#' @param observed_data Numeric matrix. Observed data values.
#' @param observed_std_errors Numeric matrix. Standard errors or related values for observed data.
#' @param permuted_data Numeric matrix. Permuted data values.
#' @param permuted_std_errors Numeric matrix. Standard errors or related values for
#' permuted data.
#' @param verbose Logical. If `TRUE`, progress messages will be displayed.
#'
#' @details
#' The function calculates overlaps for a range of smoothing constants and
#' identifies the optimal set of parameters by maximizing a z-score-based
#' metric, which compares the overlap of observed data to permuted data.
#' It computes overlap matrices for both observed and permuted data and 
#' returns the optimal parameters based on the highest z-score.
#'
#' @return A list containing the optimal parameters:
#' \itemize{
#'   \item \code{optimal_smoothing_constant}: Optimal smoothing constant or 1 
#'   for SLR (a1).
#'   \item \code{use_smoothing_flag}: SLR flag (1 if smoothing constant is 
#'   optimal, 0 if SLR is optimal) (a2).
#'   \item \code{optimal_top_n}: Optimal number of top values to consider 
#'   for overlap (k).
#'   \item \code{optimal_reproducibility}: Optimal overlap value (R).
#'   \item \code{optimal_z_score}: Optimal z-score (Z).
#'   \item \code{z_score_table}: Matrix of z-scores for all 
#'   evaluated parameters (ztable).
#' }
#'
#'
#'
Optimizing <- function(niter, smoothing_constants, top_n_values, 
                       observed_data, observed_std_errors, 
                       permuted_data, permuted_std_errors, verbose) {
  if (verbose) {
    message("Optimizing smoothing constant and method selection")
  }
  num_smoothing_strategies <- length(smoothing_constants)
  total_strategies <- num_smoothing_strategies + 1
  num_top_n_options <- length(top_n_values)
  num_features <- nrow(observed_data)
  
  observed_overlap_means <- matrix(
    nrow = total_strategies, 
    ncol = num_top_n_options,
    dimnames = list(c(smoothing_constants, "slr"), top_n_values)
  )
  
  permuted_overlap_means <- matrix(
    nrow = total_strategies, 
    ncol = num_top_n_options,
    dimnames = list(c(smoothing_constants, "slr"), top_n_values)
  )
  
  overlap_std_devs <- matrix(
    nrow = total_strategies, 
    ncol = num_top_n_options,
    dimnames = list(c(smoothing_constants, "slr"), top_n_values)
  )
  
  for (strategy_idx in seq_len(num_smoothing_strategies)) {
    current_smoothing <- smoothing_constants[strategy_idx]
    
    observed_bootstrap <- matrix(0, nrow = niter, ncol = num_top_n_options)
    permuted_bootstrap <- matrix(0, nrow = niter, ncol = num_top_n_options)
    
    overlap_results <- calOverlaps(
      observed_data, observed_std_errors, 
      permuted_data, permuted_std_errors, 
      num_features, 
      as.integer(top_n_values), 
      num_top_n_options,
      current_smoothing, 
      as.integer(niter), 
      observed_bootstrap, 
      permuted_bootstrap
    )
    
    observed_overlap_means[strategy_idx, ] <- 
      colMeans(overlap_results[["overlaps"]])
    
    permuted_overlap_means[strategy_idx, ] <- 
      colMeans(overlap_results[["overlaps_P"]])
    
    mean_values <- observed_overlap_means[strategy_idx, ]
    deviations <- t(overlap_results[["overlaps"]]) - mean_values
    sum_squared_deviations <- rowSums(deviations^2)
    overlap_std_devs[strategy_idx, ] <- 
      sqrt(sum_squared_deviations / (niter - 1))
  }
  
  slr_strategy_idx <- total_strategies
  
  observed_bootstrap_slr <- matrix(0, nrow = niter, ncol = num_top_n_options)
  permuted_bootstrap_slr <- matrix(0, nrow = niter, ncol = num_top_n_options)
  
  overlap_results_slr <- calOverlaps_slr(
    observed_data, permuted_data, 
    num_features, 
    as.integer(top_n_values), 
    num_top_n_options,
    as.integer(niter), 
    observed_bootstrap_slr, 
    permuted_bootstrap_slr
  )
  
  # Compute statistics for SLR strategy
  observed_overlap_means[slr_strategy_idx, ] <- 
    colMeans(overlap_results_slr[["overlaps"]])
  
  permuted_overlap_means[slr_strategy_idx, ] <- 
    colMeans(overlap_results_slr[["overlaps_P"]])
  
  mean_values_slr <- observed_overlap_means[slr_strategy_idx, ]
  deviations_slr <- t(overlap_results_slr[["overlaps"]]) - mean_values_slr
  sum_squared_deviations_slr <- rowSums(deviations_slr^2)
  overlap_std_devs[slr_strategy_idx, ] <- 
    sqrt(sum_squared_deviations_slr / (niter - 1))
  
  z_score_matrix <- (observed_overlap_means - permuted_overlap_means) / 
    overlap_std_devs
  
  finite_z_scores <- is.finite(z_score_matrix)
  max_z_value <- max(z_score_matrix[finite_z_scores])
  optimal_position <- which(
    z_score_matrix == max_z_value, 
    arr.ind = TRUE
  )
  
  if (nrow(optimal_position) > 1) {
    optimal_position <- optimal_position[1, ]
  }
  
  optimal_strategy_row <- optimal_position[1]
  optimal_top_n_col <- optimal_position[2]
  
  if (optimal_strategy_row <= num_smoothing_strategies) {
    optimal_smoothing_constant <- as.numeric(
      rownames(observed_overlap_means)[optimal_strategy_row]
    )
    use_smoothing_flag <- 1
  } else {
    optimal_smoothing_constant <- 1
    use_smoothing_flag <- 0
  }
  
  optimal_top_n <- as.numeric(
    colnames(observed_overlap_means)[optimal_top_n_col]
  )
  
  optimal_reproducibility <- observed_overlap_means[
    optimal_strategy_row, 
    optimal_top_n_col
  ]
  
  optimal_z_score <- z_score_matrix[
    optimal_strategy_row, 
    optimal_top_n_col
  ]
  
  gc()
  
  return(list(
    optimal_smoothing_constant = optimal_smoothing_constant,
    use_smoothing_flag = use_smoothing_flag,
    optimal_top_n = optimal_top_n,
    optimal_reproducibility = optimal_reproducibility,
    optimal_z_score = optimal_z_score,
    z_score_table = z_score_matrix
  ))
}