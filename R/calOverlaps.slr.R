#' @export

calOverlaps.slr <- function(D, pD, nrow, N, N_len, B, overlaps, overlaps_P) {

  # sort based on one vector and operate on another
  sort2_1R <- function(a, b) {
    order_a <- order(a, b, decreasing = TRUE)
    a <- a[order_a]
    b <- b[order_a]
    list(a = a, b = b)
  }


  idx_b <- seq_len(B)
  idx_offset <- B

  D <- abs(D)

  for (b in idx_b) {
    # Calculate res1, res2, pres1, pres2
    res1 <- abs(D[((b - 1) * nrow + 1):(b * nrow)])
    res2 <- abs(D[((b + idx_offset - 1) * nrow + 1):((b + idx_offset) * nrow)])
    pres1 <- abs(pD[((b - 1) * nrow + 1):(b * nrow)])
    pres2 <- abs(pD[((b + idx_offset - 1) * nrow + 1):((b + idx_offset) * nrow)])

    # Overlap calculation for res1/res2
    sorted_res <- sort2_1R(res1, res2)
    res1 <- sorted_res$a
    res2 <- sorted_res$b
    r3_res <- sort(res2, decreasing = TRUE)

    for (i in seq_len(N_len)) {
      N_i <- N[i]
      sum_overlap <- sum(res2[1:N_i] >= r3_res[N_i])
      overlaps[b, i] <- sum_overlap / N_i
    }

    # Overlap calculation for pres1/pres2
    sorted_pres <- sort2_1R(pres1, pres2)
    pres1 <- sorted_pres$a
    pres2 <- sorted_pres$b
    r3_pres <- sort(pres2, decreasing = TRUE)

    for (i in seq_len(N_len)) {
      N_i <- N[i]
      sum_overlap <- sum(pres2[1:N_i] >= r3_pres[N_i])
      overlaps_P[b, i] <- sum_overlap / N_i
    }
  }

  return(list(overlaps = overlaps, overlaps_P = overlaps_P))
}


