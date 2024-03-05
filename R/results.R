rna_velo_knn <- function(results, velocity, perplexity = 70, randseed = 0, raw = FALSE) {
  # set.seed(randseed)
  counts_s <- results$counts
  pop <- results$cell_meta$pop
  depth <- results$cell_meta$depth
  
  counts_s_lg <- t(log2(counts_s + 1))
  
  if (is.null(results$velocity)) {
    stop("The result object is not produced in velocity mode.")
  }
  
  process_velocity <- function(v) {
    assertthat::assert_that(
      nrow(counts_s) == nrow(v),
      ncol(counts_s) == ncol(v)
    )
    
    future_counts_s <- counts_s + v
    future_counts_s[future_counts_s < 0] <- 0
    future_counts_s_lg <- t(log2(future_counts_s + 1))
    future_counts_s_lg - counts_s_lg
  }
  
  
  normalize_velocity <- function(v) {
    v_normalizer <- apply(v, 2, \(vi) vi^2) %>% rowSums() %>% sqrt()
    t(t(v) / v_normalizer)
  }
  
  if (raw) {
    return(
      paired_simil(velocity, results$velocity, method = "cosine")
    )
  }
  
  dist_obj <- dist(counts_s_lg)
  dist_mat <- as.matrix(dist_obj)
  n_cells <- nrow(dist_mat)
  k <- ceiling(n_cells / 50)
  
  v_knn <- process_velocity(velocity) %>%
    apply(2, \(vi)
      distMat.KernelKnn(dist_mat, TEST_indices = NULL,
                        weights_function = 'gaussian',
                        y = vi, k = k, regression = TRUE)
    ) %>%
    normalize_velocity()
  
  v_true_knn <- process_velocity(results$velocity) %>%
    apply(2, \(vi)
      distMat.KernelKnn(dist_mat, TEST_indices = NULL,
                        weights_function = 'gaussian',
                        y = vi, k = k, regression = TRUE)
    ) %>%
    normalize_velocity()
  
  sim <- paired_simil(v_knn, v_true_knn, method = "cosine")
  
  mean(sim)
}