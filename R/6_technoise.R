add_expr_noise <- function(results) {
  cat("Adding experimental noise...\n")
  data(gene_len_pool)
  gene_len <- sample(gene_len_pool, results$num_genes, replace = FALSE)
  results$counts_obs <- True2ObservedCounts(true_counts = results$counts, meta_cell = results$cell_meta,
                                            randseed = 0,
                                            protocol = "nonUMI", alpha_mean = 0.1, alpha_sd = 0.02,
                                            gene_len = gene_len, depth_mean = 1e5, depth_sd = 3e3)
  results$atacseq_obs <- True2ObservedATAC(results$atacseq_data, randseed = 0, observation_prob = 0.3, sd_frac = 0.5) 
}

divide_batches <- function(results, nbatch = 2) {
  cat("Adding batch effects...\n")
  ngene <- nrow(results$counts_obs)
  merged <- rbind(results$counts_obs, results$atacseq_obs)
  b <- DivideBatches(counts = merged, meta_cell = results$cell_meta, nbatch = nbatch, batch_effect_size = 3)
  results$counts_with_batches <- b$counts[1:ngene,]
  results$atac_with_batches <- b$counts[-(1:ngene),]
  results$cell_meta <- b$cell_meta
}
