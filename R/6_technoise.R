add_expr_noise <- function(results, ...) {
  cat("Adding experimental noise...\n")
  data(gene_len_pool)
  gene_len <- sample(gene_len_pool, results$num_genes, replace = FALSE)
  args <- list(...)
  if (length(args) > 0) {
    rna_args <- args[names(args)[!startsWith(names(args), "atac.")]]
    atac_args <- args[names(args)[startsWith(names(args), "atac.")]]
  } else {
    rna_args <- list(); atac_args <- list()  
  }
  rna_args <- default.args(rna_args, randseed = 0,
                           protocol = "nonUMI", alpha_mean = 0.1, alpha_sd = 0.02,
                           gene_len = gene_len, depth_mean = 1e5, depth_sd = 3e3)
  atac_args <- default.args(atac_args, atac.obs.prob = 0.3, atac.sd.frac = 0.5)
  rna_args$true_counts <- floor(results$counts)
  rna_args$meta_cell <- results$cell_meta
  results$counts_obs <- do.call(True2ObservedCounts, rna_args)
  
  results$atacseq_obs <- True2ObservedATAC(results$atacseq_data, randseed = args$randseed,
                                           observation_prob = atac_args$atac.obs.prob,
                                           sd_frac = atac_args$atac.sd.frac) 
}


divide_batches <- function(results, nbatch = 2, effect = 3) {
  cat("Adding batch effects...\n")
  obs <- results$counts_obs
  if (is.list(obs)) {
    obs <- obs$counts
  }
  ngene <- nrow(obs)
  merged <- rbind(obs, results$atacseq_obs)
  if ("batch" %in% names(results$cell_meta)) {
    results$cell_meta <- results$cell_meta[, !(names(results$cell_meta) %in% "batch")]
  }
  b <- DivideBatches(counts = merged, meta_cell = results$cell_meta, nbatch = nbatch, batch_effect_size = effect)
  results$counts_with_batches <- b$counts[1:ngene,]
  results$atac_with_batches <- b$counts[-(1:ngene),]
  results$cell_meta <- b$cell_meta
}
