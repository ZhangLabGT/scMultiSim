.scMultiSim_debug_result <- T

.res_opt <- function(...) {
  opts <- list2(...)
  if (!("debug" %in% opts)) {
    opts$debug = .scMultiSim_debug_result
  }
  opts
}

.res_obj <- function(res, ...) {
  arg_list <- list2(...)
  fig <- list()
  for (i in arg_list) {
    if ((".type" %in% names(i)) && i$.type == "plot") {
      fig[[i$name]] <- i$plot
    }
  }
  obj <- list(res = res, fig = fig)
  class(obj) <- c("scMultiSim_output", class(obj))
  obj
}


write_files <- function(res, dataset_name, dir = "sim", prefix = NULL, suffix = NULL) {
  if (is(res, "scMultiSim_output")) {
    res <- res$res
  }
  
  data <- list(
    counts_s = res$counts %>% t(),
    atac = res$atacseq_data %>% t(),
    meta = res$cell_meta
  )
  
  if ("unspliced_counts" %in% names(res)) {
    data$counts_u <- res$unspliced_counts %>% t()
    data$velo <- res$velocity %>% t()
  }
  
  if ("counts_with_batches" %in% names(res)) {
    data$counts_obs <- res$counts_with_batches %>% t()
  }
  
  if ("atac_with_batches" %in% names(res)) {
    data$atac_obs <- res$atac_with_batches %>% t()
  }
  
  ds_path <- file.path(dir, dataset_name)
  if (!file.exists(ds_path)) {
    dir.create(ds_path, recursive = T)
  }
  
  for (n in names(data)) {
    d <- data[[n]]
    fn <- n
    if (is.character(prefix)) {
      fn <- paste0(prefix, "_", fn)
    }
    if (is.character(suffix)) {
      fn <- paste0(fn, "_", suffix)
    }
    write.csv(d, file = file.path(ds_path, paste0(fn, ".csv")), quote = F)
  }
  
  saveRDS(res, file.path(ds_path, "res.rds"))
  
  res
}


.plt_tsne <- function(data, name, legend = "pop", tsne.seed = 0) {
  results <- get("results", envir = caller_env())
  lst(name,
      .type = "plot",
      plot = plot_tsne(
        log2(data + 1),
        results$cell_meta[[legend]],
        legend = legend, plot.name = name, rand.seed = tsne.seed
      ))
}


Bench <- list()


Bench$rna_atac_nongrn <- function(velo = T, diff.frac = 0.8, intr.noise = 1, tsne.seed = 0) {
  options_ <- .res_opt(
    rand.seed = 1,
    GRN = NA,
    num.cells = 1000,
    num.cifs = 500,
    tree = Phyla5(),
    diff.cif.fraction = diff.frac,
    intrinsic.noise = intr.noise,
    do.velocity = velo
  )
  
  results <- sim_true_counts(options_)
  add_expr_noise(results)
  divide_batches(results)
  
  .res_obj(results,
    .plt_tsne(results$counts, "True mRNA TSNE"),
    .plt_tsne(results$atacseq_data, "True ATAC TSNE"),
    .plt_tsne(results$counts_with_batches, "Observed mRNA (pop)", "pop"),
    .plt_tsne(results$counts_with_batches, "Observed mRNA (batches)", "batch"),
    .plt_tsne(results$atac_with_batches, "Observed ATAC (pop)", "pop"),
    .plt_tsne(results$atac_with_batches, "Observed ATAC (batches)", "batch")
  ) %>% invisible()
}

Bench$rna_atac <- function(velo = T, diff.frac = 0.8, intr.noise = 1, tsne.seed = 0) {
  options_ <- .res_opt(
    rand.seed = 1,
    GRN = GRN_params,
    num.cells = 1000,
    num.cifs = 500,
    tree = Phyla5(),
    diff.cif.fraction = diff.frac,
    intrinsic.noise = intr.noise,
    do.velocity = velo
  )
  
  results <- sim_true_counts(options_)
  # if (intr.noise != 1) {
    return(results)
  # }
  
  add_expr_noise(results)
  divide_batches(results)

  .res_obj(results,
    .plt_tsne(results$counts, "True mRNA TSNE"),
    .plt_tsne(results$atacseq_data, "True ATAC TSNE"),
    .plt_tsne(results$counts_with_batches, "Observed mRNA (pop)", "pop"),
    .plt_tsne(results$counts_with_batches, "Observed mRNA (batches)", "batch"),
    .plt_tsne(results$atac_with_batches, "Observed ATAC (pop)", "pop"),
    .plt_tsne(results$atac_with_batches, "Observed ATAC (batches)", "batch")
  ) %>% invisible()
}


Bench$rna_atac_discrete_3 <- function(tsne.seed = 0) {
  tree10 <- ape::read.tree(text = '(A:1.0,B:1.0,C:1.0,D:1.0:,E:1.0,F:1.0,G:1.0,H:1.0,I:1.0,J:1.0);')
  
  options_ <- .res_opt(
    rand.seed = 1,
    GRN = GRN_params,
    num.cells = 1500,
    num.cifs = 500,
    discrete.cif = T,
    tree = tree10,
    diff.cif.fraction = 0.8,
    do.velocity = F
  )
  
  results <- sim_true_counts(options_)
  add_expr_noise(results)
  divide_batches(results, nbatch = 3)
  
  res <- .res_obj(results,
    .plt_tsne(results$counts, "True mRNA TSNE"),
    .plt_tsne(results$atacseq_data, "True ATAC TSNE"),
    .plt_tsne(results$counts_with_batches, "Observed mRNA (pop)", "pop"),
    .plt_tsne(results$counts_with_batches, "Observed mRNA (batches)", "batch"),
    .plt_tsne(results$atac_with_batches, "Observed ATAC (pop)", "pop"),
    .plt_tsne(results$atac_with_batches, "Observed ATAC (batches)", "batch")
  )
  
  write_files(res, "test_dis_3batches")
  
  res %>% invisible()
}


Bench$rna_velocity <- function(tsne.seed = 0) {
  configs <- expand.grid(
    ngenes = c(100, 200, 500),
    ncells = c(500, 750, 1000),
    seed = 1:8
  ) %>% split(., seq(nrow(.)))
  
  for (conf in configs) {
    c(ngenes, ncells, seed) %<-% conf
    ds_name <- sprintf("grn_%dcells_%dgenes_%d", ncells, ngenes, seed)
    cat(ds_name %+% "\n")
    
    options_ <- .res_opt(
      rand.seed = seed,
      GRN = GRN_params_100,
      # GRN = NA,
      num.genes = ngenes,
      num.cells = ncells,
      num.cifs = 500,
      tree = Phyla3(),
      diff.cif.fraction = 0.8,
      do.velocity = F
    )
    
    results <- sim_true_counts(options_)
    
    write_files(results, ds_name)
    
    # .res_obj(
    #   res = results,
    #   fig = fig
    # ) %>% invisible()
  }
}

Bench$integration <- function(tsne.seed = 0) {
  configs <- expand.grid(
    ngenes = c(100, 200, 500),
    ncells = c(500, 750, 1000),
    seed = 1:2
  ) %>% split(., seq(nrow(.)))
  
  for (conf in configs) {
    c(ngenes, ncells, seed) %<-% conf
    ds_name <- sprintf("dis_%dcells_%dgenes_%d", ncells, ngenes, seed)
    cat(ds_name %+% "\n")
    
    options_ <- .res_opt(
      rand.seed = seed,
      GRN = GRN_params_100,
      # GRN = NA,
      num.genes = ngenes,
      num.cells = ncells,
      num.cifs = 500,
      discrete.cif = T,
      tree = Phyla5(),
      diff.cif.fraction = 0.8,
      do.velocity = F
    )
    
    results <- sim_true_counts(options_)
    
    add_expr_noise(results)
    divide_batches(results)
    
    write_files(results, ds_name)
    
    # .res_obj(
    #   res = results,
    #   fig = fig
    # ) %>% invisible()
  }
}


Bench$cci_60 <- function(seed = 0) {
  set.seed(seed)
  
  n_grn_edge <- 40
  df1 <- data.frame(
    regu = c(1:10, sample(1:10, n_grn_edge, replace = T)),
    tgt = sample(11:30, n_grn_edge + 10, replace = T))
  # assume there's no more than 10 duplicated edges
  df1 <- unique(df1)[1:n_grn_edge,]
  
  # some receptors control GRN regulator
  n_grn_sp_edge <- 12
  df2 <- data.frame(
    regu = sample(31:45, n_grn_sp_edge + 10, replace = T),
    tgt = sample(1:30, n_grn_sp_edge + 10, replace = T))
  # assume there's no more than 10 duplicated edges
  df2 <- unique(df2)[1:n_grn_sp_edge,]

  grn_params <- data.frame(
    target    = c(df1$tgt, df2$tgt),
    regulator = c(df1$regu, df2$regu)
  ) %>% unique()
  
  grn_params$effect <- runif(n_grn_edge + n_grn_sp_edge, min = 2, max = 5)

  lig_params <- data.frame(
    target    = 31:45,
    regulator = 46:60,
    effect    = 5
  )

  options_ <- .res_opt(
    rand.seed = 1,
    threads = 1,
    GRN = grn_params,
    num.genes = 60,
    num.cells = 500,
    num.cifs = 50,
    tree = Phyla1(),
    intrinsic.noise = 1,
    cci = list(
      params = lig_params,
      max.neighbors = 4,
      cell.type.interaction = "random",
      cell.type.lr.pairs = 4:6,
      step.size = 0.2
    ),
  )

  results <- sim_true_counts(options_)
}


Bench$cci <- function(.seed = 0) {
  grn_params <- GRN_params_100 
  
  lig_params <- data.frame(
    target    = c(105, 106, 107, 108),
    regulator = c(101, 102, 103, 104),
    effect    = 5
  )
  
  options_ = .res_opt(
    rand.seed = .seed,
    threads = 1,
    GRN = grn_params,
    num.genes = 110,
    num.cells = 500,
    num.cifs = 100,
    diff.cif.fraction = 0.8,
    tree = Phyla3(),
    do.velocity = F,
    intrinsic.noise = 0.5,
    cci = list(
      params = lig_params,
      max.neighbors = 4
    )
  )
  
  results <- sim_true_counts(options_)
  
  .res_obj(
    res = results,
    grn_fig = plot_gene_module_cor_heatmap(results)
  ) %>% invisible()
}

Bench$cci_propagate <- function(.seed = 0, intr.noise = 1) {
  grn_params <- data.frame(
    target    = 6:20,
    regulator = as.vector(sapply(1:5, \(i) rep(i, 3))),
    effect    = 5 
  ) 
  
  lig_params <- data.frame(
    target    = c( 2,  3,  4,  27, 28),
    regulator = c(21, 22, 23, 24, 25),
    effect    = 5
  )
  
  options_ = .res_opt(
    rand.seed = .seed,
    threads = 1,
    GRN = grn_params,
    num.genes = 30,
    num.cells = 500,
    num.cifs = 500,
    diff.cif.fraction = 0.8,
    tree = Phyla1(),
    do.velocity = F,
    intrinsic.noise = intr.noise,
    cci = list(
      params = lig_params,
      max.neighbors = 2
    )
  )
  
  results <- sim_true_counts(options_)
  
  .res_obj(
    res = results
  )
}


rna_velo_knn <- function(results, velocity, perplexity = 70, randseed = 0, raw = F) {
  set.seed(randseed)
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