.get_results_from_global <- function() {
  get("results", envir = .GlobalEnv)
}


plot_phyla <- function(tree) {
  plotTree(tree, offset = -0.5)
  tiplabels(cex = 2)
  nodelabels(cex = 2)
}


plot_hist <- function(res, log = FALSE, title = NULL) {
  if (is(res, "scMultiSim_output")) {
    res <- res$res
  }
  
  counts <- if (is.environment(res)) {
    res$counts
  } else {
    res
  }
  
  if (log) {
    hist(counts, breaks = 100) 
  } else {
    hist(log(counts + 1), breaks = 100) 
  }
}


#' Plot a (D x cells) data matrix using TSNE.
#'
#' @param labels
#' @param data 
#' @param perplexity 
#' @param label 
#' @param saving 
#' @param plotname 
#' @param rand.seed 
#'
#' @return
#' @export
#'
#' @examples
plot_tsne <- function(data, labels, perplexity = 60, legend = '', plot.name = '', save = F, rand.seed = 0,
                      continuous = F, labels2 = NULL, lim = NULL) {
  set.seed(rand.seed)
  
  data_tsne = Rtsne(t(data), perplexity = perplexity, check_duplicates = FALSE)
  if (!continuous) {
    labels <- factor(labels)
  }
  plot_tsne <- data.frame(
    label = labels,
    x = data_tsne$Y[, 1],
    y = data_tsne$Y[, 2],
    index = seq(nrow(data_tsne$Y))
  )
  p <- ggplot(plot_tsne, aes(x, y, group = index, color = index))
  
  if (is.null(labels2)) {
    p <- p + geom_point(aes(colour = .data[['label']]), shape = 20) + labs(color = legend)
  } else {
    p <- p + geom_point(aes(colour = .data[['label']], shape = factor(labels2))) +
      scale_shape_manual(values=c(4, 15, 5)) +
      labs(color = legend)
  }
  
  if (!is.null(lim))
    p <- p + xlim(lim[1], lim[2]) + ylim(lim[3], lim[4])
  
  if (is.character(save)) {
    save_path <- if (is.character(save)) {
      save
    } else {
      paste0(plot.name, '.pdf')
    }
    pdf(save_path, 5, 5)
    print(p)
    dev.off()
    return(NULL)
  } else {
    p <- p + ggtitle(plot.name)
    return(p)
  }
}


plot_grid <- function(results = .get_results_from_global()) {
  grid <- results$sim$grid
  locs <- sapply(grid$locs, \(a) a)
  data <- data.frame(
    label = "cell",
    x = locs[1,],
    y = locs[2,],
    index = seq(ncol(locs))
  )
  p <- ggplot(data, aes(x, y, group = index, color = index))
  p <- p + geom_point(aes(colour = .data[['label']]), shape = 10) + labs(color = '')
  return(p)
}


.get_gene_module_colors <- function(GRN, N, seed = 0) {
  set.seed(seed)
  
  all_colors <- setNames(rainbow(GRN$n_reg), as.character(GRN$regulators))
  
  gene_module_color_vector <- character(N$gene)
  gene_module_color_vector[(N$grn.gene + 1):N$gene] <- NA
  for (gene_index in 1:num_GRN_genes) {
    if (is.element(gene_index, regulator_ID_list)) {
      gene_module_color_vector[gene_index] <- all_colors[as.character(gene_index)]
    } else {
      gene_module_color_vector[gene_index] <- all_colors[[which.max(gene_effects_by_regulator[gene_index,])]]
    }
  }
  print(all_colors)
  return(gene_module_color_vector)
}


plot_gene_module_cor_heatmap <- function(
  results = .get_results_from_global(),
  seed = 0,
  grn.genes.only = T, save = F
) {
  set.seed(seed)
  grn <- results$.grn$params
  num_genes <- results$num_genes
  counts <- log2(results$counts + 1)
  regulator_ID_list <- sort(unique(grn[, 2]))
  target_gene_ID_list <- sort(unique(grn[, 1]))
  num_target_genes <- length(target_gene_ID_list)
  num_regulators <- length(regulator_ID_list)
  num_GRN_genes <- num_target_genes + num_regulators
  gene_module_color_vector = GetGeneModuleColors(grn, results$.grn$geff, num_genes)
  count_correlation_matrix = GetCountCorrMatrix(counts)
  
  save_path <- if (is.character(save)) {
    save
  } else {
    "GRN_gene_module_correlation_heatmap.pdf"
  }
  
  if (grn.genes.only){
    if (save != F){
      pdf(save_path, 5, 5)
    }
    heatmap.2(count_correlation_matrix[1:num_GRN_genes, 1:num_GRN_genes], scale = "none", Rowv = T, Colv = T, dendrogram = "both", distfun = dist, hclustfun = hclust, key = T, trace = "none", cexRow = 1, cexCol = 1, RowSideColors = gene_module_color_vector[1:num_GRN_genes], ColSideColors = gene_module_color_vector[1:num_GRN_genes], col = bluered(75), main = '          GRN Gene Corr by Main Regulator')
  } else {
    if (save != F){
      pdf(save_path, 5, 5)
    }
    heatmap.2(count_correlation_matrix, scale = "none", Rowv = T, Colv = T, dendrogram = "both", distfun = dist, hclustfun = hclust, key = T, trace = "none", cexRow = 1, cexCol = 1, RowSideColors = gene_module_color_vector, ColSideColors = gene_module_color_vector, col = bluered(75), main = '       Gene Corr by Main Regulator')
  }
  if (save != F) {
    dev.off()
  }
  return()
}


plot_cell_loc <- function(
    results = .get_results_from_global(),
    size = 4, show.label = F, show.arrows = T, lr.pair = 1, .cell.pop = NULL
) {
  if (is.null(.cell.pop))
    .cell.pop <- results$cell_meta$pop
  
  locs <- sapply(results$grid$locs, \(a) a)
  data <- data.frame(
    x = locs[1,],
    y = locs[2,],
    cell_pop = .cell.pop,
    cell_type = results$cell_meta$cell.type
  )
  p <- ggplot()
  p <- p +
    geom_point(aes(x = x, y = y, color = cell_type), data = data, size = size) +
    labs(color = 'Population')
  
  # inter_data <- lapply(1:ncol(locs), \(i) cbind(results$grid$get_neighbours(i, omit.NA = F), i)) %>% 
  #   do.call(rbind, .) %>% na.omit() %>%
  #   apply(1, \(row) c(locs[,row[1]], locs[,row[2]])) %>%
  #   as.matrix() %>% t() %>% as.data.frame()
  # colnames(inter_data) <- c("x", "y", "xend", "yend")
  
  inter_data <- NULL 
  ctp <- results$cci_cell_type_param
  l <- ctp[lr.pair, 1]
  r <- ctp[lr.pair, 2]
  ctp <- ctp[ctp$ligand == l & ctp$receptor == r, ]
  
  if (show.arrows) {
    for (i in 1:ncol(locs)) {
      tp1 <- results$cell_meta$cell.type.idx[i]
      nbs <- results$grid$get_neighbours(i)
      for (nb in nbs) {
        tp2 <- results$cell_meta$cell.type.idx[nb]
        if (any(ctp$ct1 == tp1 & ctp$ct2 == tp2)) {
          inter_data <- rbind(inter_data, c(locs[, i], locs[, nb]))
        }
      }
    }
    inter_data <- as.data.frame(inter_data)
    colnames(inter_data) <- c("x", "y", "xend", "yend")
    to_right <- inter_data$xend > inter_data$x
    to_left <- inter_data$xend < inter_data$x
    to_up <- inter_data$yend > inter_data$y
    to_down <- inter_data$yend < inter_data$y
    inter_data[to_right,]$x <- inter_data[to_right,]$x + 0.2
    inter_data[to_right,]$xend <- inter_data[to_right,]$xend - 0.2
    inter_data[to_left,]$x <- inter_data[to_left,]$x - 0.2
    inter_data[to_left,]$xend <- inter_data[to_left,]$xend + 0.2
    inter_data[to_up,]$y <- inter_data[to_up,]$y + 0.2
    inter_data[to_up,]$yend <- inter_data[to_up,]$yend - 0.2
    inter_data[to_down,]$y <- inter_data[to_down,]$y - 0.2
    inter_data[to_down,]$yend <- inter_data[to_down,]$yend + 0.2
    
    p <- p + geom_segment(
      aes(x = x, y = y, xend = xend, yend = yend),
      data = inter_data,
      arrow = arrow(length = unit(4, "pt"))
    )
  }
  
  if (show.label) {
    p <- p + geom_text(aes(label = as.character(seq_along(x))), size = 2, color = 'black')
  }
  p
}


plot_grn <- function(params) {
  data <- data.frame(
    from = params[,2],
    to = params[,1],
    effect = params[,3]
  )
  ids <- sort(unique(c(data$from, data$to)))
  nodes <- data.frame(id = ids)
  x <- igraph::graph_from_data_frame(data, vertices = nodes)
  igraph::V(x)$color <- ifelse(ids %in% data$from, "#f8766d", "#e5e5e5")
  plot(x,
       vertex.frame.color = "#e5e5e5", vertex.size = 10,
       edge.arrow.size = .4)
}


 gene_corr_regulator <- function(results = .get_results_from_global(), regulator) {
  grn_params <- results$.options$GRN
  regu <- grn_params[grn_params[, 2] == regulator, 1] %>% as.character()
  corr <- gene_corr(results)[regulator,] %>% sort(decreasing = T) %>% round(digits = 2) %>% .[2:length(.)]
  
  cat(sprintf("Gene correlations for %g:\n", regulator))
  for (i in names(corr)) {
    val <- corr[[i]] 
    str <- if (i %in% regu) {
      sprintf("%s(%.2f)", i, val) %>% crayon::bgBlack() %>% crayon::green()
    } else {
      i
    }
    cat(paste0(str, " "))
  }
  
  cat("\n")
}

 
 gene_corr_grn <- function(results = .get_results_from_global()) {
  counts <- log2(results$counts + 1)
  grn_params <- results$.options$GRN
  regulators <- unique(grn_params[, 2])
  total <- 0
  for (i in seq(nrow(grn_params))) {
    rg <- grn_params[i, 2]
    tg <- grn_params[i, 1]
    total <- total + cor(counts[tg,], counts[rg,], method = "spearman")
  }
  total / nrow(grn_params) 
 }


gene_corr <- function(
  results = .get_results_from_global(),
  genes = NULL
) {
  counts <- log2(results$counts + 1)
  genes <- if (is.null(genes)) seq(nrow(counts)) else genes
  ngenes <- length(genes)
   
  res <- matrix(nrow = ngenes, ncol = ngenes) 
  for (i in seq_along(genes)) {
    for (j in 1:i) {
      g1 <- genes[i] 
      g2 <- genes[j]
      res[i, j] <- res[j, i] <- cor(counts[g1,], counts[g2,])
    }
  }
  rownames(res) <- colnames(res) <- genes
  res
}


gene_corr_cci <- function(
  results = .get_results_from_global(),
  all.genes = F,
  .pair = NULL,
  .exclude.same.types = T
) {
  options <- results$.options
  ncells <- options$num.cells
  ngenes <- nrow(results$counts)
  
  if (!("cci" %in% names(options))) {
    warn("CCI is not enabled in the result object")
    return()
  }
  sp_params <- options$cci$params
  
  nb_list <- list()
  non_nb_list <- list()
  for (icell in 1:ncells) {
    nb_list[[icell]] <- nbs <- na.omit(results$grid$get_neighbours(icell))
    non_nb_list[[icell]]  <- sample(setdiff(1:ncells, nbs), size = 4)
  }
  
  if (all.genes) {
    # correlation: CCI regulator - GRN
    
    grn <- options$GRN
    
    regulators <- unique(sp_params[,2])
    cor_list <- lapply(regulators, function(rg) {
      # each target; return cor(rg, all_tg)
      sapply(1:ngenes, function(tg) {
        rg_list <- numeric()
        tg_list <- numeric()
        # each pair
        for (icell in 1:ncells) {
          nbs <- nb_list[[icell]]
          rg_cnt <- results$counts[rg, icell]
          for (nb in nbs) {
            tg_cnt <- results$counts[tg, nb]
            rg_list <- c(rg_list, rg_cnt)
            tg_list <- c(tg_list, tg_cnt)
          }
        }
        cor(rg_list, tg_list)
      }) 
    })
    
    tg_ordered <- lapply(regulators, function(rg) {
      sp_tg <- sp_params[sp_params[,2] == rg, 1]
      grn_tg <- numeric()
      for (grn_rg in sp_tg) {
        grn_tg <- cbind(grn_tg, grn[grn[,2] == grn_rg, 1])
      }
      grn_tg
    }) %>% unlist() %>% unique()
    
    tg_ordered <- c(tg_ordered, setdiff(1:ngenes, tg_ordered)) %>% as.character()
    
    cor_df <- data.frame(tg = character(), rg = character(), cor = numeric())
    
    for (rg in seq_along(cor_list)) {
      cor_df <- rbind(cor_df, data.frame(
        tg = as.character(1:ngenes),
        rg = as.character(regulators[rg]),
        cor = cor_list[[rg]]
      )) 
    }
    
    cor_df$tg <- factor(cor_df$tg, levels = tg_ordered)
    cor_df$cor <- ifelse(cor_df$cor > 0.05, cor_df$cor, 0)
    
    ggplot(cor_df, aes(x=rg, y=tg, fill=cor)) + 
      geom_raster() +
      scale_fill_viridis_c()
  } else {
    # correlation: CCI regulator - CCI target 
    
    target <- numeric()
    regulator <- numeric()
    correlation <- numeric()
    corr.non.nbs <- numeric()
    corr.non.ct <- numeric()
    res_pair <- NULL
    ctp <- results$cci_cell_type_param
    
    for (j in 1:nrow(sp_params)) {
      rg <- sp_params[j, 2]
      tg <- sp_params[j, 1]
      
      rg_list <- numeric()
      tg_list <- numeric()
      nct_tg_list <- numeric()  # cell types without cci
      nct_rg_list <- numeric()
      non_rg_list <- numeric()
      non_tg_list <- numeric()
      
      total <- 0
      for (icell in 1:ncells) {
        nbs <- nb_list[[icell]]
        rg_cnt <- results$counts[rg, icell]
        ct1 <- results$cell_meta$cell.type.idx[icell]
        
        for (nb in nbs) {
          tg_cnt <- results$counts[tg, nb]
          ct2 <- results$cell_meta$cell.type.idx[nb]
          if (.exclude.same.types && ct1 == ct2) next
          if (any(
            ctp$ligand == sprintf("gene%d", rg) &
            ctp$receptor == sprintf("gene%d", tg) &
            ctp$ct1 == ct1 & ctp$ct2 == ct2
          )) {
            rg_list <- c(rg_list, rg_cnt)
            tg_list <- c(tg_list, tg_cnt)
          } else {
            nct_rg_list <- c(nct_rg_list, rg_cnt)
            nct_tg_list <- c(nct_tg_list, tg_cnt)
          }
        }
        non_nbs <- non_nb_list[[icell]]
        for (nb in non_nbs) {
          tg_cnt <- results$counts[tg, nb]
          non_rg_list <- c(non_rg_list, rg_cnt)
          non_tg_list <- c(non_tg_list, tg_cnt)
        }
      }
      
      target <- c(target, tg)
      regulator <- c(regulator, rg)
      correlation <- c(correlation, cor(rg_list, tg_list))
      corr.non.ct <- c(corr.non.ct, cor(nct_rg_list, nct_tg_list))
      corr.non.nbs <- c(corr.non.nbs, cor(non_rg_list, non_tg_list))
      if (!is.null(.pair) && all(.pair == c(rg, tg))) {
        res_pair <- list(rg = rg_list, tg = tg_list)
      }
    }
    
    res <- data.frame(target, regulator, correlation, corr.non.nbs, corr.non.ct)
    
    plot_df <- rbind(
      data.frame(name = "Cell types w/ CCI", value = correlation),
      data.frame(name = "Cell types wo/ CCI", value = corr.non.ct),
      data.frame(name = "Non-nbs", value = corr.non.nbs))
    p <- ggplot(plot_df, aes(x = name, y = value, fill = name)) +
      geom_boxplot() + scale_y_continuous(limits = c(-0.2,0.4))
    
    if (is.null(.pair)) {
      dots_list(res, p)
    } else {
      res_pair 
    }
  }
}


.process_velocity <- function(counts_s, velocity, velocity2 = NULL, perplexity = 70) {
  assertthat::assert_that(
    nrow(counts_s) == nrow(velocity),
    ncol(counts_s) == ncol(velocity)
  )
  
  has_v2 <- !is.null(velocity2)
  if (has_v2) {
    assertthat::assert_that(
      nrow(counts_s) == nrow(velocity2),
      ncol(counts_s) == ncol(velocity2)
    )
  }
  
  future_counts_s <- counts_s + velocity
  future_counts_s[future_counts_s < 0] <- 0
  if (has_v2) {
    future_counts_s2 <- counts_s + velocity2
    future_counts_s2[future_counts_s2 < 0] <- 0
    future_counts_s2 <- t(log2(future_counts_s2 + 1))
  }
  
  counts_s <- t(log2(counts_s + 1))
  future_counts_s <- t(log2(future_counts_s + 1))

  # tsne joint
  combined <- rbind(counts_s, future_counts_s)
  if (has_v2) {
    combined <- rbind(combined, future_counts_s2)
  }
  
  n <- nrow(counts_s)
  data_tsne <- Rtsne(combined, perplexity = perplexity)
  current_counts_s_tsne <- data_tsne$Y[1:n,]
  
  future_1 <- data_tsne$Y[(n + 1):(n * 2),]
  
  get_velo <- function(future_counts_s_tsne) {
    velocity_tsne <- future_counts_s_tsne - current_counts_s_tsne
    vx_raw <- velocity_tsne[, 1]
    vy_raw <- velocity_tsne[, 2]
    
    normalize_velocity <- function(vx, vy) {
      v_normalizer <- (vx^2 + vy^2)^(1 / 2)
      vx <- vx / v_normalizer
      vy <- vy / v_normalizer
      list(vx, vy)
    }
    
    normalize_velocity(vx_raw, vy_raw) %->% c(vx_normalized, vy_normalized)
    
    dist_obj <- dist(current_counts_s_tsne)
    dist_mat <- as.matrix(dist_obj)
    n_cells <- nrow(dist_mat)
    k <- ceiling(n_cells/50)
    
    vx_knn <- distMat.KernelKnn(dist_mat, TEST_indices = NULL, weights_function = 'gaussian', y = vx_raw, k = k, regression = TRUE)
    vy_knn <- distMat.KernelKnn(dist_mat, TEST_indices = NULL, weights_function = 'gaussian', y = vy_raw, k = k, regression = TRUE)
    
    normalize_velocity(vx_knn, vy_knn) %->% c(
      vx_knn_normalized,
      vy_knn_normalized
    )
    
    dots_list(
      current_counts_s_tsne, future_counts_s_tsne,
      vx_raw, vy_raw,
      vx_normalized, vy_normalized,
      vx_knn, vy_knn,
      vx_knn_normalized, vy_knn_normalized,
      .named = T
    ) 
  }
  
  res <- get_velo(future_1)
  
  if (has_v2) {
    future_2 <- data_tsne$Y[(n * 2 + 1):(n * 3),]
    res2 <- get_velo(future_2) 
    
    for (n in names(res2)) {
      if (grepl("current", n, fixed = T)) next
      res[[paste0(n, "2")]] <- res2[[n]]
    }
  }
  
  res
}


plot_rna_velocity <- function(
  results = .get_results_from_global(),
  velocity = results$velocity,
  perplexity = 70, arrow.length = 1, save = F, randseed = 0, ...
) {
  set.seed(randseed)
  counts_s <- results$counts
  cell_pop <- results$cell_meta$pop
  if (is.null(velocity)) {
    stop("The result object is not produced in velocity mode.")
  }
  
  .process_velocity(counts_s, velocity, perplexity = perplexity) %->% c(
    current_counts_s_tsne, future_counts_s_tsne,
    vx_raw, vy_raw,
    vx_normalized, vy_normalized,
    vx_knn, vy_knn,
    vx_knn_normalized, vy_knn_normalized
  ) 
  x1 <- current_counts_s_tsne[, 1]
  y1 <- current_counts_s_tsne[, 2]
  
  args <- default.args(width = 5, height = 5, units = "in", dpi = 1000)
  plot_tsne <- data.frame(x1 = x1, y1 = y1, index = seq(length(x1)), label = cell_pop)
  
  types <- c("raw", "normalized", "knn_normalized")
  titles <- c("Raw Values", "Normalized Values", "KNN Average Normalized Values")
  
  res <- list()
  for (i in seq_along(types)) {
    .type <- types[i]
    .name <- "v_plot_" %+% .type
    .x2e <- "x2_" %+% .type
    .y2e <- "y2_" %+% .type
    
    plot_tsne[[.x2e]] <- x1 + arrow.length * get("vx_" %+% .type)
    plot_tsne[[.y2e]] <- y1 + arrow.length * get("vy_" %+% .type)
    
    .plot <- ggplot(data = plot_tsne, aes(x1, y1, group = index, color = index)) +
      geom_point(aes(colour = .data[['label']])) +
      labs(color = 'pop') +
      geom_segment(aes(x = x1, y = y1, xend = .data[[.x2e]], yend = .data[[.y2e]]),
                   arrow = arrow(length = unit(0.1, "cm")), color = "black", alpha = 1) +
      ggtitle("Velocity - " %+% titles[i])
    
    res[[.type]] <- .plot
    
    if (is.character(save)) {
      ggsave(filename = save %+% .type %+% ".pdf", plot = .plot,
             width = args$width, height = args$height, units = args$units, dpi = args$dpi,
             device = "pdf")
    }
  }
  
  res
} 

# https://github.com/dynverse/dyngen_manuscript/blob/master/package/R/paired_simil.R
paired_simil <- function(x, y, method = "spearman", margin = 1) {
  assertthat::assert_that(
    nrow(x) == nrow(y),
    ncol(x) == ncol(y)
  )
  out <- purrr::map_dbl(
    if (margin == 1) {seq_len(nrow(x))} else {seq_len(ncol(x))},
    function(i) {
      subx <- if (margin == 1) x[i,,drop = FALSE] else x[,i,drop = FALSE]
      suby <- if (margin == 1) y[i,,drop = FALSE] else y[,i,drop = FALSE]
      dynutils::calculate_similarity(subx, suby, method = method, margin = margin)[[1]]
    }
  )
  
  names(out) <- if (margin == 1) rownames(x) else colnames(x)
  
  out
}


cos.sim <- function(a,b) 
{
  dot_product = sum(a*b)
  anorm = sqrt(sum((a)^2))
  bnorm = sqrt(sum((b)^2))
  dot_product / (anorm * bnorm)
}


rna_velocity_cosine <- function(
    results = .get_results_from_global(),
    velocity,
    perplexity = 70, randseed = 0
) {
  set.seed(randseed)
  counts_s <- results$counts
  pop <- results$cell_meta$pop
  depth <- results$cell_meta$depth
  
  if (is.null(results$velocity)) {
    stop("The result object is not produced in velocity mode.")
  }
  
  .process_velocity(counts_s, results$velocity, velocity2 = velocity, perplexity) %->% c(
    current_counts_s_tsne, future_counts_s_tsne,
    vx_raw, vy_raw,
    vx_normalized, vy_normalized,
    vx_knn, vy_knn,
    vx_knn_normalized, vy_knn_normalized,
    future_counts_s_tsne2,
    vx_raw2, vy_raw2,
    vx_normalized2, vy_normalized2,
    vx_knn2, vy_knn2,
    vx_knn_normalized2, vy_knn_normalized2
  ) 
  
  x1 <- current_counts_s_tsne[, 1]
  y1 <- current_counts_s_tsne[, 2]
  
  
  corr <- sapply(seq_along(x1), function(i)
    cos.sim(
      c(vx_knn_normalized[i], vy_knn_normalized[i]),
      c(vx_knn_normalized2[i], vy_knn_normalized2[i])
      ))
  plot_data <- data.frame(
    x = x1, y = y1,
    corr = corr
  )
  p <- ggplot(plot_data, aes(x, y, color = corr)) + geom_point() + scale_colour_gradient2()
  
  return(dots_list(p, corr))
  
  cls_label <- character(length(pop))
  
  for (p in unique(pop)) {
    depth_in_pop <- depth[pop == p]
    depth_min <- min(depth_in_pop)
    depth_max <- max(depth_in_pop)
    span <- (depth_max - depth_min) / 20
    for (i in 1:20) {
      lo <- depth_min + (i - 1) * span
      hi <- depth_min + i * span
      cls_label[pop == p & depth >= lo & depth <= hi] <- paste0(p, "_", i)
    }
  }
  
  labels <- unique(cls_label)
  gt_velo_mean <- numeric()
  res_velo_mean <- numeric()
  for (lb in labels) {
    cells <- cls_label == lb
    locs <- cbind(x1[cells], y1[cells])
    center <- colMeans(locs)
    dists <- sqrt(rowSums(sweep(locs, 2, center)**2))
    weights <- dnorm(dists, 0, 3000)
    weights <- weights / sum(weights)
    gt_velo_mean <- rbind(
      gt_velo_mean,
      colSums(cbind(vx_knn_normalized[cells], vy_knn_normalized[cells]) * weights)
    )
    res_velo_mean <- rbind(
      res_velo_mean,
      colSums(cbind(vx_knn_normalized2[cells], vy_knn_normalized2[cells]) * weights)
    )
  }
  
  sapply(1:nrow(gt_velo_mean), function(i) cos.sim(gt_velo_mean[i,], res_velo_mean[i,]))
  # paired_simil(gt_velo_mean, res_velo_mean, margin = 1, method = "cosine")
}