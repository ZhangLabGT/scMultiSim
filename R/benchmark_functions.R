#' cluster-ability of a dataset 
#' 
#'how well can the cell population be reconstructed by dimensionality reduction methods
#' @param data transcriptomics matrix
#' @param meta information about cells, must contain 'pop' column
#' @param n_pc number of principal components that is fed to tsne
#' @param perplexity perplexity parameter for tsne
#' @param dims number of dimensions kept after tsne
#' @param npop number of populations
#' @param return_kmeans if true returns kmeans information
#' @param pca_scale whether or not to scale pca
#'

ClusterQuality <- function(data, meta, n_pc, perplexity, dims, npop = 5, return_kmeans = F, pca_scale = F) {
  uniqcols <- c(1:length(data[1,]))[!duplicated(t(data))]
  data <- data[, uniqcols]; meta <- meta[uniqcols, , drop = FALSE]
  uniqrows <- c(1:length(data[, 1]))[!duplicated(data)]
  data <- data[uniqrows,]
  data_pc <- prcomp(t(data), scale. = pca_scale)
  data_pc <- data_pc$x[, c(1:n_pc)]
  data_tsne = Rtsne(dims = dims, data_pc, perplexity = perplexity)
  lowdim_data <- data_tsne$Y
  assignment <- kmeans(lowdim_data, npop, nstart = 20)
  ri <- rand.index(assignment[[1]], meta$pop)
  lowdim_dist <- as.matrix(dist(lowdim_data, method = "euclidean"))
  sw <- silhouette(x = meta$pop, dist = dist(lowdim_data, method = "euclidean"))
  sw_summ <- summary(sw)
  sw_mean <- sw_summ[['avg.width']]
  minpop <- sw_summ[["clus.sizes"]] == min(sw_summ[["clus.sizes"]])
  sw_minpop <- sw_summ[["clus.avg.widths"]][minpop]
  # for all cells that are in pop2, what is the most common cluster they are in
  # for all clusters, which one has the highest prop of pop2 cells
  # for the most common pop2 cluster, what is the proportion of pop2 cells? 
  min_in_clst <- sapply(split(x = meta$pop, f = assignment[[1]]), function(X) {
    sum(X %in% c(1:npop)[minpop]) })
  clst_w_min <- sapply(split(x = meta$pop, f = assignment[[1]]), function(X) {
    sum(X %in% c(1:npop)[minpop]) / length(X) })
  res <- c(ri, sw_mean, sw_minpop, min_in_clst, clst_w_min)
  names(res) <- c("RI", "SW_mean", "SW_minpop", paste0("min_in_clst", 1:npop), paste0("clst_w_min", 1:npop))
  if (return_kmeans)
    return(list(measures = res, kmeans_res = assignment))
  else
    return(res)
}

#' rand index
#'
#' compare two clustering result: proportion of pairs of individual that share the same clustering result in both groupings
#' @param group1 first clustering
#' @param group2 second clustering 
rand.index <- function(group1, group2)
{
  x <- c(sapply(group1, function(x) { x == group1 }))
  y <- c(sapply(group2, function(x) { x == group2 }))
  same <- sum(x == y)
  ri <- same / length(x)
  return(ri)
}


#' This function assigns each gene to a gene module based on which regulator gene has the largest effect on the gene
#' @param GRN_params GRN_params is a matrix where: #    - column 1 is the target gene ID, #    - column 2 is the gene ID which acts as a transcription factor for the target (regulated) gene #    - column 3 is the effect of the column 2 gene ID on the column 1 gene ID
#' @param num_genes number of genes
#' @param gene_effects_by_regulator gene effects of each regulator gene (column) on every gene ID (row)
#' @param randseed should produce same result if all other parameters are all the same
GetGeneModuleColors <- function(GRN_params, gene_effects_by_regulator, num_genes, randseed=0) {
  set.seed(randseed)
  regulator_ID_list <- sort(unique(GRN_params[, 2]))
  target_gene_ID_list <- sort(unique(GRN_params[, 1]))
  num_target_genes <- length(target_gene_ID_list)
  num_regulators <- length(regulator_ID_list)
  num_GRN_genes <- num_target_genes + num_regulators

  all_colors <- rainbow(num_regulators)
  names(all_colors) <- as.character(regulator_ID_list)
  gene_module_color_vector <- character(num_genes)
  gene_module_color_vector[(num_GRN_genes + 1):num_genes] <- NA
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

#' This function finds the correlation between every pair of genes
#' @param counts rna seq counts
GetCountCorrMatrix <- function(counts) {
count_correlation_matrix <- cor(t(counts), method = "spearman")
if (any(is.na(count_correlation_matrix))) {
  print('some genes have no counts across all cells; the correlation with these genes will be set to 0 in the heatmap')
  count_correlation_matrix = replace(count_correlation_matrix, which(is.na(count_correlation_matrix)), 0)
}
  return(count_correlation_matrix)
}

#' This function assigns each gene to a gene module based on which regulator gene has the largest effect on the gene
#' @param counts rna seq counts
#' @param GRN_params GRN_params is a matrix where: #    - column 1 is the target gene ID, #    - column 2 is the gene ID which acts as a transcription factor for the target (regulated) gene #    - column 3 is the effect of the column 2 gene ID on the column 1 gene ID
#' @param gene_effects_by_regulator gene effects of each regulator gene (column) on every gene ID (row)
#' @param num_genes number of genes
#' @param saving if the plot should be saved into a file
#' @param GRN_genes_only if the plot shot be generated only for the GRN genes
#' @param randseed should produce same result if all other parameters are all the same
PlotGeneModuleCorrelationHeatmap <- function(counts, GRN_params, gene_effects_by_regulator, num_genes, GRN_genes_only = F, saving = F, randseed=0) {
  set.seed(randseed)
  regulator_ID_list <- sort(unique(GRN_params[, 2]))
  target_gene_ID_list <- sort(unique(GRN_params[, 1]))
  num_target_genes <- length(target_gene_ID_list)
  num_regulators <- length(regulator_ID_list)
  num_GRN_genes <- num_target_genes + num_regulators
  gene_module_color_vector = GetGeneModuleColors(GRN_params, gene_effects_by_regulator, num_genes)
  count_correlation_matrix = GetCountCorrMatrix(counts)
  if (GRN_genes_only){
    if (saving){
      pdf("GRN_gene_module_correlation_heatmap.pdf", 5, 5)
    }
    heatmap.2(count_correlation_matrix[1:num_GRN_genes, 1:num_GRN_genes], scale = "none", Rowv = T, Colv = T, dendrogram = "both", distfun = dist, hclustfun = hclust, key = T, trace = "none", cexRow = 1, cexCol = 1, RowSideColors = gene_module_color_vector[1:num_GRN_genes], ColSideColors = gene_module_color_vector[1:num_GRN_genes], col = bluered(75), main = '          GRN Gene Corr by Main Regulator')
  } else {
    if (saving){
      pdf("gene_module_correlation_heatmap.pdf", 5, 5)
    }
    heatmap.2(count_correlation_matrix, scale = "none", Rowv = T, Colv = T, dendrogram = "both", distfun = dist, hclustfun = hclust, key = T, trace = "none", cexRow = 1, cexCol = 1, RowSideColors = gene_module_color_vector, ColSideColors = gene_module_color_vector, col = bluered(75), main = '       Gene Corr by Main Regulator')
  }
    if (saving){
      dev.off()
    }
}

#' This function gets the average correlation between genes which are most co-regulated by a regulator gene
#' @param counts rna seq counts
#' @param GRN_params GRN_params is a matrix where: #    - column 1 is the target gene ID, #    - column 2 is the gene ID which acts as a transcription factor for the target (regulated) gene #    - column 3 is the effect of the column 2 gene ID on the column 1 gene ID
#' @param gene_effects_by_regulator gene effects of each regulator gene (column) on every gene ID (row)
#' @param num_genes number of genes
GetGeneModuleCorrelation <- function(counts, GRN_params, gene_effects_by_regulator, num_genes) {
  regulator_ID_list <- sort(unique(GRN_params[, 2]))
  target_gene_ID_list <- sort(unique(GRN_params[, 1]))
  num_target_genes <- length(target_gene_ID_list)
  num_regulators <- length(regulator_ID_list)
  num_GRN_genes <- num_target_genes + num_regulators
  cor_gene_module_pairs <- numeric()
  gene_module_color_vector = GetGeneModuleColors(GRN_params, gene_effects_by_regulator, num_genes)
  is_NA = F
  for (module_gene_index1 in seq_along(gene_module_color_vector[1:num_GRN_genes])) {
    for (module_gene_index2 in seq_along(gene_module_color_vector[1:num_GRN_genes])) {
      if ((module_gene_index1 != module_gene_index2) && (gene_module_color_vector[module_gene_index1] == gene_module_color_vector[module_gene_index2])) {
        cor_gene_module_pairs_next = suppressWarnings(cor(counts[module_gene_index1,], counts[module_gene_index2,], method = "spearman"))
        if (!is.na(cor_gene_module_pairs_next)) {
          cor_gene_module_pairs <- c(cor_gene_module_pairs, cor_gene_module_pairs_next)
        } else {
          is_NA = T
        }
      }
    }
  }
  if (is_NA == T){
    print('some genes have no counts across all cells; the correlation with these genes is ignored in the mean within-gene-module correlation calculation')
  }
  gene_module_correlation <- mean(cor_gene_module_pairs)
  return(gene_module_correlation)
}

#' This function gets the average correlation rna seq counts and region effect on genes for genes which are only associated with 1 chromatin region
#' @param counts rna seq counts
#' @param atacseq_data atac seq data
#' @param region2gene a 0 1 coupling matrix between regions and genes of shape (nregions) x (num_genes), where a value of 1 indicates the gene is affected by a particular region
#' @export
Get_1region_ATAC_correlation <- function(counts, atacseq_data, region2gene) {
  target_genes <- which(colSums(region2gene > 0) == 1)
  ATAC_1region_correlation <- numeric()
  for (gene_index in target_genes) {
    region <- which(region2gene[, gene_index] > 0)
    correlation <- suppressWarnings(cor(atacseq_data[region,], counts[gene_index,], method = "spearman"))
    ATAC_1region_correlation <- c(ATAC_1region_correlation, correlation)
  }
  ATAC_1region_correlation <- mean(ATAC_1region_correlation, na.rm = T)
  return(ATAC_1region_correlation)
}

#' This function gets the average correlation rna seq counts and chromatin region effect on genes
#' @param counts rna seq counts
#' @param atacseq_data atac seq data
#' @param num_genes number of genes
#' @export
Get_ATAC_correlation <- function(counts, atacseq_data, num_genes) {
  ATAC_correlation <- numeric()
  for (gene_index in 1:num_genes) {
    ATAC_correlation <- c(ATAC_correlation, suppressWarnings(cor(atacseq_data[gene_index,], counts[gene_index,], method = "spearman")))
  }
  ATAC_correlation <- mean(ATAC_correlation, na.rm = T)
  return(ATAC_correlation)
}

#' Perform tSNE 
#' 
#' This function takes the ground truth (meta), the expresssion matrix (data), and performs tSNE, saves the result to a jpeg file with user specified name
#' @param meta simulation parameters
#' @param data expression matrix
#' @param plotname the name of the jpeg file
#' @param label the column name of the meta data that the points needs to be colored by
#' @param perplexity the perplexity parameter used in the tsne visualization
#' @param saving whether or not to save the plot
#' @param randseed should produce same result if all other parameters are all the same
#' @import Rtsne
#' @import ggplot2
PlotTsne <- function(meta, data, perplexity = 60, label, saving = T, plotname, randseed=0) {
  set.seed(randseed)
  data_tsne = Rtsne(t(data), perplexity = perplexity, check_duplicates = FALSE)
  plot_tsne <- data.frame(label = factor(meta), x = data_tsne$Y[, 1], y = data_tsne$Y[, 2], index = seq(nrow(data_tsne$Y)))
  p <- ggplot(plot_tsne, aes(x, y, group = index, color = index))

  p <- p + geom_point(aes(colour = .data[['label']]), shape = 20) + labs(color = label)

  if (saving == T) { ggsave(p, filename = paste0(plotname, '.pdf'), device = 'pdf', width = 5, height = 5) }
  if (saving == F) { p <- p + ggtitle(plotname) }
  return(list(plot_tsne, p))
}

#' Perform PCA  
#' 
#' This function takes the ground truth (meta), the expresssion matrix (data), and performs tSNE, saves the result to a jpeg file with user specified name
#' @param meta simulation parameters
#' @param data expression matrix
#' @param plotname the name of the jpeg file
#' @param label the column name of the meta data that the points needs to be colored by
#' @param discrete whether or not the population is discrete versus continuous
#' @param saving whether or not to save the plot
PlotPCA <- function(meta, data, plotname, label, discrete = T, saving = F) {
  uniqcols <- c(1:length(data[1,]))[!duplicated(t(data))]
  data <- data[, uniqcols]; meta <- meta[uniqcols,]
  uniqrows <- c(1:length(data[, 1]))[!duplicated(data)]
  data <- data[uniqrows,]
  data_pc = prcomp(t(data))
  if (discrete == T) {
    plot_pca = data.frame(meta, label = factor(meta[, label]), x = data_pc$x[, 1], y = data_pc$x[, 2])
  }else {
    plot_pca = data.frame(meta, label = meta[, label], x = data_pc$x[, 1], y = data_pc$x[, 2])
  }
  p <- ggplot(plot_pca, aes(x, y))
  p <- p + geom_point()
  p <- p +
    geom_point(aes(colour = plot_pca[['label']])) +
    labs(color = label)
  if (saving == T) { ggsave(p, filename = plotname, device = 'jpeg', width = 5, height = 4) }
  return(list(plot_pca, p))
}


#' calculate area under curve
#' @param x_vec vector of x values
#' @param y_vec vector of y values
cal_AUC <- function(x_vec, y_vec) {
  xtemp <- x_vec[2:length(x_vec)] - x_vec[1:(length(x_vec) - 1)]
  ytemp <- y_vec[2:length(y_vec)] - y_vec[1:(length(y_vec) - 1)]
  return(sum(xtemp * (y_vec[1:(length(y_vec) - 1)] + ytemp / 2)))
}

#' calculate sensitivity and specificity
#' @param isDE_gold true values
#' @param isDE_pred predicted values
sens_and_spec <- function(isDE_gold, isDE_pred) {
  TP <- sum(isDE_gold & isDE_pred, na.rm = T)
  TN <- sum(!isDE_gold & !isDE_pred, na.rm = T)
  FP <- sum(!isDE_gold & isDE_pred, na.rm = T)
  FN <- sum(isDE_gold & !isDE_pred, na.rm = T)
  sensi <- TP / (TP + FN); speci <- TN / (TN + FP)
  return(c(sensi, speci))
}

#' convert a decimal number to k based number
#' 
#' @param k the base number
#' @param nbits the length (number of bits) in the target string
#' @param dec_number the decimal number to convert
#' @return a vector with values as the converted string at each position
Dec2k_based <- function(k, nbits, dec_number) {
  res_code <- numeric(nbits)
  temp <- dec_number
  for (iexp in seq((nbits - 1), 0, -1)) {
    ratio <- temp / (k^iexp)
    if (ratio >= 1) {
      res_code[nbits - iexp] <- floor(ratio)
      temp <- temp - floor(ratio) * (k^iexp)
    } else { res_code[nbits - iexp] <- 0 }
  }
  return(res_code)
}

#' convert a k based number to decimal number
#' @param k number base
#' @param basek_vec input vector to convert
basek2decimal <- function(k, basek_vec) {
  nbits <- length(basek_vec)
  tosum <- numeric(nbits)
  for (ibit in 1:nbits) {
    tosum[ibit] <- basek_vec[ibit] * (k^(nbits - ibit))
  }
  return(sum(tosum))
}

#' calculate adjusted rand index between labels from a clustering algorithm and known labels
#' @param cluster_res results from a cluster algorithm
#' @param label known label of the cells
cal_ARI <- function(cluster_res, label) {
  ri_all <- adj.rand.index(cluster_res, label)
  ri_pop <- sapply(sort(unique(label)), function(i) {
    adj.rand.index(cluster_res, label == i)
  })
  ri <- c(ri_all, ri_pop)
  return(ri)
}

KNNP <- function(data, label, mutual = T) {
  knn <- get.knn(data, k = 20)
  if (mutual == T) {
    knn <- lapply(c(1:length(data[, 1])), function(i) {
      X = knn[[1]][i,]
      mut = sapply(X, function(x) {
        i %in% knn[[1]][x,]
      })
      return(X[mut])
    })
  }else {
    knn <- lapply(c(1:length(data[, 1])), function(i) {
      X = knn[[1]][i,]
      return(X) })
  }
  purity <- sapply(c(1:length(knn)), function(i) {
    correct <- label[i] == label[knn[[i]]]
    return(sum(correct) / length(correct))
  })
  return(purity)
}

#' calculate knn purity of a data space with given labels
#' @param data input data
#' @param label known label of the cells
cal_KNNP <- function(data, label) {
  knn_purity <- KNNP(data, label)
  kp_all <- mean(knn_purity, na.rm = T)
  kp_pop <- sapply(sort(unique(label)), function(i) {
    mean(knn_purity[label == i], na.rm = T)
  })
  kp <- c(kp_all, kp_pop)
  return(kp)
}

#' calculate pseudotime correlation between our true pseudotime and predicted pseudotime
#' true pseudotime is a list where each element corresponds to an lineage
#' @param true_time true pseudotime
#' @param pred_time predicted pseudotime
cor_pseudotime <- function(true_time, pred_time) {
  cor_vec <- sapply(1:length(true_time), function(i) {
    return(cor(pred_time[paste0("cell", true_time[[i]][, 1])], true_time[[i]][, 2], method = "spearman"))
  })
  return(mean(abs(cor_vec)))
}
