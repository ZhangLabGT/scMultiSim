#' Visualize the simulated velocity, colored using the simulated time
#' @param counts_s spliced rna counts
#' @param velocity true rna velocity
#' @param cell_time cell pseudotime
#' @param cell_pop cell population metadata
#' @param width width of the figure
#' @param height height of the figure
#' @param units select from "cm", "mm" and "in", default "in"
#' @param dpi dpi value of the figure
#' @param arrow.length number of time units into the future to calculate future cells based on velocity
#' @param perplexity perplexity used in tsne visualization
#' @param saving if the plot should be saved into a file
#' @param randseed should produce same result if all other parameters are all the same
#' @return None
#' @import ggplot2 KernelKnn
plotVelo <- function(counts_s, velocity, cell_time, cell_pop, width = 20*2/3, height = 10, units = "in", dpi = 1000, perplexity = 70, arrow.length = 1, saving = F, randseed=0) {
  set.seed(randseed)
  future_counts_s <- counts_s + velocity
  future_counts_s[which(future_counts_s < 0)] <- 0

  counts_s = t(log2(counts_s + 1))
  future_counts_s = t(log2(future_counts_s + 1))

  # tsne joint
  counts_s_current_and_future <- list(counts_s, future_counts_s)
  counts_s_current_and_future = do.call(rbind, counts_s_current_and_future)

  data_tsne = Rtsne(counts_s_current_and_future, perplexity = perplexity)
  current_counts_s_tsne = data_tsne$Y[1:(nrow(data_tsne$Y) / 2),]
  future_counts_s_tsne = data_tsne$Y[(nrow(data_tsne$Y) / 2 + 1):nrow(data_tsne$Y),]

  x1 = current_counts_s_tsne[, 1]
  y1 = current_counts_s_tsne[, 2]

  velocity_tsne = future_counts_s_tsne - current_counts_s_tsne
  vx_raw = velocity_tsne[, 1]
  vy_raw = velocity_tsne[, 2]

  normalize_velocity <- function(vx, vy) {
    v_normalizer = (vx^2 + vy^2)^(1 / 2)
    vx = vx / v_normalizer
    vy = vy / v_normalizer
    return(list(vx, vy))
  }

  v_raw = normalize_velocity(vx_raw, vy_raw)
  vx_normalized = v_raw[[1]]
  vy_normalized = v_raw[[2]]

  dist_obj = dist(current_counts_s_tsne)
  dist_mat = as.matrix(dist_obj)
  n_cells = nrow(dist_mat)
  k = ceiling(n_cells/50)

  vx_knn = distMat.KernelKnn(dist_mat, TEST_indices = NULL, weights_function = 'gaussian', y = vx_raw, k = k, regression = TRUE)
  vy_knn = distMat.KernelKnn(dist_mat, TEST_indices = NULL, weights_function = 'gaussian', y = vy_raw, k = k, regression = TRUE)

  v_knn_normalized = normalize_velocity(vx_knn, vy_knn)
  vx_knn_normalized = v_knn_normalized[[1]]
  vy_knn_normalized = v_knn_normalized[[2]]

  x2_raw = x1 + arrow.length * vx_raw
  y2_raw = y1 + arrow.length * vy_raw

  x2_normalized = x1 + arrow.length * vx_normalized
  y2_normalized = y1 + arrow.length * vy_normalized

  x2_knn_normalized = x1 + arrow.length * vx_knn_normalized
  y2_knn_normalized = y1 + arrow.length * vy_knn_normalized

  index = seq(nrow(data_tsne$Y) / 2)

  plot_tsne <- data.frame(x1 = x1, y1 = y1, x2_raw = x2_raw, y2_raw = y2_raw, x2_normalized = x2_normalized, y2_normalized = y2_normalized, x2_knn_normalized = x2_knn_normalized, y2_knn_normalized = y2_knn_normalized, index = index, label = cell_pop)

  v_plot_raw <- ggplot(data = plot_tsne, aes(x1, y1, group = index, color = index)) +
    geom_point(aes(colour = .data[['label']])) +
    labs(color = 'pop') +
    geom_segment(aes(x = x1, y = y1, xend = x2_raw, yend = y2_raw), arrow = arrow(length = unit(0.1, "cm")), color = "black", alpha = 1) +
    ggtitle('Velocity - Raw Values')

  v_plot_normalized <- ggplot(data = plot_tsne, aes(x1, y1, group = index, color = index)) +
    geom_point(aes(colour = .data[['label']])) +
    labs(color = 'pop') +
    geom_segment(aes(x = x1, y = y1, xend = x2_normalized, yend = y2_normalized), arrow = arrow(length = unit(0.1, "cm")), color = "black", alpha = 1) +
    ggtitle('Velocity - Values Normalized ')

  v_plot_knn_normalized <- ggplot(data = plot_tsne, aes(x1, y1, group = index, color = index)) +
    geom_point(aes(colour = .data[['label']])) +
    labs(color = 'pop') +
    geom_segment(aes(x = x1, y = y1, xend = x2_knn_normalized, yend = y2_knn_normalized), arrow = arrow(length = unit(0.1, "cm")), color = "black", alpha = 1) +
    ggtitle('Velocity - KNN Average Values Normalized ')

  if (saving == T){
    ggsave(filename = 'velocity raw.pdf', plot = v_plot_raw, device = "pdf", path = "./", height = height, width = width, units = units, dpi = dpi)
    ggsave(filename = 'velocity normalized.pdf', plot = v_plot_normalized, device = "pdf", path = "./", height = height, width = width, units = units, dpi = dpi)
    ggsave(filename = 'velocity knn normalized.pdf', plot = v_plot_knn_normalized, device = "pdf", path = "./", height = height, width = width, units = units, dpi = dpi)
  }
  return(list(v_plot_raw = v_plot_raw, v_plot_normalized = v_plot_normalized, v_plot_knn_normalized = v_plot_knn_normalized))
}


#' Calculate the coefficient of variation of a vector x
#' @param x input vector x
cv <- function(x) { return(sd(x) / mean(x)) }

#' Calculate the fano factor of a vector x
#' @param x input vector x
fano <- function(x) { return((sd(x))^2 / mean(x)) }

#' Calculate the percentage of non-zero values in a vector x
#' @param x input vector x
percent_nonzero <- function(x) { return(sum(x > 0) / length(x)) }


#' how well do dots fit to the diagonal (for evaluating our qqplots)
#' @param x a vector on x axis
#' @param y a vector on y axis with the same length as x
#' @param nbins number of bins to use
dist2diag <- function(x, y, nbins) {
  min_val <- min(c(x, y)); max_val <- max(c(x, y))
  # make nbins equal interval bins
  bin_width <- (max_val - min_val) / nbins
  binned_xy <- sapply(1:nbins, function(ibin) {
    bin_min <- min_val + (ibin - 1) * bin_width
    bin_max <- bin_min + bin_width
    in_bin <- which(x >= bin_min & x <= bin_max)
    return(c(mean(x[in_bin]), mean(y[in_bin])))
  })
  return(mean(abs(binned_xy[1,] - binned_xy[2,]), na.rm = T))
}


#' Get the best matched parameter
#'
#' This function matches a real dataset to a database of summary information of simulated datasets, plots a qqplot for user to visualize similarity between their dataset and the simulated dataset, and suggests parameters to use in the simulation
#' @param tech 'nonUMI','UMI' the match database are constructed based on reasonable values for each technology
#' @param counts expression matrix
#' @param plotfilename output name for qqplot
#' @param n_optimal number of top parameter configurations to return
#' @param depth_range if one knows the rough range of sequencing depth, it can be input here.
#' @param alpha_range if one knows the rough range of mRNA capture efficiency, it can be input here.
#' @param idx_set set of indices
#' @return three set of best matching parameters that was used to simulate the best matching dataset to the experimental dataset
BestMatchParams <- function(tech, counts, plotfilename, n_optimal = 3,
                            depth_range = c(-Inf, Inf), alpha_range = c(-Inf, Inf), idx_set = NULL) {
  #counts <- counts[rowSums(counts>0)>10, ]
  mean_exprs <- quantile(rowMeans(counts + 1, na.rm = T), seq(0, 1, 0.002))
  sd_exprs <- quantile(apply(counts, 1, sd), seq(0, 1, 0.002), na.rm = T)
  percent0 <- quantile(apply(counts, 1, percent_nonzero), seq(0, 1, 0.002))

  tempdata <- read.table(system.file("extdata/grid_summary", sprintf("mean_bins_%s.txt", tech), package = "SymSim"),
                         stringsAsFactors = F)
  mean_bins <- unname(as.matrix(tempdata))
  tempdata <- read.table(system.file("extdata/grid_summary", sprintf("sd_bins_%s.txt", tech), package = "SymSim"),
                         stringsAsFactors = F)
  sd_bins <- unname(as.matrix(tempdata))
  tempdata <- read.table(system.file("extdata/grid_summary", sprintf("nonzero_bins_%s.txt", tech), package = "SymSim"),
                         stringsAsFactors = F)
  nonzero_bins <- unname(as.matrix(tempdata))
  load(system.file("extdata/grid_summary", sprintf("sim_params_%s.RData", tech), package = "SymSim"))
  chosen_params <- which(sim_params$depth_mean >= depth_range[1] &
                           sim_params$depth_mean <= depth_range[2] &
                           sim_params$alpha_mean >= alpha_range[1] &
                           sim_params$alpha_mean <= alpha_range[2])
  if (!is.null(idx_set)) {
    chosen_params <- intersect(chosen_params, idx_set)
  }
  if (length(chosen_params) > 0) {
    grid_summary <- list(mean_bins[chosen_params,], nonzero_bins[chosen_params,], sd_bins[chosen_params,])
  } else { stop("Error in the depth range") }
  exp_summary <- list(mean_exprs, percent0, sd_exprs)

  dists <- lapply(c(1:3), function(i) {
    if (i %in% c(1, 3)) {
      dist <- apply(grid_summary[[i]], 1, function(X) {
        mean(abs(log10(X) - log10(exp_summary[[i]]))) + dist2diag(log10(X), log10(exp_summary[[i]]), nbins = 20)
      })
    } else {
      dist <- apply(grid_summary[[i]], 1, function(X) {
        mean(abs(X - exp_summary[[i]])) + dist2diag(X, exp_summary[[i]], nbins = 20)
      })
    }
    return(dist)
  })

  dists <- do.call(cbind, dists)
  dists <- rowSums(dists)
  sorted_dists <- sort.int(dists, index.return = T)
  best_match <- sorted_dists$ix[1:n_optimal]

  best_params <- lapply(chosen_params[best_match], function(X) { sim_params[X,] })
  plotnames <- c('log10(mean)', 'percent_nonzero', 'log10(sd)')

  if (!is.na(plotfilename)) {
    pdf(file = sprintf("%s.pdf", plotfilename), 10, 23)
    par(mfrow = c(6, 3))
    for (i in c(1:n_optimal)) {
      for (k in c(1:3)) {
        bin1 = grid_summary[[k]][best_match[i],]
        bin2 = exp_summary[[k]]
        if (k %in% c(1, 3)) { bin1 <- log(base = 10, bin1); bin2 <- log(base = 10, bin2) }
        plot(bin1, bin2, pch = 16, xlab = 'simulated values', ylab = 'experimental values', main = paste('No.', i, 'best', 'match', plotnames[k]))
        lines(c(-10, 10), c(-10, 10), col = 'red')
      }
    }
    dev.off()
  }

  par(mfrow = c(1, 3))
  for (k in c(1:3)) {
    bin1 = grid_summary[[k]][best_match[1],]
    bin2 = exp_summary[[k]]
    if (k %in% c(1, 3)) { bin1 <- log(base = 10, bin1); bin2 <- log(base = 10, bin2) }
    plot(bin1, bin2, pch = 16, xlab = 'simulated values', ylab = 'experimental values', main = plotnames[k])
    lines(c(-10, 10), c(-10, 10), col = 'red')
  }
  best_params <- do.call(rbind, best_params)
  # best_params <- best_params[, c("gene_effects_sd", "gene_effect_prob", "nevf", "Sigma",
  #                                "alpha_mean", "alpha_sd", "depth_mean", "depth_sd")]
  best_params$dist <- sorted_dists$x[1:n_optimal]
  return(best_params)
}

#' Getting logged expression distribution
#'
#' Prepares for plotting the Count Heatmap
#' @param counts the expression matrix
#' @param log_count_bins a vector of the cut-offs for the histogram
#' @return a matrix where the rows are the genes and columns are the number of samples within a count category
LogDist <- function(counts, log_count_bins) {
  log_dist = apply(log(counts + 1, base = 10), 1, function(x) {
    if (sum(is.na(x)) != length(x)) {
      dist0 = sum(x == 0)
      range_c = log_count_bins
      count = table(cut(x[x > 0], range_c))
      dist = c(dist0, count) / (sum(count) + dist0)
      return(dist) }else {
      return(NA)
    }
  })
  log_dist = t(log_dist)
  return(log_dist)
}


#' Plotting logged expression distribution
#'
#' takes an expression matrix and makes a 2D histogram on the log scale where each row is a gene and the number of samples in a bin is shown by the intensity of the color
#' @param log_real the logged distribution of count distribution, obtained through the LogDist function
#' @param mean_counts the average expression for each gene, used for sorting purpose
#' @param given_ord the given order of genes
#' @param zeropropthres the genes with zeroproportion greater than this number is not plotted (default to 0.8)
#' @param filename the name of the output plot. Will not be used if saving=F.
#' @param saving if the plot should be saved into a file
#' @param data_name a string which is included in the title of the plot to describe the data used
PlotCountHeatmap <- function(log_real, mean_counts, data_name, given_ord = NA, zeropropthres = 0.8,
                             filename = NA, saving = F) {
  mean_counts = mean_counts[log_real[, 1] < zeropropthres]
  log_real = log_real[log_real[, 1] < zeropropthres,]
  colnames(log_real)[1] = '.0'
  plot_real = melt(log_real)
  plot_real$freq = plot_real$value
  genenames = rownames(log_real)
  if (is.null(genenames)) {
    genenames = as.character(c(1:length(log_real[, 1])))
    rownames(log_real) = genenames
  }
  if (is.na(given_ord[1])) {
    cluster_ord <- hclust(dist(log_real, method = "euclidean"), method = "ward.D")$order
    ord = order(cut(log(mean_counts + 1), 30), order(cluster_ord))
  }else { ord <- given_ord }
  plot_real$Gene <- factor(plot_real$X1, levels = genenames[ord])
  p <- ggplot(plot_real, aes(X2, Gene)) +
    geom_tile(aes(fill = freq)) +
    scale_fill_gradient(low = "white", high = "black", trans = 'identity') +
    ggtitle(sprintf('distribution of mRNA counts of %s', data_name)) +
    labs(colour = 'Percentage of Cells', x = 'log10(Count) bins', y = 'Genes') +
    scale_y_discrete(breaks = NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if (saving == T) { ggsave(filename, dev = 'jpeg', width = 8, height = 8) }else { return(list(ord, p)) }
}


#' Plotting the histograms of kon,koff,s values
#'
#' plot colored histograms of parameters
#' @param params a matrix of 3 columns, the first one is kon, the second is koff and the third is s
#' @param samplename the prefix of the plot, the suffix is '.params_dist.jpeg'
#' @param saving if the plot should be saved to file
#' @return make a plot of three histograms
PlotParamHist<-function(params,samplename,saving=F){
  df <- data.frame(kon = log(base=10,params[,1]),koff=log(base=10,params[,2]),s=log(base=10,params[,3]))
  df <- melt(df)
  p1 <- ggplot(df,aes(x=value)) +
    geom_histogram(data=subset(df,variable == 'kon'),aes(y = ..density..), binwidth=density(df$value)$bw) +
    geom_density(data=subset(df,variable == 'kon'),fill="red", alpha = 0.2)
  p2 <- ggplot(df,aes(x=value)) +
    geom_histogram(data=subset(df,variable == 'koff'),aes(y = ..density..), binwidth=density(df$value)$bw) +
    geom_density(data=subset(df,variable == 'koff'),fill="green", alpha = 0.2)
  p3 <- ggplot(df,aes(x=value)) +
    geom_histogram(data=subset(df,variable == 's'),aes(y = ..density..), binwidth=density(df$value)$bw) +
    geom_density(data=subset(df,variable == 's'),fill="blue", alpha = 0.2)
  if(saving==T){ggsave(paste(samplename,'.params_dist.jpeg',sep=''),plot=arrangeGrob(p1, p2, p3, ncol=1),device='jpeg')}else{p <- arrangeGrob(p1, p2, p3, ncol=1)}
  return(p)
}


#' rescale2range
#'
#' Subfunction for Plotting FNR. Rescale the values in vec such that the lagest is n, and the smallest is 1.
#' @param vec input vector
#' @param n the largest integer to scale the vector to (the smallest is 1)
rescale2range <- function(vec, n) {
  a <- (n - 1) / (max(vec) - min(vec))
  return(a * vec + (1 - a * min(vec)))
}

#' Plotting PCA results (PC1 and PC2)
#' @param PCAres the PCA results
#' @param col_vec a vector to specify the colors for each point
#' @param figuretitle title for the plot
plotPCAbasic <- function(PCAres, col_vec, figuretitle) {
  variance_perc <- 100 * (PCAres$sdev)^2 / sum((PCAres$sdev)^2)
  plot(PCAres$x[, 1], PCAres$x[, 2], col = col_vec, pch = 20,
       xlab = sprintf("PC1 %4.2f%%", variance_perc[1]),
       ylab = sprintf("PC2 %4.2f%%", variance_perc[2]),
       main = figuretitle)
}

#' arrange multiple plots from ggplot2 to one figure.
#' From http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#' ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
#' @param cols Number of columns in layout
#' @param layout A matrix specifying the layout. If present, 'cols' is ignored.
#' @param ... plots to make a list of
#' @param plotlist a list of plots
multiplot <- function(..., plotlist = NULL, cols = 1, layout = NULL) {

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots / cols)),
                     ncol = cols, nrow = ceiling(numPlots / cols))
  }

  if (numPlots == 1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#' plot the histogram of number of reads per UMI
#' @param hist_res contains information for creating histogram
#' @param plot_title plot title
plot_nreads_UMI <- function(hist_res, plot_title) {
  toplot <- hist_res$counts / sum(hist_res$counts[2:length(hist_res$counts)]) * 100;
  names(toplot) <- hist_res$mids
  xpos <- hist_res$mids[2:length(hist_res$mids)]
  plotres <- plot(hist_res$mids, toplot, log = "x", col = adjustcolor("white", alpha.f = 1),
                  xlab = "Reads/molecule", ylab = "Fraction (%)", main = plot_title, ylim = c(1, max(toplot[2:length(toplot)])))
  rect(xleft = xpos - 0.5, ybottom = 0, xright = xpos + 0.5, ytop = toplot[2:length(toplot)], col = gray(0.5), border = NA)
}


#' function QQplot2RealData; plot both with ggplot2 and with basic plot
#' @param real_data gene expression matrix of experimental data
#' @param sim_data gene expression matrix of simulated data
#' @param expressed_prop the minimum proportion of cells a gene should be expressed in
#' @param plot_file_name the pdf file name for the output plots
#' @param print_data_dir the directory to `print` the source data for plots; if NA then do not print
#' @import ggplot2
QQplot2RealData <- function(real_data, sim_data, expressed_prop, plot_file_name, print_data_dir = NA) {
  real_data <- real_data[rowSums(real_data > 0) > floor(dim(real_data)[2] * expressed_prop),]
  sim_data <- sim_data[rowSums(sim_data > 0) > floor(dim(sim_data)[2] * expressed_prop),]

  mean_sim <- quantile(rowMeans(sim_data + 1, na.rm = T), seq(0, 1, 0.002))
  sd_sim <- quantile(apply(sim_data, 1, sd), seq(0, 1, 0.002), na.rm = T)
  nonzero_sim <- quantile(apply(sim_data, 1, percent_nonzero), seq(0, 1, 0.002))

  mean_real <- quantile(rowMeans(real_data + 1, na.rm = T), seq(0, 1, 0.002))
  sd_real <- quantile(apply(real_data, 1, sd), seq(0, 1, 0.002), na.rm = T)
  nonzero_real <- quantile(apply(real_data, 1, percent_nonzero), seq(0, 1, 0.002))

  sim_summary <- list(mean_sim, nonzero_sim, sd_sim)
  real_summary <- list(mean_real, nonzero_real, sd_real)
  d_data <- lapply(c(1:3), function(k) {
    bin1 = sim_summary[[k]]
    bin2 = real_summary[[k]]
    if (k %in% c(1, 3)) {
      bin1 = log(base = 10, bin1); bin2 = log(base = 10, bin2)
    }
    data.frame(simulation = bin1, experimental = bin2, summary = c('Mean', 'Percent Nonzero', 'SD')[k], tech = 'nonUMI')
  })
  p_d = lapply(d_data, function(X) {
    ggplot(X, aes(x = simulation, y = experimental)) +
      labs(x = "", y = "") +
      geom_point(alpha = 0.6) +
      geom_abline(intercept = 0, slope = 1, col = 'red')
  })

  pdf(plot_file_name, 10.5, 3.5)
  par(mfrow = c(1, 3))
  plot(log10(mean_sim), log10(mean_real), pch = 16, xlab = 'simulated data', ylab = 'real data', main = "log10(mean+1)")
  abline(a = 0, b = 1, col = "red")
  plot(nonzero_sim, nonzero_real, pch = 16, xlab = 'simulated data', ylab = 'real data', main = "percent_nonzero")
  abline(a = 0, b = 1, col = "red")
  plot(log10(sd_sim), log10(sd_real), pch = 16, xlab = 'simulated data', ylab = 'real data', main = "log10(sd)")
  abline(a = 0, b = 1, col = "red")
  par(mfrow = c(1, 1))
  multiplot(
    p_d[[1]] +
      labs(x = "", y = "") +
      theme(text = element_text(size = 15)),
    p_d[[2]] +
      labs(x = "", y = "") +
      theme(text = element_text(size = 15)),
    p_d[[3]] +
      labs(x = "", y = "") +
      theme(text = element_text(size = 15)), cols = 3)
  dev.off()
  print(sprintf("QQ-plots printed into file %s.", plot_file_name))
  if (!is.na(print_data_dir)) {
    write.table(data.frame(simulated = log10(mean_sim), experimental = log10(mean_real)),
                sprintf("%s/log10meanplus1.txt", print_data_dir), quote = F, row.names = F)
    write.table(data.frame(simulated = nonzero_sim, experimental = nonzero_real),
                sprintf("%s/percent_nonzero.txt", print_data_dir), quote = F, row.names = F)
    write.table(data.frame(simulated = log10(sd_sim), experimental = log10(sd_real)),
                sprintf("%s/log10sd.txt", print_data_dir), quote = F, row.names = F)
  }
  return(NULL)
}


#' retrieve genes' differential expression information
#' it outputs three measures for every gene:
#' nDiffEVF: the number of DiffEVFs used for each gene
#' logFC_theoretical: log2 fold change based on kinetic parameters
#' wil.p_true_counts: adjusted wilcoxon p-value based on true counts
#' @param true_counts_res the output of function SimulateTrueCounts()
#' @param popA the first population to be compared with (usually a number)
#' @param popB the second population to be compared with
#' @param num_genes number of genes
getDEgenes <- function(true_counts_res, popA, popB, num_genes) {
  meta_cell <- true_counts_res$cell_meta
  meta_gene <- true_counts_res$gene_effects
  popA_idx <- which(meta_cell$pop == popA)
  popB_idx <- which(meta_cell$pop == popB)
  ngenes <- num_genes

  DEstr <- sapply(strsplit(colnames(meta_cell)[which(grepl("evf", colnames(meta_cell)))], "_"), "[[", 2)
  param_str <- sapply(strsplit(colnames(meta_cell)[which(grepl("evf", colnames(meta_cell)))], "_"), "[[", 1)
  n_useDEevf <- sapply(1:ngenes, function(igene) {
    return(sum(abs(meta_gene[[1]][igene, DEstr[which(param_str == "kon")] == "DE"]) - 0.001 > 0) +
             sum(abs(meta_gene[[2]][igene, DEstr[which(param_str == "koff")] == "DE"]) - 0.001 > 0) +
             sum(abs(meta_gene[[3]][igene, DEstr[which(param_str == "s")] == "DE"]) - 0.001 > 0))
  })

  kon_mat <- true_counts_res$kinetic_params[[1]]
  koff_mat <- true_counts_res$kinetic_params[[2]]
  s_mat <- true_counts_res$kinetic_params[[3]]

  logFC_theoretical <- sapply(1:ngenes, function(igene)
    return(log2(mean(s_mat[igene, popA_idx] * kon_mat[igene, popA_idx] / (kon_mat[igene, popA_idx] + koff_mat[igene, popA_idx])) /
                  mean(s_mat[igene, popB_idx] * kon_mat[igene, popB_idx] / (kon_mat[igene, popB_idx] + koff_mat[igene, popB_idx])))))

  true_counts <- true_counts_res$counts
  true_counts_norm <- t(t(true_counts) / colSums(true_counts)) * 10^6

  wil.p_true_counts <- sapply(1:ngenes, function(igene)
    return(wilcox.test(true_counts_norm[igene, popA_idx], true_counts_norm[igene, popB_idx])$p.value))

  wil.adjp_true_counts <- p.adjust(wil.p_true_counts, method = 'fdr')

  return(list(nDiffEVF = n_useDEevf, logFC_theoretical = logFC_theoretical, wil.p_true_counts = wil.p_true_counts))
}

#' retrieve the information of cells on continuous trajectories
#' Outputs a data frame, where each row corresponds to a cell
#' Each cell has information "pseudotime" (distance from root) and "branch" (on which branch is the cell)
#' @param cell_meta the cell meta information stored in the output of SimulateTrueCounts() or True2ObservedCounts()
getTrajectoryGenes <- function(cell_meta) {
  temp <- cell_meta[, 2:3]
  colnames(temp) <- c("branch", "pseudotime")
  rownames(temp) <- cell_meta[, 1]
  return(temp)
}

#' loads an RData file, and returns it
#' @param fileName RData file name
loadRData <- function(fileName){
    load(fileName)
    get(ls()[ls() != "fileName"])
}