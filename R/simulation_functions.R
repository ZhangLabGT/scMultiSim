#############################################################
# Master Equation Related Functions
#############################################################

#' @param num_cells number of cells
#' @param unregulated_to_regulated_gene_ratio ratio of unregulated genes to regulated genes.  The final number of genes is rounded up to the nearest 10
#' @param intrinsic_noise a 0 to 1 value which is the weight assigned to the random sample from the Beta-Poisson distribution, where the weight of the Beta-Poisson mean value is given a weight of 1 minus the intrinsic noise value.
#' @param Sigma parameter of the std of evf values within the same population
#' @param num_evfs number of evfs
#' @param diffEVF_fraction fraction of evfs which are diffEVFs
#' @param phyla the cell developmental tree if chosing 'discrete' or 'continuous' evf type. Can either be generated randomly (using pbtree(nclusters) function from phytools package) or read from newick format file using the ape package
#' @param beta splicing rate of each gene
#' @param d degradation rate of each gene
#' @param atac_effect a 0 to 1 value which is the weight assigned to the influence of chromatin accessability data on gene expression
#' @param num_cycles for generating velocity data, the number of cycles run before sampling the gene expression of a cell
#' @param cycle_length for generating velocity data, a factor which is multiplied by the expected time to transition from kon to koff and back to to form the the length of a cycle
#' @param randomseed should produce same result if all other parameters are all the same
#' @param do_velocity set to TRUE to incorporate the RNA velocity simulation module to simulate unspliced and spliced RNA sequence counts and ground truth velocity
#' @param metrics set to TRUE to calculate metrics assessing gene expression clustering, the influence of the GRN and chromatin accessibility data on gene expression, and RNA velocity
#' @return
#' @export
SimulateTrueCounts <- function(num_cells = 1000, unregulated_to_regulated_gene_ratio = 0.1, num_evfs = 500, diffEVF_fraction = 0.9, Sigma = 0.1,
                               atac_effect = 0.5, beta = 0.4, d = 1, num_cycles = 3, cycle_length = 1, intrinsic_noise = 1, randseed = 0, do_velocity = FALSE, metrics = FALSE, phyla = Phyla5()) {

  num_diffEVF <- num_evfs * diffEVF_fraction
  vary <- "s"; evf_type <- "continuous"
  evf_center <- 1; impulse <- F; geffect_mean <- 0; bimod <- 0; gene_effect_prob <- 0.3; gene_effects_sd <- 1
  param_realdata <- "zeisel.imputed"; scale_s <- 1; prop_hge <- 0.015; mean_hge <- 5; SE <- F
  set.seed(randseed)
  seed <- sample(1:1e5, size = 2)
  param_names <- c("kon", "koff", "s")
  load("scMultiSim/data/dens_nonzero.RData")  # this is the density function of log(x+1), where x is the non-zero values for ATAC-SEQ data

  regulator_ID_list <- sort(unique(target_gene_GRN_params[, 2]))
  target_gene_ID_list <- sort(unique(target_gene_GRN_params[, 1]))

  num_target_genes <- length(target_gene_ID_list)
  num_regulators <- length(regulator_ID_list)
  num_GRN_genes <- num_target_genes + num_regulators
  num_genes <- ceiling((num_GRN_genes + num_GRN_genes * unregulated_to_regulated_gene_ratio) / 10) * 10

  # get evf
  evf_all <- ContinuousEVF(phyla, num_cells, n_nd_evf = num_evfs - num_diffEVF, n_de_evf = num_diffEVF,
                           evf_center = evf_center, vary = vary, impulse = impulse,
                           Sigma = Sigma, plotting = F, seed = seed[1], num_regulators = num_regulators)

  gene_effects_by_regulator <- RegulatorGeneEffects(num_genes = num_genes, randseed = seed[2], target_gene_GRN_params = target_gene_GRN_params,
                                                    regulator_ID_list = regulator_ID_list, target_gene_ID_list = target_gene_ID_list)

  # gene_effects are the gene effect values, and have shape (kinetic parameter) x (num_genes) x (num_evfs + num_regulators)
  gene_effects <- GeneEffects(num_genes = num_genes, num_evfs = num_evfs, randseed = seed[2], prob = gene_effect_prob, geffect_mean = geffect_mean,
                              geffect_sd = gene_effects_sd, num_diffEVF = num_diffEVF,
                              regulator_ID_list = regulator_ID_list, target_gene_ID_list = target_gene_ID_list, gene_effects_by_regulator = gene_effects_by_regulator)

  # the scATAC-Seq data cannot be generated randomly. Using the SymSim framework,
  # we are generating the scATAC-Seq data using a similar idea as the scRNA-Seq data,
  # where we consider EVFs as cell IDs and region IDs similar to gene effects.

  nregions <- 3 * num_genes # regions 1 to nregion are considered sequentially located on the genome
  cumulative_probability <- cumsum(c(0.1, 0.5, 0.4)) # the probability that a gene is regulated by respectively 0, 1, 2 consecutive regions

  # generate region2gene, a 0 1 coupling matrix between regions and genes of shape (nregions) x (num_genes), where a value
  # of 1 indicates the gene is affected by a particular region
  region2gene <- Region2Gene(num_genes, nregions, cumulative_probability, randseed)

  # region_effects is like the gene_effect for the ATAC-Seq feature
  region_effect <- RegionEffects(nregions = nregions, nevf = num_evfs, randseed = seed[2], prob = gene_effect_prob,
                                 effect_mean = geffect_mean, effect_sd = gene_effects_sd)

  # to generate the EVF for ATAC-SEQ, we use the same EVF generation process as for s but with a different random seed, as s had DiffEVFs
  evf_ATAC <- ContinuousEVF(phyla, num_cells, n_nd_evf = num_evfs - num_diffEVF, n_de_evf = num_diffEVF,
                            evf_center = evf_center, vary = vary, impulse = impulse,
                            Sigma = Sigma, plotting = F, seed = seed[1] + 1, num_regulators = num_regulators)[[1]][[3]]

  # first generate scATAC-Seq data, only after this we can generate scRNA-Seq data
  atacseq_data <- ATAC_SEQ(evf_ATAC, region_effect, dens_nonzero, randseed)

  # get kon and koff & scale
  data(param_realdata.zeisel.imputed)
  match_params[, 1:3] <- log(base = 10, match_params[, 1:3])
  match_params_den <- lapply(1:3, function(i) {
    density(match_params[, i], n = 2000)
  })

  params <- Get_params(gene_effects, evf_all[[1]][c(1, 2)], match_params_den[c(1, 2)], bimod, scale_s = scale_s,
                       atacseq_data = t(atacseq_data), region2gene = region2gene, atac_effect = atac_effect)

  result <- RNA_SEQ(evf_all, phyla, num_genes, num_cells, num_regulators, evf_center, gene_effects, gene_effects_by_regulator, match_params_den,
                    params, intrinsic_noise, regulator_ID_list, scale_s, Sigma, do_velocity, beta, d, num_cycles, cycle_length)

  counts <- result$counts_s
  unspliced_counts <- result$counts_u
  cell_time <- result$cell_time
  velocity <- result$velocity

  cell_meta <- cbind(cellid=paste('cell',seq(1,num_cells),sep='_'),evf_all[[2]],evf_all[[1]])

  #################### code for assessment metrics ####################
  if (metrics) {

    if (do_velocity) {
      plotVelo(filename = "velocity.pdf", result, meta = evf_all[[2]], beta = beta, d = d, num_cycles = num_cycles, Sigma = Sigma, num_cells = num_cells, diffEVF_fraction = diffEVF_fraction, randseed = randseed, num_evfs = num_evfs, arrow.length = 1, phyla = phyla, neutral = evf_all[[3]][1:num_cells,])
    }

    # tsne true spliced counts
    PlotTsne(meta = evf_all[[2]], data = log2(counts + 1), label = 'pop', saving = T, plotname = 'Tsne True Spliced Counts')
    # tsne true unspliced counts
    PlotTsne(meta = evf_all[[2]], data = log2(result$counts_u + 1), label = 'pop', saving = T, plotname = 'Tsne True Unspliced Counts')
    # tsne atacseq
    PlotTsne(meta = evf_all[[2]], data = log2(atacseq_data + 1), label = 'pop', saving = T, plotname = 'Tsne True ATAC-seq Data')

    count_correlation_matrix <- cor(t(counts), method = "spearman")
    if (any(is.na(count_correlation_matrix))){
      print('some genes have no spliced counts across all cells; the correlation with these genes is set to 0 in the heatmap and is ignored in the mean within-gene-module correlation calculation')
      count_correlation_matrix = replace(count_correlation_matrix, which(is.na(count_correlation_matrix)), 0)
    }

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

    # jpeg("results_heatmap.jpeg", 600, 600)
    pdf("results_heatmap.pdf", 5, 5)
    aa <- heatmap.2(count_correlation_matrix, scale = "none", Rowv = T, Colv = T, dendrogram = "both", distfun = dist, hclustfun = hclust, key = T, trace = "none", cexRow = 1, cexCol = 1, RowSideColors = gene_module_color_vector, ColSideColors = gene_module_color_vector, col = bluered(75), main = '       Gene Corr by Main Regulator')
    dev.off()

    # jpeg("results_GRNheatmap.jpeg", 600, 600)
    pdf("results_GRNheatmap.pdf", 5, 5)
    bb <- heatmap.2(count_correlation_matrix[1:num_GRN_genes, 1:num_GRN_genes], scale = "none", Rowv = T, Colv = T, dendrogram = "both", distfun = dist, hclustfun = hclust, key = T, trace = "none", cexRow = 1, cexCol = 1, RowSideColors = gene_module_color_vector[1:num_GRN_genes], ColSideColors = gene_module_color_vector[1:num_GRN_genes], col = bluered(75), main = '          GRN Gene Corr by Main Regulator')
    dev.off()

    # output to outfile.txt the gene module correlation, RNA-ATAC correlation for genes effected by only 1 region, and the RNA-ATAC correlation for all genes
    outfile_header <- c("g_mod_cor", "ATAC1r_cor", "ATAC_cor")
    correlation_out_title = 'correlation results.txt'

    cat(outfile_header, '\n', file = correlation_out_title, sep = '\t', append = FALSE)

    cor_gene_module_pairs <- numeric()
    for (module_gene_index1 in seq_along(gene_module_color_vector[1:num_GRN_genes])) {
      for (module_gene_index2 in seq_along(gene_module_color_vector[1:num_GRN_genes])) {
        if ((module_gene_index1 != module_gene_index2) && (gene_module_color_vector[module_gene_index1] == gene_module_color_vector[module_gene_index2])) {
          cor_gene_module_pairs_next = suppressWarnings(cor(counts[module_gene_index1,], counts[module_gene_index2,], method = "spearman"))
          if (!is.na(cor_gene_module_pairs_next)){
          cor_gene_module_pairs <- c(cor_gene_module_pairs, cor_gene_module_pairs_next)
          }
        }
      }
    }
    gene_module_correlation <- mean(cor_gene_module_pairs)

    target_genes <- which(colSums(region2gene > 0) == 1)
    ATAC_1region_correlation <- numeric()
    for (gene_index in target_genes) {
      region <- which(region2gene[, gene_index] > 0)
      correlation <- suppressWarnings(cor(atacseq_data[region,], counts[gene_index,], method = "spearman"))
      ATAC_1region_correlation <- c(ATAC_1region_correlation, correlation)
    }
    ATAC_1region_correlation <- mean(ATAC_1region_correlation, na.rm = T)

    ATAC_correlation <- numeric()
    for (gene_index in 1:num_genes) {
      ATAC_correlation <- c(ATAC_correlation, suppressWarnings(cor(atacseq_data[gene_index,], counts[gene_index,], method = "spearman")))
    }
    ATAC_correlation <- mean(ATAC_correlation, na.rm = T)

    outfile_line <- c(gene_module_correlation, ATAC_1region_correlation, ATAC_correlation)
    cat(outfile_line, '\n', file = correlation_out_title, sep = '\t', append = TRUE)
  }


  return(list(counts=counts, cell_meta=cell_meta, unspliced_counts=unspliced_counts, velocity=velocity, cell_time=cell_time, atacseq_data=atacseq_data, num_genes=num_genes))
}

#' Generate true transcript counts for linear structure
#' @param kinet_params kinetic parameters, include k_on, k_off, s and beta
#' @param start_state the starting state: on or off of each gene
#' @param start_s spliced count of the root cell in the branch
#' @param start_u unspliced count of the root cell in the branch
#' @param randpoints1 the value which evf mean is generated from
#' @param ncells1 number of cells in the branch
#' @param ngenes number of genes
#' @param beta_vec splicing rate of each gene
#' @param d_vec degradation rate of each gene
#' @return a list of 4 elements, the first element is true counts, second is the gene level meta information, the third is cell level meta information, including a matrix of evf and a vector of cell identity, and the fourth is the parameters kon, koff and s used to simulation the true counts
#' @import phytools
#' @export
gen_1branch <- function(kinet_params, start_state, start_s, start_u, randpoints1, ncells1, ngenes, beta_vec, d_vec, cycle_length_factor, cell) {

  # totaltime equals to total numeber of cells
  totaltime <- ncells1

  # store unspliced count
  counts_u1 <- matrix(0, ngenes, ncells1)
  # store spliced count
  counts_s1 <- matrix(0, ngenes, ncells1)
  # store state matrix
  state_mat <- matrix(0, ngenes, ncells1)

  # store the kinetic values, include k_on, k_off, s and beta

  kinet_params <- list(k_on = kinet_params$k_on, k_off = kinet_params$k_off, s = kinet_params$s,
                        beta = beta_vec, d = d_vec)

  xs <- list()
  ys <- list()
  which_cells <- list()

  for (igene in 1:ngenes) {
    k_on <- kinet_params$k_on[igene]
    k_off <- kinet_params$k_off[igene]
    s <- kinet_params$s[igene]
    beta <- kinet_params$beta[igene]
    d <- kinet_params$d[igene]

    cycle_length <- 1 / k_on + 1 / k_off
    min_wtime <- min(1 / k_on, 1 / k_off)
    npart <- max(ceiling(cycle_length / min_wtime) * cycle_length_factor, cycle_length * cycle_length_factor)
    stepsize <- cycle_length / npart

    nsteps <- ceiling(totaltime / stepsize)
    if (is.null(start_s)){
      x <- numeric(nsteps); x[1] <- s * k_on / (k_on + k_off) / beta
      y <- numeric(nsteps); y[1] <- s * k_on / (k_on + k_off) / d
    } else {
      x <- numeric(nsteps)
      y <- numeric(nsteps)
      x[1] <- start_u[igene]
      y[1] <- start_s[igene]
    }

    # which_cell stores with length nsteps stores the cell correspond to current time step.
    which_cell <- numeric(nsteps); which_cell[1] <- 1
    curr_time <- numeric(nsteps);
    p_table <- matrix(0, 2, 2)
    p_table[1, 2] <- 1 / (npart * (1 / k_off / cycle_length)); p_table[1, 1] <- 1 - p_table[1, 2];
    p_table[2, 1] <- 1 / (npart * (1 / k_on / cycle_length)); p_table[2, 2] <- 1 - p_table[2, 1]

    curr_state <- numeric(nsteps); curr_state[1] <- start_state[igene] # 1 means on state, 2 means off state
    t <- 1

    while (TRUE) {
      t <- t + 1
      curr_time[t] <- curr_time[t - 1] + stepsize
      if ((curr_time[t] - which_cell[t - 1]) > 0) {
        which_cell[t] <- which_cell[t - 1] + 1
        if (which_cell[t] > ncells1) {
          break
        }

        # check npart code here, on p_table, update the step size
        cycle_length <- 1 / k_on + 1 / k_off
        min_wtime <- min(1 / k_on, 1 / k_off)
        npart <- max(ceiling(cycle_length / min_wtime) * cycle_length_factor, cycle_length * cycle_length_factor)
        stepsize <- cycle_length / npart
        p_table <- matrix(0, 2, 2)
        p_table[1, 2] <- 1 / (npart * (1 / k_off / cycle_length)); p_table[1, 1] <- 1 - p_table[1, 2];
        p_table[2, 1] <- 1 / (npart * (1 / k_on / cycle_length)); p_table[2, 2] <- 1 - p_table[2, 1]
      } else { which_cell[t] <- which_cell[t - 1] }

      if (runif(1, 0, 1) > p_table[curr_state[t - 1], 1]) {
        curr_state[t] <- 2
      } else { curr_state[t] <- 1 }

      if (curr_state[t] == 1) {
        x[t] <- x[t - 1] + s * stepsize - beta * x[t - 1] * stepsize
      } else {
        x[t] <- x[t - 1] - beta * x[t - 1] * stepsize
      }
      if (x[t] < 0) { x[t] <- 0 }
      y[t] <- y[t - 1] + beta * x[t - 1] * stepsize - d * y[t - 1] * stepsize
      if (y[t] < 0) { y[t] <- 0 }

    }
    xs[[igene]] <- x; ys[[igene]] <- y; which_cells[[igene]] <- which_cell
    # extract value for each cell
    for (icell in 1:ncells1) {
      all_idx <- which(which_cell == icell)
      closest <- which.min(abs((curr_time[all_idx] - (icell - 1)) - randpoints1[icell]))
      counts_u1[igene, icell] <- x[all_idx[closest]]
      counts_s1[igene, icell] <- y[all_idx[closest]]
      state_mat[igene, icell] <- curr_state[all_idx[closest]]
    }
  }

  cell_time <- randpoints1 + (0:(ncells1 - 1))

  # calculate the ground truth velocity
  velo_mat <- beta_vec * counts_u1 - d_vec * counts_s1

  dynamics <- list()
  dynamics[["unspliced_steps"]] <- xs
  dynamics[["spliced_steps"]] <- ys
  dynamics[["which_cells"]] <- which_cells

  return(list(counts_u = counts_u1, counts_s = counts_s1, kinet_params = kinet_params, state_mat = state_mat, cell_time = cell_time, velocity = velo_mat, dynamics = dynamics))
}

#' Getting RNA_SEQ true counts by cell one by one taking into account the GRN conditions of the previous cell
#'
#' @param num_cells number of cells
#' @param num_genes number of genes
#' @param num_regulators number of regulator genes
#' @param regulator_ID_list list of the IDs of regulator genes
#' @param gene_effects a list of three matrices (generated using the GeneEffects function),
#' @param gene_effects_by_regulator a list of three matrices (generated using the GeneEffects function),
#' @param params the parameters for simulating gene expression from EVf and gene effects
#' @param intrinsic_noise a 0 to 1 value which is the weight assigned to the random sample from the Beta-Poisson distribution, where the weight of the Beta-Poisson mean value is given a weight of 1 minus the intrinsic noise value.
#' @param Sigma parameter of the std of evf values within the same population
#' @param phyla the cell developmental tree if chosing 'discrete' or 'continuous' evf type. Can either be generated randomly (using pbtree(nclusters) function from phytools package) or read from newick format file using the ape package
#' @param scale_s a factor to scale the s parameter, which is used to tune the size of the actual cell (small cells have less number of transcripts in total)
#' @param evf_all a list of two object, one is the evf, and the other is a dataframe indicating the branch each cell comes from (pop) and its depth in the tree
#' @return generate a true RNA SEQ count matrix, which has the shape (genes) x (cells)

RNA_SEQ <- function(evf_all, phyla, num_genes, num_cells, num_regulators, evf_center, gene_effects, gene_effects_by_regulator, match_params_den,
                    params, intrinsic_noise, regulator_ID_list, scale_s, Sigma, do_velocity, beta, d, num_cycles, cycle_length) {
  # set.seed(randseed)
  # cells one by one
  neutral <- evf_all[[3]][1:num_cells,]
  edges <- phyla$edge
  root <- setdiff(edges[, 1], edges[, 2])
  evf_reg <- matrix(0, nrow = num_regulators, ncol = num_cells)
  evf_start <- rnorm(num_regulators, evf_center, Sigma)
  cellcount <- 1

  num_cycles = num_cycles
  d_vec <- rnorm(n = num_genes, mean = d, sd = 0.1) # degradation rate of the size ngenes * 1
  beta_vec <- rnorm(n = num_genes, mean = beta, sd = 0.1) # splicing rate of the size ngenes * 1
  counts_u <- matrix(nrow = num_genes, ncol = num_cells)
  counts_s <- matrix(nrow = num_genes, ncol = num_cells)
  state_mat <- matrix(nrow = num_genes, ncol = num_cells)
  cell_time <- vector(length = num_cells)
  velocity <- matrix(nrow = num_genes, ncol = num_cells)
  # generate the promoter states of the root cell with dimension ngenes * 1
  root_state <- sample(c(1, 2), size = num_genes, replace = TRUE)

  # for each edge in the tree
  for (ii in seq(1:dim(edges)[1])) {

    # get current parent and child
    parent <- edges[ii, 1]
    child <- edges[ii, 2]
    # find the set of cell index numbers for cells which belong to the edge (have the parent and child associated with the edge)
    edge_cell_idxs <- which(neutral[, 1] == parent & neutral[, 2] == child, arr.ind = T)
    # num of cells on this edge
    ncells_edge <- length(edge_cell_idxs)

    # first evf on this edge # curr_evf is the one for previous cell, donot store it
    if (parent == root) {
      curr_evf <- evf_start
    }else {
      grandparent <- edges[edges[, 2] == parent, 1]
      parent_last_cell_idx <- max(which(neutral[, 1] == grandparent & neutral[, 2] == parent, arr.ind = T))
      curr_evf <- evf_reg[, parent_last_cell_idx]
    }
    if (all(curr_evf == 0)) { print('no evf') }

    counts_edge <- matrix(0, nrow = num_genes, ncol = ncells_edge)
    evf_edge <- matrix(0, nrow = num_regulators, ncol = ncells_edge)

    # for each cell on the edge
    for (jj in seq(1:ncells_edge)) {
      # get s for a cell
      s_cell <- evf_all[[1]][[3]][cellcount,] %*% t(gene_effects[[3]]) + curr_evf %*% t(gene_effects_by_regulator)
      # scale
      temp <- alply(s_cell, 1, function(Y) { Y })
      values <- do.call(c, temp)
      ranks <- rank(values)
      sorted <- sort(SampleDen(nsample = max(ranks), den_fun = match_params_den[[3]]))
      s_cell <- matrix(data = sorted[ranks], ncol = length(s_cell))
      s_cell <- t(apply(s_cell, 2, function(x) { x <- 10^x })) * scale_s
      s_cell <- t(s_cell)

      if (do_velocity == FALSE) {
        # get counts from beta poisson model
        counts_cell <- sapply(1:num_genes, function(i) {
          y <- rbeta(1, params[[1]][i, cellcount], params[[2]][i, cellcount])
          x <- rpois(1, y * s_cell[i])
          y_mean <- params[[1]][i, cellcount] / (params[[1]][i, cellcount] + params[[2]][i, cellcount])
          x_mean <- y_mean * s_cell[i]
          return(intrinsic_noise * x + (1 - intrinsic_noise) * x_mean)
        })
      } else {
        # get spliced/unspliced counts and velocity from kinetic model
        s_cell <- as.vector(s_cell)
        cell_idx = edge_cell_idxs[[jj]]

        kinetic_params_all <- list()
        kinetic_params_all$k_on <- params[[1]][, cellcount]
        kinetic_params_all$k_off <- params[[2]][, cellcount]
        kinetic_params_all$s <- s_cell

        if (jj == 1) {
          if (parent == root) {
            start_s <- NULL
            start_u <- NULL
            start_state <- root_state
            start_cell_time <- 0
            cycles = 15
          } else {
            start_s <- counts_s[, parent_last_cell_idx]
            start_u <- counts_u[, parent_last_cell_idx]
            start_state <- state_mat[, parent_last_cell_idx]
            start_cell_time <- cell_time[parent_last_cell_idx]
            cycles = num_cycles
          }
        } else {
          last_cell_idx = edge_cell_idxs[[jj - 1]]
          start_s <- counts_s[, last_cell_idx]
          start_u <- counts_u[, last_cell_idx]
          start_state <- state_mat[, last_cell_idx]
          start_cell_time <- cell_time[last_cell_idx]
          cycles = num_cycles
        }
        randpoints <- runif(n = cycles, min = 0, max = 1)

        result <- gen_1branch(kinet_params = kinetic_params_all, start_state = start_state, start_u = start_u, start_s = start_s, cycle_length_factor = cycle_length,
                              randpoints1 = randpoints, ncells1 = cycles, ngenes = num_genes, beta_vec = beta_vec, d_vec = d_vec, cell = cellcount)
        counts_u[, cell_idx] <- result$counts_u[, cycles]
        counts_s[, cell_idx] <- result$counts_s[, cycles]
        state_mat[, cell_idx] <- result$state_mat[, cycles]
        cell_time[cell_idx] <- result$cell_time[[cycles]] + start_cell_time
        velocity[, cell_idx] <- result$velocity[, cycles]

        counts_cell <- counts_s[, cell_idx]
      }
      # counts of regulators
      counts_reg <- counts_cell[regulator_ID_list]
      evf_thiscell <- counts_reg / (counts_reg + mean(counts_cell))
      counts_edge[, jj] <- counts_cell
      evf_edge[, jj] <- evf_thiscell
      curr_evf <- evf_thiscell
      cellcount <- cellcount + 1
    }
    counts_s[, edge_cell_idxs] <- counts_edge
    evf_reg[, edge_cell_idxs] <- evf_edge
  }
  return(list(counts_u = counts_u, counts_s = counts_s, cell_time = cell_time, velocity = velocity))

  #####################################################################################################################################
}

#' Getting ATAC SEQ data
#'
#' @param evf_ATAC EVF values for the ATAC features
#' @param region_effect like the gene_effect for the ATAC-Seq feature
#' @param randomseed should produce same result if all other parameters are all the same
#' @return ATAC SEQ data
ATAC_SEQ <- function(evf_ATAC, region_effect, dens_nonzero, randseed) {
  set.seed(randseed)
  params_region <- evf_ATAC %*% t(region_effect)

  # now rescale these values to a distribution of real scATAC-Seq dataset
  p_zero <- 0.8 # the proportion of 0s we see in the data

  temp <- alply(params_region, 1, function(Y) { Y })
  values <- do.call(c, temp)
  nzeros <- floor(length(values) * p_zero); nnonzero <- length(values) - nzeros
  ranks <- rank(values)
  mapped_data <- numeric(length(values))
  sorted_nonzeros <- sort(SampleDen(nsample = nnonzero, den_fun = dens_nonzero))
  mapped_data[which(ranks > nzeros)] <- 2^(sorted_nonzeros[ranks[ranks > nzeros] - nzeros]) - 1
  atacseq_data <- matrix(data = mapped_data, ncol = length(params_region[1,]), byrow = T)
  atacseq_data <- t(atacseq_data)
  return(atacseq_data)
}

#' Getting Region2Gene matrices
#'
# generate region2gene, a 0 1 coupling matrix between regions and genes of shape (nregions) x (num_genes), where a value
# of 1 indicates the gene is affected by a particular region
#' @param num_genes number of genes
#' @param nregions # regions 1 to nregion are considered sequentially located on the genome
#' @param cumulative_probability # the probability that a gene is regulated by respectively 0, 1, 2 consecutive regions
#' @param randomseed should produce same result if all other parameters are all the same
#' @return region2gene, of shape (nregions) x (num_genes)
Region2Gene <- function(num_genes, nregions, cumulative_probability, randseed) {
  set.seed(randseed)
  region2gene <- matrix(0, nregions, num_genes)
  random_uniforms <- runif(num_genes) # generate a vector of random uniform values of length num_genes to pick the number of regions effecting each gene

  for (gene_index in 1:num_genes) {
    # if a gene is regulated by 1 region, pick a random region and set the region2gene index corresponding to the region-gene pair to 1
    if (random_uniforms[gene_index] >= cumulative_probability[1] & random_uniforms[gene_index] < cumulative_probability[2]) {
      region2gene[ceiling(runif(1, min = 0, max = nregions)), gene_index] <- 1
      # if a gene is regulated by 2 regions, pick a random region and set the region2gene index corresponding to the region-gene pair to 1 for the selected region and the next region
    } else if (random_uniforms[gene_index] >= cumulative_probability[2]) {
      startpos <- ceiling(runif(1, min = 0, max = nregions - 1))
      region2gene[startpos:(startpos + 1), gene_index] <- c(1, 1)
    }
  }
  return(region2gene)
}

#' Getting RegulatorGeneEffects matrices
#'
#' @param num_genes number of genes
#' @param target_gene_GRN_params
#' @param randomseed should produce same result if all other parameters are all the same
#' @param regulator_ID_list list of the IDs of regulator genes
#' @param target_gene_ID_list list of the IDs of target genes
#' @return a list of 3 matrices, each of dimension ngenes * nevf
RegulatorGeneEffects <- function(num_genes, randseed, target_gene_GRN_params, regulator_ID_list, target_gene_ID_list) {
  set.seed(randseed)
  #################### calculate gene effect values of the regulators on the target (regulated) genes ####################
  num_target_genes <- length(target_gene_ID_list)
  num_regulators <- length(regulator_ID_list)
  num_GRN_genes <- num_target_genes + num_regulators

  gene_effects_by_regulator <- matrix(0L, nrow = num_GRN_genes, ncol = num_regulators)
  for (row_index in seq_along(target_gene_GRN_params[, 2])) {
    target_gene <- target_gene_GRN_params[row_index, 1]
    regulator_gene <- target_gene_GRN_params[row_index, 2]
    regulator_gene_effect <- target_gene_GRN_params[row_index, 3]
    gene_effects_by_regulator[target_gene, which(regulator_ID_list %in% regulator_gene)] <- regulator_gene_effect
  }

  if (num_genes > nrow(gene_effects_by_regulator)) {
    gene_effects_by_regulator <- rbind(gene_effects_by_regulator, matrix(0, nrow = num_genes - nrow(gene_effects_by_regulator), ncol = num_regulators))
  }
  return(gene_effects_by_regulator)
}

#' Getting GeneEffects matrices 
#'
#' This function randomly generates the effect size of each evf on the dynamic expression parameters
#' @param num_genes number of genes
#' @param num_evfs number of evfs
#' @param randomseed (should produce same result if ngenes, nevf and randseed are all the same)
#' @param prob the probability that the effect size is not 0
#' @param geffect_mean the mean of the normal distribution where the non-zero effect sizes are dropped from 
#' @param geffect_sd the standard deviation of the normal distribution where the non-zero effect sizes are dropped from 
#' @return a list of 3 matrices, each of dimension ngenes * nevf
GeneEffects <- function(num_genes, num_evfs, randseed, prob, geffect_mean, geffect_sd,
                        num_diffEVF, regulator_ID_list, target_gene_ID_list, gene_effects_by_regulator) {
  set.seed(randseed)
  num_target_genes <- length(target_gene_ID_list)
  num_regulators <- length(regulator_ID_list)
  num_GRN_genes <- num_target_genes + num_regulators
  #################### calculate initial gene effect values with original Symsim approach ####################
  gene_effects <- lapply(c('kon', 'koff', 's'), function(param) {
    effect <- lapply(c(1:num_genes), function(i) {
      nonzero <- sample(size = num_evfs, x = c(0, 1), prob = c((1 - prob), prob), replace = T)
      nonzero[nonzero != 0] = rnorm(sum(nonzero), mean = geffect_mean, sd = geffect_sd)
      return(nonzero)
    })
    return(do.call(rbind, effect))
  })

  #################### add the effect of regulator gene effects to gene_effects ####################
  # the regulators should use only DiffEVFs:
  regulator_diff_EVF <- matrix(0, num_regulators, num_diffEVF)
  for (i in seq_along(regulator_ID_list)) {
    chosen <- sample(1:num_diffEVF, replace = F)
    regulator_diff_EVF[i, chosen[1:2]] <- 2 # gene effect of 2 added to two random differential evfs for each master regulator gene
  }
  gene_effects[[3]][(regulator_ID_list),] <- cbind(matrix(0, num_regulators, num_evfs - num_diffEVF), regulator_diff_EVF)

  # for every target gene, it should use the same gene_effect vector as its regulators
  # if a gene has multiple regulators, its gene effects will be the combination of that of the regulators.
  # we do it for the gene-effect for s only.
  for (target_gene_ID in (target_gene_ID_list)) {
    gene_effects[[3]][target_gene_ID,] <- (gene_effects_by_regulator[target_gene_ID,] %*% gene_effects[[3]][regulator_ID_list,]) /
      sum(gene_effects_by_regulator[target_gene_ID,] > 0) /
      2
    gene_effects[[3]][target_gene_ID, which(abs(gene_effects[[3]][target_gene_ID,]) < 0.2)] <- 0
  }

  # increase the scale of gene-effects for the new genes which do not have regulators.
  if (num_genes > num_GRN_genes) {
    gene_effects[[3]][(num_GRN_genes + 1):num_genes,] <- gene_effects[[3]][(num_GRN_genes + 1):num_genes,] * 3
  }

  # add gene effect value columns to the k_on and k_off kinetic parameter matrices.  The number of columns added is equal to the number of regulators,
  # so that k_on and k_off have the same shape as the rate of transcription s, since s had columns added equal to the number of regulators representing
  # the effect of the regulators on the transcription rate s.  Since the regulators do not affect k_on and k_off, the added columns are zero-valued
  gene_effects[[1]] <- cbind(gene_effects[[1]], matrix(0, ncol = num_regulators, nrow = num_genes))
  gene_effects[[2]] <- cbind(gene_effects[[2]], matrix(0, ncol = num_regulators, nrow = num_genes))

  return(gene_effects)
}

#' Getting RegionEffect matrix
#'
#' @param nregions number of regions
#' @param nevf number of evfs
#' @param randomseed (should produce same result if nregions, nevf and randseed are all the same)
#' @param prob the probability that the effect size is not 0
#' @param effect_mean the mean of the normal distribution where the non-zero effect sizes are dropped from
#' @param effect_sd the standard deviation of the normal distribution where the non-zero effect sizes are dropped from
#' @return a matrix of dimension nregions * nevf
RegionEffects <- function(nregions, nevf, randseed, prob, effect_mean, effect_sd) {
  set.seed(randseed)
  region_effects <- lapply(c(1:nregions), function(i) { # region_effects is like the gene_effect for the ATAC-Seq feature
    nonzero <- sample(size = nevf, x = c(0, 1), prob = c((1 - prob), prob), replace = T)
    nonzero[nonzero != 0] = rnorm(sum(nonzero), mean = effect_mean, sd = effect_sd)
    return(nonzero)
  })
  region_effects <- do.call(rbind, region_effects)
  return(region_effects)
}

#' sample from smoothed density function
#' @param nsample number of samples needed
#' @param den_fun density function estimated from density() from R default
SampleDen <- function(nsample, den_fun) {
  probs <- den_fun$y / sum(den_fun$y)
  bw <- den_fun$x[2] - den_fun$x[1]
  bin_id <- sample(size = nsample, x = c(1:length(probs)), prob = probs, replace = T)
  counts <- table(bin_id)
  sampled_bins <- as.numeric(names(counts))
  samples <- lapply(c(1:length(counts)), function(j) {
    runif(n = counts[j], min = (den_fun$x[sampled_bins[j]] - 0.5 * bw), max = (den_fun$x[sampled_bins[j]] + 0.5 * bw))
  })
  samples <- do.call(c, samples)
  return(samples)
}

#' Getting the parameters for simulating gene expression from EVf and gene effects
#'
#' This function takes gene_effect and EVF, take their dot product and scale the product to the correct range
#' by using first a logistic function and then adding/dividing by constant to their correct range
#' @param gene_effects a list of three matrices (generated using the GeneEffects function),
#' each corresponding to one kinetic parameter. Each matrix has nevf columns, and ngenes rows.
#' @param evf a vector of length nevf, the cell specific extrinsic variation factor
#' @param match_param_den the fitted parameter distribution density to sample from
#' @param bimod the bimodality constant
#' @param scale_s a factor to scale the s parameter, which is used to tune the size of the actual cell (small cells have less number of transcripts in total)
#' @return params a matrix of ngenes * 3
#' @examples
#' Get_params()

# params <- Get_params(gene_effects,evf_res[[1]],match_params_den,bimod,scale_s=scale_s,
# atacseq_data, region2gene)
Get_params <- function(gene_effects, evf, match_param_den, bimod, scale_s, atacseq_data, region2gene, atac_effect) {
  # For koff and s, treat all genes in the same way
  # For kon, consider the effect of ATAC-Seq
  params <- lapply(1:length(evf), function(iparam) { evf[[iparam]] %*% t(gene_effects[[iparam]]) })
  scaled_params <- lapply(c(1:length(evf)), function(i) {
    X <- params[[i]]
    temp <- alply(X, 1, function(Y) { Y })
    values <- do.call(c, temp)
    ranks <- rank(values)
    sorted <- sort(SampleDen(nsample = max(ranks), den_fun = match_param_den[[i]]))
    temp3 <- matrix(data = sorted[ranks], ncol = length(X[1,]), byrow = T)
    return(temp3)
  })

  # calculate kon. Re-calculate kon for the genes which are regulated by regions
  # first, normalize for each column in region2gene
  colidx <- which(colSums(region2gene) > 1)
  for (icol in colidx) {
    region2gene[, icol] <- region2gene[, icol] / sum(region2gene[, icol])
  }
  kon_atac <- atacseq_data %*% region2gene

  mask_matrix <- matrix(0, dim(kon_atac)[1], dim(kon_atac)[2])
  colidx <- which(colSums(region2gene) > 0)
  mask_matrix[, colidx] <- matrix(1, dim(kon_atac)[1], length(colidx))
  idx_consider <- which(mask_matrix > 0)
  # replace the 0 values by values from the original kon values
  zero_idx <- which(kon_atac == 0 & mask_matrix == 1)
  temp <- scaled_params[[1]][zero_idx]
  temp <- temp - min(temp)
  temp <- temp * (min(kon_atac[kon_atac > 0]) / 2 / max(temp))
  kon_atac[zero_idx] <- temp
  consider_values <- kon_atac[idx_consider]
  ranks <- atac_effect * rank(consider_values) + (1 - atac_effect) * rank(scaled_params[[1]][idx_consider])
  sorted <- sort(SampleDen(nsample = max(ranks), den_fun = match_param_den[[1]]))
  scaled_params[[1]][idx_consider] <- sorted[ranks]

  # adjust parameters with the bimod parameter
  bimod_perc <- 1
  ngenes <- dim(scaled_params[[1]])[2]; bimod_vec <- numeric(ngenes)
  bimod_vec[1:ceiling(ngenes * bimod_perc)] <- bimod
  bimod_vec <- c(rep(bimod, ngenes / 2), rep(0, ngenes / 2))
  scaled_params[[1]] <- apply(t(scaled_params[[1]]), 2, function(x) { x <- 10^(x - bimod_vec) })
  scaled_params[[2]] <- apply(t(scaled_params[[2]]), 2, function(x) { x <- 10^(x - bimod_vec) })
  if (length(evf) == 3) {
    scaled_params[[3]] <- t(apply(scaled_params[[3]], 2, function(x) { x <- 10^x })) * scale_s
  }
  return(scaled_params)
}


#' Getting the parameters for simulating gene expression from EVf and gene effects
#'
#' This function takes gene_effect and EVF, take their dot product and scale the product to the correct range
#' by using first a logistic function and then adding/dividing by constant to their correct range
#' @param evf a vector of length nevf, the cell specific extrinsic variation factor
#' @param gene_effects a list of three matrices (generated using the GeneEffects function),
#' each corresponding to one kinetic parameter. Each matrix has nevf columns, and ngenes rows.
#' @param param_realdata the fitted parameter distribution to sample from
#' @param bimod the bimodality constant
#' @return params a matrix of ngenes * 3
#' @examples
#' Get_params()
Get_params2 <- function(gene_effects, evf, bimod, ranges) {
  params <- lapply(gene_effects, function(X) { evf %*% t(X) })
  scaled_params <- lapply(c(1:3), function(i) {
    X <- params[[i]]
    temp <- apply(X, 2, function(x) { 1 / (1 + exp(-x)) })
    temp2 <- temp * (ranges[[i]][2] - ranges[[i]][1]) + ranges[[i]][1]
    return(temp2)
  })
  scaled_params[[1]] <- apply(scaled_params[[1]], 2, function(x) { x <- 10^(x - bimod) })
  scaled_params[[2]] <- apply(scaled_params[[2]], 2, function(x) { x <- 10^(x - bimod) })
  scaled_params[[3]] <- apply(scaled_params[[3]], 2, function(x) { x <- abs(x) })
  scaled_params <- lapply(scaled_params, t)
  return(scaled_params)
}

#' @export
get_prob <- function(glength) {
  if (glength >= 1000) { prob <- 0.7 } else {
    if (glength >= 100 & glength < 1000) { prob <- 0.78 }
    else if (glength < 100) { prob <- 0 }
  }
  return(prob)
}

#' This function simulates the amplification, library prep, and the sequencing processes.
#' @param true_counts_1cell the true transcript counts for one cell (one vector)
#' @param protocol a string, can be "nonUMI" or "UMI"
#' @param rate_2cap the capture efficiency for this cell
#' @param gene_len gene lengths for the genes/transcripts, sampled from real human transcript length
#' @param amp_bias amplification bias for each gene, a vector of length ngenes
#' @param rate_2PCR PCR efficiency, usually very high
#' @param nPCR1 the number of PCR cycles
#' @param LinearAmp if linear amplification is used for pre-amplification step, default is FALSE
#' @param LinearAmp_coef the coeficient of linear amplification, that is, how many times each molecule is amplified by
#' @param N_molecules_SEQ number of molecules sent for sequencing; sequencing depth
#' @return read counts (if protocol="nonUMI") or UMI counts (if protocol="UMI)
#' @export
amplify_1cell <- function(true_counts_1cell, protocol, rate_2cap, gene_len, amp_bias,
                          rate_2PCR, nPCR1, nPCR2, LinearAmp, LinearAmp_coef, N_molecules_SEQ) {
  ngenes <- length(gene_len)
  if (protocol == "nonUMI") { data(len2nfrag) } else
    if (protocol == "UMI") { } else
    { stop("protocol input should be nonUMI or UMI") }
  inds <- vector("list", 2)
  print('true_counts_1cell')
  print(true_counts_1cell)
  # expand the original vector and apply capture efficiency
  # maintain a transcript index vector: which transcript the molecule belongs to
  expanded_res <- expand2binary(c(true_counts_1cell, 1))
  expanded_vec <- expanded_res[[1]]
  trans_idx <- expanded_res[[2]]
  print('expanded_vec')
  print(expanded_vec)
  inds[[1]] <- which(expanded_vec > 0)
  expanded_vec <- expanded_vec[inds[[1]]]
  print('inds[[1]]')
  print(inds[[1]])
  print('trans_idx 768')
  print(trans_idx)
  trans_idx <- trans_idx[inds[[1]]]
  print('trans_idx 771')
  print(trans_idx)

  captured_vec <- expanded_vec; captured_vec[runif(length(captured_vec)) > rate_2cap] <- 0
  if (sum(captured_vec[1:(length(captured_vec) - 1)]) < 1) { return(rep(0, ngenes)) }
  captured_vec[length(captured_vec)] <- 1

  inds[[2]] <- which(captured_vec > 0);
  print('captured_vec')
  print(captured_vec)
  captured_vec <- captured_vec[inds[[2]]]
  print('captured_vec')
  print(captured_vec)
  print('inds2')
  print(inds[[2]])
  print('trans_idx')
  print(trans_idx)
  trans_idx <- trans_idx[inds[[2]]]
  print('trans_idx')
  print(trans_idx)
  print('rate_2PCR')
  print(rate_2PCR)
  print('amp_bias')
  print(amp_bias)
  amp_rate <- c((rate_2PCR + amp_bias[trans_idx[1:(length(trans_idx) - 1)]]), 1)
  # pre-amplification:
  if (LinearAmp) {
    PCRed_vec <- captured_vec * LinearAmp_coef
  } else {
    temp <- runif(length(captured_vec)) < amp_rate
    print('org temp')
    print(temp)
    temp <- temp * 2 + captured_vec - temp
    for (iPCR in 2:nPCR1) {
      eff <- runif(length(temp)) * amp_rate
      print('amprate')
      print(amp_rate)
      print('temp')
      print(temp)
      v1 <- temp * (1 - eff)
      round_down <- (v1 - floor(v1)) < runif(length(v1))
      v1[round_down] <- floor(v1[round_down])
      v1[!round_down] <- ceiling(v1[!round_down])
      temp <- v1 + 2 * (temp - v1)
    }
    PCRed_vec <- temp
  }

  if (protocol == "nonUMI") { # add fragmentation step here
    temp_vec <- PCRed_vec
    for (i in seq(2, 1, -1)) {
      temp_vec1 <- numeric(); temp_vec1[inds[[i]]] <- temp_vec;
      temp_vec <- temp_vec1; temp_vec[is.na(temp_vec)] <- 0
    }
    recovered_vec <- temp_vec[1:(length(temp_vec) - 1)]
    amp_mol_count = numeric(ngenes);
    GI = c(0, cumsum(true_counts_1cell));
    for (i in which(true_counts_1cell > 0)) {
      x = recovered_vec[(GI[i] + 1):GI[i + 1]]
      amp_mol_count[i] = sum(x)
    }

    # for every copy of each transcript, convert it into number of fragments
    frag_vec <- numeric(ngenes)
    for (igene in which(amp_mol_count > 0)) {
      frag_vec[igene] <- sum(sample(len2nfrag[as.character(gene_len[igene]),],
                                    amp_mol_count[igene], replace = TRUE)) }
    # another 8 rounds of amplification to the fragments (fragmentation bias gets amplified)
    for (iPCR in 1:2) {
      frag_vec <- frag_vec + sapply(frag_vec, function(x) rbinom(n = 1, x, prob = rate_2PCR))
    }
    for (iPCR in 3:nPCR2) {
      frag_vec <- frag_vec + round(frag_vec * rate_2PCR)
    }
    SEQ_efficiency = N_molecules_SEQ / sum(frag_vec)
    if (SEQ_efficiency >= 1) { read_count <- frag_vec } else {
      read_count <- sapply(frag_vec, function(Y) { rbinom(n = 1, size = Y, prob = SEQ_efficiency) }) }
    return(read_count)
  } else if (protocol == "UMI") {

    prob_vec <- sapply(gene_len[trans_idx[1:(length(trans_idx) - 1)]], get_prob)
    # fragmentation:
    frag_vec <- sapply(1:(length(PCRed_vec) - 1), function(igene)
    { return(rbinom(n = 1, size = PCRed_vec[igene], prob = prob_vec[igene])) })

    # another 10 rounds of amplification to the fragments (fragmentation bias gets amplified)
    for (iPCR in 1:2) {
      frag_vec <- frag_vec + sapply(frag_vec, function(x) rbinom(n = 1, x, prob = rate_2PCR))
    }

    frag_vec <- round(frag_vec * (1 + rate_2PCR)^(nPCR2 - 1))

    SEQ_efficiency <- N_molecules_SEQ / sum(frag_vec)
    if (SEQ_efficiency >= 1) { sequenced_vec <- frag_vec } else {
      sequenced_vec <- sapply(frag_vec, function(Y) { rbinom(n = 1, size = Y, prob = SEQ_efficiency) }) }

    temp_vec <- c(sequenced_vec, 1)
    for (i in seq(2, 1, -1)) {
      temp_vec1 <- numeric(); temp_vec1[inds[[i]]] <- temp_vec;
      temp_vec <- temp_vec1; temp_vec[is.na(temp_vec)] <- 0
    }
    recovered_vec <- temp_vec[1:(length(temp_vec) - 1)]

    UMI_counts = numeric(ngenes);
    GI = c(0, cumsum(true_counts_1cell));
    for (i in which(true_counts_1cell > 0)) {
      x = recovered_vec[(GI[i] + 1):GI[i + 1]];
      UMI_counts[i] = sum(x > 0);
    }

    return(list(UMI_counts, sequenced_vec, sum(frag_vec > 0)))
  }
}

#' sample from truncated normal distribution
#' @param a the minimum value allowed
#' @param b the maximum value allowed
rnorm_trunc <- function(n, mean, sd, a, b) {
  vec1 <- rnorm(n, mean = mean, sd = sd)
  beyond_idx <- which(vec1 < a | vec1 > b)
  if (length(beyond_idx) > 0) { # for each value < rate_2cap_lb
    substi_vec <- sapply(1:length(beyond_idx), function(i) {
      while (TRUE) {
        temp <- rnorm(1, mean = mean, sd = sd)
        if (temp > a | temp > b) { break } }
      return(temp) })
    vec1[beyond_idx] <- substi_vec
  }
  return(vec1)
}

#' Compute value of impulse function given parameters.
#' Enforces lower bound on value of function to avoid numerical
#' errors during model fitting.
#' @param vecImpulseParam (numeric vector number of impulse model parameters)
#' \{beta, h0, h1, h2, t1, t2\}
#' Vector of impulse model parameters.
#' @param vecTimepoints (numeric vector length number of time points)
#' Time points to be evaluated.
#'
#' @return vecImpulseValue (vec number of vecTimepoints)
#' Model values for given time points.
#'
#' @author David Sebastian Fischer
evalImpulse <- function(vecImpulseParam, vecTimepoints) {
  # beta is vecImpulseParam[1] h0 is vecImpulseParam[2] h1 is
  # vecImpulseParam[3] h2 is vecImpulseParam[4] t1 is vecImpulseParam[5]
  # t2 is vecImpulseParam[6]
  vecImpulseValue <- sapply(vecTimepoints, function(t) {
    (1 / vecImpulseParam[3]) *
      (vecImpulseParam[2] + (vecImpulseParam[3] - vecImpulseParam[2]) *
        (1 / (1 + exp(-vecImpulseParam[1] * (t - vecImpulseParam[5]))))) *
      (vecImpulseParam[4] + (vecImpulseParam[3] - vecImpulseParam[4]) *
        (1 / (1 + exp(vecImpulseParam[1] * (t - vecImpulseParam[6])))))
  })
  vecImpulseValue[vecImpulseValue < 10^(-10)] <- 10^(-10)

  return(vecImpulseValue)
}


#' Creating an example tree with 5 tips
#' @param plotting True for plotting the tree on console, False for no plot
#' @return a tree object
#' @export
Phyla5 <- function(plotting = F) {
  phyla <- rtree(2)
  phyla <- compute.brlen(phyla, 1)
  tip <- compute.brlen(phyla, 1)
  phyla <- bind.tree(phyla, tip, 1)
  phyla <- bind.tree(phyla, tip, 2)
  phyla <- bind.tree(phyla, tip, 2)
  phyla <- compute.brlen(phyla, c(1, 1, 1, 1, 1, 0.2, 0.2, 3))
  edges <- cbind(phyla$edge, phyla$edge.length)
  edges <- cbind(c(1:length(edges[, 1])), edges)
  connections <- table(c(edges[, 2], edges[, 3]))
  root <- as.numeric(names(connections)[connections == 2])
  tips <- as.numeric(names(connections)[connections == 1])
  phyla$tip.label <- as.character(tips)
  if (plotting == T) {
    plot(phyla, show.tip.label = F, lwd = 2)
    tiplabels(cex = 2)
    nodelabels(cex = 2)
  }
  return(phyla)
}

#' Creating an example tree with 3 tips
#' @param plotting True for plotting the tree on console, False for no plot
#' @return a tree object
#' @export
Phyla3 <- function(plotting = F) {
  # par(mfrow=c(2,2))
  phyla <- rtree(2)
  phyla <- compute.brlen(phyla, 1)
  tip <- compute.brlen(phyla, 1)
  phyla <- bind.tree(phyla, tip, 1)
  phyla <- compute.brlen(phyla, c(1, 1, 1, 2))
  edges <- cbind(phyla$edge, phyla$edge.length)
  edges <- cbind(c(1:length(edges[, 1])), edges)
  connections <- table(c(edges[, 2], edges[, 3]))
  root <- as.numeric(names(connections)[connections == 2])
  tips <- as.numeric(names(connections)[connections == 1])
  phyla$tip.label <- as.character(tips)

  if (plotting == T) {
    plot(phyla, show.tip.label = F, lwd = 2)
    tiplabels(cex = 2)
    nodelabels(cex = 2)
  }
  return(phyla)
}

#' Generating EVFs for cells sampled along the trajectory of cell development
#' @param phyla tree for cell developement
#' @param ncells number of cells
#' @param n_nd_evf Number of EVFs that do not have an impulse signal
#' @param n_de_evf Number of EVFs with an impulse signal
#' @param impulse if the impluse model should be used instead of Brownian motion
#' @param evf_center the mean of Gaussain function where the non-Diff EVFs are sampled from
#' @param vary which parameters are affected by Diff-EVFs. Can be "kon", "koff", "s", "all", "except_kon", "except_koff", "except_s". Suggestions are "all" or "s"
#' @param tip The leaf that the path with impulse lead to
#' @param Sigma The standard deviation of the brownian motion of EVFs changing along the tree
#' @param plotting Whether to plot the trajectory or not
#' @param plotname The string to be used in the output file name
#' @param seed the random seed
#' @return a list of two object, one is the evf, and the other is a dataframe indicating the branch each cell comes from (pop) and its depth in the tree (depth)
ContinuousEVF <- function(phyla, ncells, n_nd_evf, n_de_evf, impulse = F, evf_center = 1, vary = 's',
                          Sigma, plotting = T, plotname = 'cont_evf.pdf', seed, num_regulators) {
  set.seed(seed)
  edges <- cbind(phyla$edge, phyla$edge.length)
  edges <- cbind(c(1:length(edges[, 1])), edges)
  connections <- table(c(edges[, 2], edges[, 3]))
  root <- as.numeric(names(connections)[connections == 2])
  tips <- as.numeric(names(connections)[connections == 1])
  internal <- as.numeric(names(connections)[connections == 3])
  if (vary == 'all') {
    N_DE_evfs = c(n_de_evf, n_de_evf, n_de_evf)
    N_ND_evfs = c(n_nd_evf, n_nd_evf, n_nd_evf)
  }else if (vary == 'kon') {
    N_DE_evfs = c(n_de_evf, 0, 0)
    N_ND_evfs = c(n_nd_evf, n_de_evf + n_nd_evf, n_de_evf + n_nd_evf)
  }else if (vary == 'koff') {
    N_DE_evfs = c(0, n_de_evf, 0)
    N_ND_evfs = c(n_de_evf + n_nd_evf, n_nd_evf, n_de_evf + n_nd_evf)
  }else if (vary == 's') {
    N_DE_evfs = c(0, 0, n_de_evf)
    N_ND_evfs = c(n_nd_evf + n_de_evf, n_nd_evf + n_de_evf, n_nd_evf)
  }else if (vary == 'except_kon') {
    N_DE_evfs = c(0, n_de_evf, n_de_evf)
    N_ND_evfs = c(n_nd_evf + n_de_evf, n_nd_evf, n_nd_evf)
  }else if (vary == 'except_koff') {
    N_DE_evfs = c(n_de_evf, 0, n_de_evf)
    N_ND_evfs = c(n_nd_evf, n_de_evf + n_nd_evf, n_nd_evf)
  }else if (vary == 'except_s') {
    N_DE_evfs = c(n_de_evf, n_de_evf, 0)
    N_ND_evfs = c(n_nd_evf, n_nd_evf, n_nd_evf + n_de_evf)
  }
  neutral <- SampleSubtree(root, 0, evf_center, edges, ncells)
  param_names <- c("kon", "koff", "s")
  evfs <- lapply(c(1:3), function(parami) {
    nd_evf <- lapply(c(1:N_ND_evfs[parami]), function(ievf) {
      rnorm(ncells, evf_center, Sigma)
    })
    nd_evf <- do.call(cbind, nd_evf)
    if (N_DE_evfs[parami] != 0) {
      #if there is more than 1 de_evfs for the parameter we are looking at
      if (impulse == T) {
        pdf(file = plotname, width = 15, height = 5)
        tip <- rep(tips, ceiling(N_DE_evfs[parami] / length(tips)))
        de_evf <- lapply(c(1:N_DE_evfs[parami]), function(evf_i) {
          impulse <- ImpulseEVFpertip(phyla, edges, root, tips, internal, neutral, tip[evf_i], Sigma, evf_center = evf_center)
          if (plotting == T) { PlotRoot2Leave(impulse, tips, edges, root, internal) }
          re_order <- match(
            apply(neutral[, c(1:3)], 1, function(X) { paste0(X, collapse = '_') }),
            apply(impulse[, c(1:3)], 1, function(X) { paste0(X, collapse = '_') }))
          return(impulse[re_order,])
        })
        dev.off()
      }else {
        de_evf <- lapply(c(1:N_DE_evfs[parami]), function(evf_i) {
          SampleSubtree(root, 0, evf_center, edges, ncells, neutral = neutral)
        })
      }
      de_evf <- lapply(de_evf, function(X) { X[, 4] })
      de_evf <- do.call(cbind, de_evf)
      de_evf <- de_evf[c(1:ncells),]
      evfs <- cbind(nd_evf, de_evf)
      colnames(evfs) <- c(
        paste(param_names[parami], rep('nonDE', length(nd_evf[1,])), c(1:length(nd_evf[1,])), sep = '_'),
        paste(param_names[parami], rep('DE', length(de_evf[1,])), c(1:length(de_evf[1,])), sep = '_'))
    }else {
      evfs <- nd_evf
      colnames(evfs) <- paste(param_names[parami], rep('nonDE', length(nd_evf[1,])), c(1:length(nd_evf[1,])), sep = '_')
    }
    return(evfs)
  })
  meta <- data.frame(pop = apply(neutral[, c(1:2)], 1, function(X) { paste0(X, collapse = '_') }), depth = neutral[, 3])
  evf_all <- list(evfs, meta[c(1:ncells),], neutral)

  # add extra columns to make evf and gene_effects the same length
  evf_nd_reg <- lapply(1:2, function(parami) {
    evf_small <- lapply(1:num_regulators, function(ievf) {
      rnorm(ncells, evf_center, Sigma)
    })
    evf_small <- do.call(cbind, evf_small)
    colnames(evf_small) <- paste(param_names[parami], rep('reg', length(evf_small[1,])), seq_along(evf_small[1,]), sep = '_')
    return(evf_small)
  })
  evf_all[[1]][[1]] <- cbind(evf_all[[1]][[1]], evf_nd_reg[[1]])
  evf_all[[1]][[2]] <- cbind(evf_all[[1]][[2]], evf_nd_reg[[2]])

  return(evf_all)
  # note previously the number of sampled evfs and meta isn't necessarily ncells?
}

#' Generating EVFs for cells sampled from tip populations from a tree
#' @param phyla tree for cell developement
#' @param ncells_total number of cells from all populations
#' @param min_popsize size of the rarest population
#' @param i_minpop to specify which population has the smallest size
#' @param Sigma The standard deviation of the brownian motion of EVFs changing along the tree
#' @param n_nd_evf number of non-Diff EVFs
#' @param n_de_evf number of Diff EVFs
#' @param vary which parameters are affected by Diff-EVFs. Can be "kon", "koff", "s", "all", "except_kon", "except_koff", "except_s". Suggestions are "all" or "s"
#' @param evf_center the value used to generated evf means. Suggested value is 1
#' @param seed the random seed
#' @return a list of two object, one is the evf, and the other is a dataframe indicating the population each cell comes from (pop)
DiscreteEVF <- function(phyla, ncells_total, min_popsize, i_minpop, Sigma, n_nd_evf, n_de_evf,
                        vary, evf_center, seed) {
  set.seed(seed)
  npop <- length(phyla$tip.label)
  # set the number of cells in each population: first give each population min_popsize cells
  # then randomly distribute the rest of cells to all populations except the smallest one
  ncells_pop <- rep(min_popsize, npop)
  if (ncells_total < min_popsize * npop) {
    stop("The size of the smallest population is too big for the total number of cells") }
  larger_pops <- setdiff(1:npop, i_minpop)
  ncells_pop[larger_pops] <- floor((ncells_total - min_popsize) / length(larger_pops))
  leftover <- ncells_total - sum(ncells_pop)
  if (leftover > 0) {
    temp <- sample(larger_pops, leftover, replace = F); ncells_pop[temp] <- ncells_pop[temp] + 1
  }

  vcv_evf_mean <- vcv.phylo(phyla, cor = T)
  param_names <- c("kon", "koff", "s")
  if (vary == 'all') {
    N_DE_evfs = c(n_de_evf, n_de_evf, n_de_evf)
    N_ND_evfs = c(n_nd_evf, n_nd_evf, n_nd_evf)
  }else if (vary == 'kon') {
    N_DE_evfs = c(n_de_evf, 0, 0)
    N_ND_evfs = c(n_nd_evf, n_de_evf + n_nd_evf, n_de_evf + n_nd_evf)
  }else if (vary == 'koff') {
    N_DE_evfs = c(0, n_de_evf, 0)
    N_ND_evfs = c(n_de_evf + n_nd_evf, n_nd_evf, n_de_evf + n_nd_evf)
  }else if (vary == 's') {
    N_DE_evfs = c(0, 0, n_de_evf)
    N_ND_evfs = c(n_nd_evf + n_de_evf, n_nd_evf + n_de_evf, n_nd_evf)
  }else if (vary == 'except_kon') {
    N_DE_evfs = c(0, n_de_evf, n_de_evf)
    N_ND_evfs = c(n_nd_evf + n_de_evf, n_nd_evf, n_nd_evf)
  }else if (vary == 'except_koff') {
    N_DE_evfs = c(n_de_evf, 0, n_de_evf)
    N_ND_evfs = c(n_nd_evf, n_de_evf + n_nd_evf, n_nd_evf)
  }else if (vary == 'except_s') {
    N_DE_evfs = c(n_de_evf, n_de_evf, 0)
    N_ND_evfs = c(n_nd_evf, n_nd_evf, n_nd_evf + n_de_evf)
  }

  if (sum(N_DE_evfs) < 5) { warning("The number of DE evfs is less than 5; in the case of a small number of DE evfs, the structure of generated data
	                       may not closely follow the input tree. One can either increase nevf or percent_DEevf to avoid this warning.") }

  evfs <- lapply(1:3, function(iparam) {
    if (N_ND_evfs[iparam] > 0) {
      pop_evf_nonDE <- lapply(c(1:npop), function(ipop) {
        evf <- sapply(c(1:(N_ND_evfs[iparam])), function(ievf) { rnorm(ncells_pop[ipop], evf_center, Sigma) })
        return(evf)
      })
      pop_evf_nonDE <- do.call(rbind, pop_evf_nonDE)
      colnames(pop_evf_nonDE) <- rep('nonDE', N_ND_evfs[iparam])
    } else { pop_evf_nonDE <- NULL }
    if (N_DE_evfs[iparam] > 0) {
      pop_evf_mean_DE <- mvrnorm(N_DE_evfs[iparam], rep(evf_center, npop), vcv_evf_mean)
      pop_evf_DE <- lapply(c(1:npop), function(ipop) {
        evf <- sapply(c(1:N_DE_evfs[iparam]), function(ievf) { rnorm(ncells_pop[ipop], pop_evf_mean_DE[ievf, ipop], Sigma) })
        return(evf)
      })
      pop_evf_DE <- do.call(rbind, pop_evf_DE)
      colnames(pop_evf_DE) <- rep('DE', N_DE_evfs[iparam])
    } else { pop_evf_DE <- NULL }

    evfs_per_param <- cbind(pop_evf_nonDE, pop_evf_DE)
    colnames(evfs_per_param) <- sprintf("%s_%s_evf%d", param_names[iparam], colnames(evfs_per_param),
                                        1:(N_ND_evfs[iparam] + N_DE_evfs[iparam]))
    return(evfs_per_param)
  })
  meta <- data.frame(pop = do.call(c, lapply(c(1:npop), function(i) { rep(i, ncells_pop[i]) })))
  return(list(evfs, meta))
}


#' Simulate observed count matrix given technical biases and the true counts
#' @param true_counts gene cell matrix
#' @param meta_cell the meta information related to cells, will be combined with technical cell level information and returned
#' @param protocol a string, can be "nonUMI" or "UMI"
#' @param alpha_mean the mean of rate of subsampling of transcripts during capture step, default at 10 percent efficiency
#' @param alpha_sd the std of rate of subsampling of transcripts
#' @param lenslope amount of length bias
#' @param nbins number of bins for gene length
#' @param gene_len a vector with lengths of all genes
#' @param amp_bias_limit range of amplification bias for each gene, a vector of length ngenes
#' @param rate_2PCR PCR efficiency, usually very high, default is 0.8
#' @param nPCR1 the number of PCR cycles in "pre-amplification" step, default is 16
#' @param nPCR2 the number of PCR cycles used after fragmentation.
#' @param LinearAmp if linear amplification is used for pre-amplification step, default is FALSE
#' @param LinearAmp_coef the coeficient of linear amplification, that is, how many times each molecule is amplified by
#' @param depth_mean mean of sequencing depth
#' @param depth_sd std of sequencing depth
#' @import SummarizedExperiment
#' @export
True2ObservedCounts <- function(true_counts, meta_cell, protocol, alpha_mean = 0.1, alpha_sd = 0.002,
                                gene_len, depth_mean, depth_sd, lenslope = 0.02, nbins = 20,
                                amp_bias_limit = c(-0.2, 0.2),
                                rate_2PCR = 0.8, nPCR1 = 16, nPCR2 = 10, LinearAmp = F, LinearAmp_coef = 2000) {
  print('true_counts[, 1]')
  print(true_counts[, 1])
  ngenes <- dim(true_counts)[1]; ncells <- dim(true_counts)[2]
  amp_bias <- cal_amp_bias(lenslope, nbins, gene_len, amp_bias_limit)
  rate_2cap_lb <- 0.0005; depth_lb <- 200 # lower bound for capture efficiency and sequencing depth
  rate_2cap_vec <- rnorm_trunc(n = ncells, mean = alpha_mean, sd = alpha_sd, a = rate_2cap_lb, b = 1)
  depth_vec <- rnorm_trunc(n = ncells, mean = depth_mean, sd = depth_sd, a = depth_lb, b = Inf)
  observed_counts <- lapply(c(1:ncells), function(icell) {
    amplify_1cell(true_counts_1cell = true_counts[, icell], protocol = protocol,
                  rate_2cap = rate_2cap_vec[icell], gene_len = gene_len, amp_bias = amp_bias,
                  rate_2PCR = rate_2PCR, nPCR1 = nPCR1, nPCR2 = nPCR2, LinearAmp = LinearAmp,
                  LinearAmp_coef = LinearAmp_coef, N_molecules_SEQ = depth_vec[icell])
  })

  meta_cell2 <- data.frame(alpha = rate_2cap_vec, depth = depth_vec, stringsAsFactors = F)
  meta_cell <- cbind(meta_cell, meta_cell2)

  if (protocol == "UMI") {
    UMI_counts <- do.call(cbind, lapply(observed_counts, "[[", 1))
    nreads_perUMI <- lapply(observed_counts, "[[", 2)
    nUMI2seq <- sapply(observed_counts, "[[", 3)
    observed_counts <- UMI_counts
  } else
    observed_counts <- do.call(cbind, observed_counts)

  if (protocol == "UMI") { return(list(counts = observed_counts, cell_meta = meta_cell, nreads_perUMI = nreads_perUMI,
                                       nUMI2seq = nUMI2seq))
  } else
    return(list(counts = observed_counts, cell_meta = meta_cell))
}

#' Divide the observed counts into multiple batches by adding batch effect to each batch
#' @param observed_counts_res the output from True2ObservedCounts
#' @param nbatch number of batches
#' @param batch_effect_size amount of batch effects. Larger values result in bigger differences between batches. Default is 1.
#' @export
DivideBatches <- function(observed_counts_res, nbatch, batch_effect_size = 1) {
  ## add batch effects to observed counts
  # use different mean and same sd to generate the multiplicative factor for different gene in different batch
  observed_counts <- observed_counts_res[["counts"]]
  meta_cell <- observed_counts_res[["cell_meta"]]
  ncells <- dim(observed_counts)[2]; ngenes <- dim(observed_counts)[1]
  batchIDs <- sample(1:nbatch, ncells, replace = TRUE)
  meta_cell2 <- data.frame(batch = batchIDs, stringsAsFactors = F)
  meta_cell <- cbind(meta_cell, meta_cell2)

  mean_matrix <- matrix(0, ngenes, nbatch)
  gene_mean <- rnorm(ngenes, 0, 1)
  temp <- lapply(1:ngenes, function(igene) {
    return(runif(nbatch, min = gene_mean[igene] - batch_effect_size, max = gene_mean[igene] + batch_effect_size))
  })
  mean_matrix <- do.call(rbind, temp)

  batch_factor <- matrix(0, ngenes, ncells)
  for (igene in 1:ngenes) {
    for (icell in 1:ncells) {
      batch_factor[igene, icell] <- rnorm(n = 1, mean = mean_matrix[igene, batchIDs[icell]], sd = 0.01)
    }
  }
  observed_counts <- round(2^(log2(observed_counts) + batch_factor))
  return(list(counts = observed_counts, cell_meta = meta_cell))
}

#' Simulate technical biases 
#' @param lenslope amount of length bias. This value sould be less than 2*amp_bias_limit[2]/(nbins-1)
#' @param nbins number of bins for gene length
#' @param gene_len transcript length of each gene
#' @param amp_bias_limit range of amplification bias for each gene, a vector of length ngenes
cal_amp_bias <- function(lenslope, nbins, gene_len, amp_bias_limit) {

  ngenes <- length(gene_len)
  len_bias_bin <- (-c(1:nbins)) * lenslope
  len_bias_bin <- len_bias_bin - median(len_bias_bin)
  if (max(len_bias_bin) > amp_bias_limit[2]) {
    stop("The lenslope parameter is too large.")
  }
  max_rand_bias <- amp_bias_limit[2] - max(len_bias_bin)

  rand_bias <- rnorm(ngenes, mean = 0, sd = max_rand_bias)
  rand_bias[rand_bias > max_rand_bias] <- max_rand_bias
  rand_bias[rand_bias < -max_rand_bias] <- -max_rand_bias
  #rand_bias <- runif(ngenes, -max_rand_bias,  max_rand_bias)

  binsize <- floor(ngenes / nbins)
  genes_in_bins <- vector("list", nbins)
  bin4genes <- numeric(ngenes)
  for (ibin in 1:(nbins - 1)) {
    genes_in_bins[[ibin]] <- order(gene_len)[((ibin - 1) * binsize + 1):(ibin * binsize)]
    bin4genes[genes_in_bins[[ibin]]] <- ibin
  }
  genes_in_bins[[nbins]] <- order(gene_len)[((nbins - 1) * binsize + 1):ngenes]
  bin4genes[genes_in_bins[[nbins]]] <- nbins

  len_bias <- numeric(ngenes); len_bias <- len_bias_bin[bin4genes]
  amp_bias <- rand_bias + len_bias
  return(amp_bias)
}

#' expand transcript counts to a vector of binaries of the same length of as the number of transcripts
#' @param true_counts_1cell number of transcript in one cell
expand2binary <- function(true_counts_1cell) {
  print('true_counts_1cell')
  print(true_counts_1cell)
  expanded_vec <- rep(1, sum(true_counts_1cell))
  trans_idx <- sapply(which(true_counts_1cell > 0),
                      function(igene) { return(rep(igene, true_counts_1cell[igene])) })
  trans_idx <- unlist(trans_idx)
  return(list(expanded_vec, trans_idx))
}
