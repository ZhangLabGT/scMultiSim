library(zeallot)

#############################################################
# Master Equation Related Functions
#############################################################

#' Generate true results according to simulation parameters
#' @param GRN_params GRN_params is a matrix where: #    - column 1 is the target gene ID, #    - column 2 is the gene ID which acts as a transcription factor for the target (regulated) gene #    - column 3 is the effect of the column 2 gene ID on the column 1 gene ID
#' @param num_cells total number of cells from all populations
#' @param unregulated_to_regulated_gene_ratio ratio of unregulated genes to regulated genes.  The final number of genes is rounded up to the nearest 10
#' @param intrinsic_noise a 0 to 1 value which is the weight assigned to the random sample from the Beta-Poisson distribution, where the weight of the Beta-Poisson mean value is given a weight of 1 minus the intrinsic noise value.
#' @param Sigma parameter of the std of evf values within the same cell type, controls heterogeneity for each cell type
#' @param num_evfs number of EVFs for each kinetic parameter
#' @param diffEVF_fraction fraction of evfs which are differential evfs between cell types
#' @param phyla a tree which defines relationship between populations
#' @param beta splicing rate of each gene
#' @param d degradation rate of each gene
#' @param num_cycles for generating velocity data, the number of cycles run before sampling the gene expression of a cell
#' @param cycle_length for generating velocity data, a factor which is multiplied by the expected time to transition from kon to koff and back to to form the the length of a cycle
#' @param atac_effect a 0 to 1 value which is the weight assigned to the influence of chromatin accessability data on gene expression
#' @param randseed should produce same result if all other parameters are all the same
#' @param do_velocity set to TRUE to incorporate the RNA velocity simulation module to simulate unspliced and spliced RNA sequence counts and ground truth velocity
#' @param nregions_distribution the probability that a gene is regulated by respectively 0, 1, ..., (length(nregions_distribution) - 1) consecutive regions
#' @param p_zero the proportion of 0s we see in the ATAC-seq data
#' @param vary which kinetic parameters have differential evfs. Can be "all", "kon", "koff", "s", "except_kon", "except_koff", "except_s"
#' @param evf_center the value which evf mean is generated from (default=1)
#' @param impulse when generating continuous populations, use the impulse model or not. Default is FALSE
#' @param bimod adjusts the bimodality of gene expression, thus controlling intrinsic variation
#' @param geffect_mean the mean of gene effect size
#' @param gene_effect_prob probability of non-zero values in the gene effect vectors
#' @param gene_effects_sd controls differences between genes
#' @param reffect_mean the mean of region effect size
#' @param region_effect_prob probability of non-zero values in the region effect vectors
#' @param region_effects_sd controls differences between regions
#' @param scale_s the cell size parameter in (0,1). Use smaller value for cell types known to be small (like naive cells)
#' @return a list of eleven elements: 1. A matrix containing the true (spliced) transcript counts, 2. Gene level meta information, 3. Cell level meta information, including a matrix of EVFs and a vector of cell identity (for example, the population the cell belongs to), 4. The parameters k~on~, k~off~ and s used to simulation the true counts, 5. A matrix containing the true unspliced transcript counts, 6. The true RNA velocity information for each cell, 7. The pseudotime at which the cell counts were generated, 8. The scATAC-seq data, 9. The matrix region2gene, a 0 1 matrix of shape (nregions) x (num_genes), where a value of 1 indicates the gene is affected by the accessibility of the particular chromatin region, 10. Gene effects of each regulator gene (column) on every gene ID (row), 11. The number of genes used in the experiment
SimulateTrueCounts_spatial <- function(
  GRN_params, num_cells = 1000, unregulated_to_regulated_gene_ratio = 0.1, num_evfs = 500, diffEVF_fraction = 0.9, Sigma = 0.1, atac_effect = 0.5,
  beta = 0.4, d = 1, num_cycles = 3, cycle_length = 1, intrinsic_noise = 1, randseed = 0, do_velocity = FALSE, phyla = Phyla5(),
  nregions_distribution = cumsum(c(0.1, 0.5, 0.4)), p_zero = 0.8, vary = "s", evf_center = 1, impulse = F, bimod = 0, geffect_mean = 0,
  gene_effect_prob = 0.3, gene_effects_sd = 1, reffect_mean = 0, region_effect_prob = 0.3, region_effects_sd = 1, scale_s = 1,
  dyn_grn_params = NULL, spatial_params = NULL, mg_kon = NULL
) {
  num_diffEVF <- num_evfs * diffEVF_fraction
  set.seed(randseed)
  seed <- sample(1:1e5, size = 2)
  data(dens_nonzero)  # this is the density function of log(x+1), where x is the non-zero values for ATAC-SEQ data
  
  GRN_params <- NormalizeGRNParams(GRN_params)
  target_gene_ID_list <- sort(unique(GRN_params[, 1]))
  regulator_ID_list <- sort(unique(GRN_params[, 2]))
  
  num_target_genes <- length(target_gene_ID_list)
  num_regulators <- length(regulator_ID_list)
  num_GRN_genes <- num_target_genes + num_regulators
  num_genes <- ceiling((num_GRN_genes + num_GRN_genes * unregulated_to_regulated_gene_ratio) / 10) * 10
  
  # Spatial expr
  spatial_pr <- parse_spatial_params(spatial_params, num_genes)
  
  # get all possible paths in the tree; save them in `paths`
  phyla_info <- phyla_nodes(phyla)
  root <- phyla_info$root
  paths <- list()
  getPaths <- function(node, path) {
    path <- c(path, node)
    succ_edges <- phyla$edge[which(phyla$edge[, 1] == node),]
    if (nrow(succ_edges) == 0) {
      # append the path to `paths`
      paths <<- append(paths, list(path))
    } else {
      for (succ in succ_edges[, 2]) {
        getPaths(succ, path)
      }
    }
  }
  if (nrow(phyla$edge) == 1) {
    paths <- list(list(phyla$edge[1], phyla$edge[2]))
  } else {
    getPaths(root, list())
  }
  
  path_abs_len <- sapply(seq_along(paths), function(i) {
    path <- paths[[i]]
    len <- 0
    for (j in 1:(length(path)-1)) {
      parent <- path[[j]]
      child <- path[[j + 1]]
      len <- len + phyla_info$edges[phyla_info$edges[, 2] == parent & phyla_info$edges[, 3] == child, 4]
    }
    len
  })
  total_ncell <- ceiling((num_cells - 2) / max(path_abs_len) * sum(phyla$edge.length))
  if (total_ncell < num_cells) { total_ncell <- num_cells }
  
  gene_effects_by_regulator <- RegulatorGeneEffects(num_genes = num_genes, randseed = seed[2], GRN_params = GRN_params,
                                                    regulator_ID_list = regulator_ID_list, target_gene_ID_list = target_gene_ID_list)
  
  # the scATAC-Seq data cannot be generated randomly. Using the SymSim framework,
  # we are generating the scATAC-Seq data using a similar idea as the scRNA-Seq data,
  # where we consider EVFs as cell IDs and region IDs similar to gene effects.
  
  nregions <- length(nregions_distribution) * num_genes # regions 1 to nregion are considered sequentially located on the genome
  
  # region_effects is like the gene_effect for the ATAC-Seq feature, dimension nregions * nevf
  region_effect <- RegionEffects(nregions = nregions, nevf = num_evfs, randseed = seed[2], prob = region_effect_prob,
                                 effect_mean = reffect_mean, effect_sd = region_effects_sd)
  
  # to generate the EVF for ATAC-SEQ, we use the same EVF generation process as for s but with a different random seed, as s had DiffEVFs
  evf_ATAC_all <- ContinuousEVF(phyla, total_ncell, n_nd_evf = num_evfs - num_diffEVF, n_de_evf = num_diffEVF,
                                evf_center = evf_center, vary = vary, impulse = impulse,
                                Sigma = Sigma, plotting = F, seed = seed[1] + 1,
                                num_regulators = num_regulators + spatial_pr$num_regulators * spatial_pr$num_max_neighbours)
  atac_neutral <- evf_ATAC_all[[3]][1:total_ncell,]
  
  # get the length of each edge in the truncated atac_neutral, each row: (parent, child, len)
  edge_len <- t(apply(unique(atac_neutral[, 1:2]), 1, function(edge) {
    c(edge,
      sum(atac_neutral[, 1] == edge[1] & atac_neutral[, 2] == edge[2], na.rm = T))
  }))
  
  path_len <- list()
  path_evf_ATAC <- lapply(seq_along(paths), function(i) {
    path <- paths[[i]]
    len <- 0
    result <- NULL
    for (j in 1:(length(path)-1)) {
      parent <- path[[j]]
      child <- path[[j + 1]]
      idx <- atac_neutral[, 1] == parent & atac_neutral[, 2] == child
      result <- rbind(result, evf_ATAC_all[[1]][[3]][idx, ])
      len <- len + edge_len[edge_len[, 1] == parent & edge_len[, 2] == child, 3]
    }
    path_len[[i]] <<- len
    result
  })
  
  # first generate scATAC-Seq data, only after this we can generate scRNA-Seq data
  atacseq_data <- lapply(seq_along(paths), function(i) {
    ATAC_SEQ(path_evf_ATAC[[i]], region_effect, dens_nonzero, p_zero, randseed)
  })
  
  # get kon and koff & scale
  data(param_realdata.zeisel.imputed)
  match_params[, 1:3] <- log(base = 10, match_params[, 1:3])
  match_params_den <- lapply(1:3, function(i) {
    density(match_params[, i], n = 2000)
  })
  
  # generate region2gene, a 0 1 coupling matrix between regions and genes of shape (nregions) x (num_genes), where a value
  # of 1 indicates the gene is affected by a particular region
  region2gene <- Region2Gene(num_genes, nregions, nregions_distribution, randseed)
  
  # gene_effects are the gene effect values, and have shape (kinetic parameter) x (num_genes) x (num_evfs + num_regulators)
  gene_effects <- GeneEffects(num_genes = num_genes, num_evfs = num_evfs, randseed = seed[2], prob = gene_effect_prob, geffect_mean = geffect_mean,
                              geffect_sd = gene_effects_sd, num_diffEVF = num_diffEVF,
                              regulator_ID_list, target_gene_ID_list, gene_effects_by_regulator = gene_effects_by_regulator )
  
  for (i in 1:2) {
    gene_effects[[i]] <- cbind(gene_effects[[i]], spatial_pr$effects)
  }
  
  getParams <- function(evf_all, pathi) {
    atac <- atacseq_data[[pathi]]
    params <- Get_params(gene_effects, evf_all[1:2], match_params_den[1:2], bimod, scale_s = scale_s,
                         atacseq_data = t(atac), region2gene = region2gene, atac_effect = atac_effect)
    
    return(list(
      gene_effects = gene_effects,
      gene_effects_by_regulator = gene_effects_by_regulator,
      params = params
    ))
  }
  
  grn <- list(geff = gene_effects_by_regulator, regulators = regulator_ID_list, target_genes = target_gene_ID_list)
  
  evf_params <- list(
    evf_center = evf_center,
    Sigma = Sigma,
    n_nd_evf = num_evfs - num_diffEVF,
    n_de_evf = num_diffEVF
  )
  
  result <- RNA_SEQ_spatial(
    phyla, num_genes, num_cells, total_ncell, num_regulators, paths, path_len,
    evf_params,
    match_params_den, intrinsic_noise, scale_s,
    do_velocity, beta, d, num_cycles, cycle_length,
    grn, getParams, spatial_pr, seed, mg_kon)
  
  counts <- result$counts_s
  unspliced_counts <- result$counts_u
  cell_time <- result$cell_time
  velocity <- result$velocity
  params <- result$params
  gene_effects <- result$gene_effects
  
  cell_meta <- cbind(cellid = paste('cell', seq(1, num_cells), sep = '_'), result$meta)
  
  dots_list(counts, gene_effects, cell_meta, kinetic_params = params,
            unspliced_counts, velocity,
            cell_time, atacseq_data, region2gene,
            gene_effects_by_regulator, num_genes,
            dyn_grn = grn, grid = result$grid, counts_hist = result$counts_hist,
            curr_evf = result$curr_evf, curr_ligand_evf = result$curr_ligand_evf,
            cell_evf = result$cell_evf, geff = result$geff, s_hist=result$s_hist)
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
#' @param cycle_length_factor for generating velocity data, a factor which is multiplied by the expected time to transition from kon to koff and back to to form the the length of a cycle
#' @param cell the cell number currently having counts generated
#' @return a list of 4 elements, the first element is true counts, second is the gene level meta information, the third is cell level meta information, including a matrix of evf and a vector of cell identity, and the fourth is the parameters kon, koff and s used to simulation the true counts
#' @import phytools
gen_1branch_spatial <- function(kinet_params, start_state, start_s, start_u, randpoints1, ncells1, ngenes, beta_vec, d_vec, cycle_length_factor, cell) {
  
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
    if (is.null(start_s)) {
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
    curr_time <- numeric(nsteps)
    p_table <- matrix(0, 2, 2)
    p_table[1, 2] <- 1 / (npart * (1 / k_off / cycle_length)); p_table[1, 1] <- 1 - p_table[1, 2]
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
        p_table[1, 2] <- 1 / (npart * (1 / k_off / cycle_length)); p_table[1, 1] <- 1 - p_table[1, 2]
        p_table[2, 1] <- 1 / (npart * (1 / k_on / cycle_length)); p_table[2, 2] <- 1 - p_table[2, 1]
      } else {
        which_cell[t] <- which_cell[t - 1]
      }
      
      if (runif(1, 0, 1) > p_table[curr_state[t - 1], 1]) {
        curr_state[t] <- 2
      } else {
        curr_state[t] <- 1
      }
      
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
      counts_u1[igene, icell] <- as.integer(x[all_idx[closest]])
      counts_s1[igene, icell] <- as.integer(y[all_idx[closest]])
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
#' @param gene_effects_by_regulator gene effects of each regulator gene (column) on every gene ID (row)
#' @param params the parameters for simulating gene expression from EVf and gene effects
#' @param intrinsic_noise a 0 to 1 value which is the weight assigned to the random sample from the Beta-Poisson distribution, where the weight of the Beta-Poisson mean value is given a weight of 1 minus the intrinsic noise value.
#' @param Sigma parameter of the std of evf values within the same population
#' @param phyla the cell developmental tree if chosing 'discrete' or 'continuous' evf type. Can either be generated randomly (using pbtree(nclusters) function from phytools package) or read from newick format file using the ape package
#' @param scale_s a factor to scale the s parameter, which is used to tune the size of the actual cell (small cells have less number of transcripts in total)
#' @param evf_all a list of two object, one is the evf, and the other is a dataframe indicating the branch each cell comes from (pop) and its depth in the tree
#' @param do_velocity set to TRUE to incorporate the RNA velocity simulation module to simulate unspliced and spliced RNA sequence counts and ground truth velocity
#' @param beta splicing rate of each gene
#' @param d degradation rate of each gene
#' @param num_cycles for generating velocity data, the number of cycles run before sampling the gene expression of a cell
#' @param cycle_length for generating velocity data, a factor which is multiplied by the expected time to transition from kon to koff and back to to form the the length of a cycle
#' @param evf_center mean value of normal distribution used to generate evf start
#' @param match_params_den density parameters that data is matched to
#' @return generate a true RNA SEQ count matrix, which has the shape (genes) x (cells)

RNA_SEQ_spatial <- function(
  phyla, num_genes, num_cells, num_total_cells, num_regulators,
  paths, path_len,
  evf_params,
  match_params_den, intrinsic_noise, scale_s,
  do_velocity, beta, d, num_cycles, cycle_length,
  grn, geff_getter, lig, seed, mg_kon = NULL
) {
  num_ligand_cif <- lig$num_regulators * lig$num_max_neighbours
  num_extra_cif <- num_regulators + num_ligand_cif
  
  root <- phyla_nodes(phyla)$root
  n_paths <- length(paths)
  n_steps <- max(unlist(path_len))
  cell_path <- sample(1:n_paths, num_cells, replace = T)
  
  cat(sprintf('total steps: %d\n', n_steps))
  
  grid <- CreateSpatialGrid(num_cells)
  
  # get evf
  cat('get evfs...\n')
  cell_evf <- ContinuousEVF_spatial(
    phyla, num_cells, cell_path, paths, path_len,
    num_total_cells, evf_params$n_nd_evf, evf_params$n_de_evf,
    impulse = F, evf_center = evf_params$evf_center, vary = 's',
    evf_params$Sigma, seed, num_extra_cif
  )
  
  # GRN
  cat('get params...\n')
  geff <- lapply(1:num_cells, function(i) geff_getter(cell_evf[[i]]$evfs, cell_path[[i]]))
  
  
  # cells one by one
  #evf_reg <- matrix(0, nrow = num_regulators, ncol = num_cells)
  curr_evf <- lapply(1:num_cells, function(i) rnorm(num_regulators, evf_params$evf_center, evf_params$Sigma))
  curr_ligand_evf <- lapply(1:num_cells, function(i) rnorm(lig$num_regulators, evf_params$evf_center, evf_params$Sigma))
  ## cellcount <- 1
  ## params[[3]] <- matrix(nrow = nrow(params[[1]]), ncol = ncol(params[[1]]))
  meta <- data.frame(pop = character(), depth = double())
  
  num_cycles <- num_cycles
  d_vec <- rnorm(n = num_genes, mean = d, sd = 0.1) # degradation rate of the size ngenes * 1
  beta_vec <- rnorm(n = num_genes, mean = beta, sd = 0.1) # splicing rate of the size ngenes * 1
  counts_u <- matrix(nrow = num_genes, ncol = num_cells)
  counts_s <- matrix(nrow = num_genes, ncol = num_cells)
  counts_hist = array(NA, dim = c(num_cells, n_steps, num_genes))
  state_mat <- matrix(nrow = num_genes, ncol = num_cells)
  cell_time <- vector(length = num_cells)
  velocity <- matrix(nrow = num_genes, ncol = num_cells)
  # generate the promoter states of the root cell with dimension ngenes * 1
  root_state <- sample(c(1, 2), size = num_genes, replace = TRUE)
  
  s_hist = data.frame()
  
  cat('STEP = ')
  for (t in 1:n_steps) {
    if (t > num_cells) { next }
    cat(sprintf('%d..', t))
    # add new cell
    grid$allocate(t)
    # there are t cells at current step
    if (t == 498) {
      #browser()
    }
    for (icell in 1:t) {
      step <- t - icell + 1
      if (step > nrow(cell_evf[[icell]]$evfs[[3]])) {
        next
      }
      geff_cell <- geff[[icell]]
      gene_effects <- geff_cell$gene_effects
      gene_effects_by_regulator <- geff_cell$gene_effects_by_regulator
      params <- geff_cell$params
      # get s for a cell
      neighbours <- grid$get_neighbours(icell)
      ligand_cif <- double(num_ligand_cif)
      for (i in seq_along(neighbours)) {
        nb <- neighbours[i]
        if (is.na(nb)) { next }
        base <- (i - 1) * lig$num_regulators
        for (j in seq_along(lig$regulators)) {
          # !!! what's the value here
          ligand_cif[base + j] <- curr_ligand_evf[[nb]][j]
        }
      }
      # if (t == 10) browser()
      extra_cif <- c(curr_evf[[icell]], ligand_cif)
      s_cell <- cell_evf[[icell]]$evfs[[3]][step,] %*% t(gene_effects[[3]]) +
        extra_cif %*% t(cbind(gene_effects_by_regulator, lig$effects))
      # scale
      s0 = s_cell[101]
      temp <- alply(s_cell, 1, function(Y) { Y })
      values <- do.call(c, temp)
      ranks <- rank(values)
      sorted <- sort(SampleDen(nsample = max(ranks), den_fun = match_params_den[[3]]))
      s_cell <- matrix(data = sorted[ranks], ncol = length(s_cell))
      s1 = s_cell[101]
      s_cell <- t(apply(s_cell, 2, function(x) { x <- 10^x })) * scale_s
      s_cell <- t(s_cell)
      # counts
      counts_cell <- sapply(1:num_genes, function(i) {
        # k_on <- if (i == 105) mg_kon else params[[1]][i, step]
        k_on <- params[[1]][i, step]
        k_off <- params[[2]][i, step]
        y <- rbeta(1, k_on, k_off)
        x <- rpois(1, y * s_cell[i])
        y_mean <- k_on / (k_on + k_off)
        x_mean <- y_mean * s_cell[i]
        return(intrinsic_noise * x + (1 - intrinsic_noise) * x_mean)
      })
      
      s_hist <- rbind(s_hist, list(
        n1= extra_cif[7], n2= extra_cif[9],n3= extra_cif[11],n4= extra_cif[13],
        cif_sum = sum(extra_cif[c(7,9,11,13)]),
        other_sum = sum(sapply(1:t, function(i) { if (i %in% neighbours) 0 else curr_ligand_evf[[i]][1] })),
        s = s_cell[101],
        s_self = (cell_evf[[icell]]$evfs[[3]][step,] %*% t(gene_effects[[3]]))[101],
        s0 = s0, s1 = s1,
        s_extra = (extra_cif %*% t(cbind(gene_effects_by_regulator, lig$effects)))[101],
        k_on = params[[1]][101, step], k_off = params[[1]][105, step],
        count = counts_cell[101]
        ))
      
      counts_reg <- counts_cell[c(grn$regulators, lig$regulators)]
      evf_thiscell <- counts_reg / (counts_reg + mean(counts_cell))
      
      counts_s[, icell] <- counts_cell
      counts_hist[icell, t, ] <- counts_cell
      #evf_reg[, icell] <- evf_thiscell
      n_regu <- length(grn$regulators)
      curr_evf[[icell]] <- evf_thiscell[1:n_regu]
      curr_ligand_evf[[icell]] <- evf_thiscell[(n_regu+1):(n_regu+lig$num_regulators)]
      meta[icell, ] <- cell_evf[[icell]]$meta[step, ]
    }
  }
  cat('\n')
  
  return(list(
    counts_u = counts_u, counts_s = counts_s,
    cell_time = cell_time, velocity = velocity, params = params, meta = meta,
    grid = grid, counts_hist = counts_hist,
    curr_evf = curr_evf, curr_ligand_evf = curr_ligand_evf,
    cell_evf = cell_evf, geff = geff, s_hist=s_hist
  ))
}


#' @param phyla tree for cell developement
#' @param ncells number of cells
#' @param n_nd_evf Number of EVFs that do not have an impulse signal
#' @param n_de_evf Number of EVFs with an impulse signal
#' @param impulse if the impluse model should be used instead of Brownian motion
#' @param evf_center the mean of Gaussain function where the non-Diff EVFs are sampled from
#' @param vary which parameters are affected by Diff-EVFs. Can be "kon", "koff", "s", "all", "except_kon", "except_koff", "except_s". Suggestions are "all" or "s"
#' @param Sigma The standard deviation of the brownian motion of EVFs changing along the tree
#' @param plotting Whether to plot the trajectory or not
#' @param plotname The string to be used in the output file name
#' @param seed the random seed
#' @param num_regulators number of regulator genes
#' @return a list of two object, one is the evf, and the other is a dataframe indicating the branch each cell comes from (pop) and its depth in the tree (depth)
ContinuousEVF_spatial <- function(
  phyla, num_cells, path_i, paths, path_len,
  ncells, n_nd_evf, n_de_evf, impulse = F, evf_center = 1, vary = 's',
  Sigma, seed, num_regulators
) {
  set.seed(seed)
  phyla_info <- phyla_nodes(phyla)
  
  # num of de/nd evf for k_on/k_off/s
  n_evf <- n_de_evf + n_nd_evf
  is_vary <- switch(
    vary,
    'all'         = c(T, T, T),
    'kon'         = c(T, F, F),
    'koff'        = c(F, T, F),
    's'           = c(F, F, T),
    'except_kon'  = c(F, T, T),
    'except_koff' = c(T, F, T),
    'except_s'    = c(T, T, F)
  )
  N_DE_evfs <- sapply(is_vary, function(x) ifelse(x, n_de_evf, 0))
  N_ND_evfs <- n_evf - N_DE_evfs
  
  neutral <- SampleSubtree(phyla_info$root, 0, evf_center, phyla_info$edges, ncells, neutral = NA)
  neutral <- neutral[1:ncells,]
  param_names <- c("kon", "koff", "s")
  
  de_evfs <- list(NULL, NULL, NULL)
  
  lapply(1:num_cells, function(celli) {
    pathi <- path_i[[celli]]  
    path <- paths[[pathi]]
    ln <- length(path)
    ncells_on_path <- path_len[[pathi]]
    
    # each cell
    evfs <- lapply(1:3, function(parami) {
      param_name <- param_names[parami]
      n_nd_evf <- N_ND_evfs[parami]
      n_de_evf <- N_DE_evfs[parami]
      
      nd_evf <- lapply(1:n_nd_evf, function(ievf) {
        rnorm(ncells_on_path, evf_center, Sigma)
      })
      nd_evf <- do.call(cbind, nd_evf)
      
      # nd_evf only
      if (n_de_evf == 0) {
        evfs <- nd_evf
        colnames(evfs) <- paste(param_name, rep('nonDE', length(nd_evf[1,])), c(1:length(nd_evf[1,])), sep = '_')
        return(evfs) 
      }
      
      # there is more than 1 de_evfs for the parameter we are looking at
      if (is.null(de_evfs[[parami]])) {
        de_evf_all <- lapply(1:n_de_evf, function(evf_i) {
          # supply neutral to have the same t_sample values for all cells
          SampleSubtree(phyla_info$root, 0, evf_center, phyla_info$edges, ncells, neutral = neutral)
        })
        de_evfs[[parami]] <<- de_evf_all
      } else {
        de_evf_all <- de_evfs[[parami]]
      }
      
      de_evf <- lapply(de_evf_all, function(x) x[x[, 1] %in% path[1:ln-1] & x[, 2] %in% path[2:ln], 4])
      de_evf <- do.call(cbind, de_evf)
      
      # only take ncells rows since SampleSubTree's result is not in the exact size
      # de_evf <- de_evf[c(1:ncells),]
      if (ncells_on_path != nrow(de_evf)) {
        stop('ncells_on_path is not consistent')
      }
      
      evfs <- cbind(nd_evf, de_evf)
      colnames(evfs) <- c(
        paste(param_name, rep('nonDE', n_nd_evf), c(1:n_nd_evf), sep = '_'),
        paste(param_name, rep('DE', n_de_evf), c(1:n_de_evf), sep = '_')
      )
      
      return(evfs)
    })
    
    neutral <- neutral[neutral[, 1] %in% path[1:ln-1] & neutral[, 2] %in% path[2:ln], ]
    
    meta <- data.frame(
      pop = apply(neutral[, c(1:2)], 1, function(X) { paste0(X, collapse = '_') }),
      depth = neutral[, 3]
    )
    evf_all <- list(evfs = evfs, meta = meta[c(1:ncells_on_path),], neutral = neutral, num_nd = n_nd_evf)
    
    # add extra columns to make evf and gene_effects the same length
    for (parami in 1:2) {
      evf_small <- lapply(1:num_regulators, function(ievf) {
        rnorm(ncells_on_path, evf_center, Sigma)
      })
      evf_small <- do.call(cbind, evf_small)
      colnames(evf_small) <- paste(param_names[parami], rep('reg', length(evf_small[1,])), seq_along(evf_small[1,]), sep = '_')
      evf_all$evfs[[parami]] <- cbind(evf_all$evfs[[parami]], evf_small)
    }
    
    return(evf_all)
  })
}


SamplePath <- function(path, depth, anc_state, edges, ncells, neutral = NA) {
  for (j in 1:(length(path)-1)) {
    parent <- path[[j]]
    children <- path[[j + 1]]
    edge <- edges[edges[, 2] == parent & edges[, 3] == children,]
    t_sample <- neutral[neutral[, 1] == edge[2] & neutral[, 2] == edge[3], 3]
    result1 <- SampleEdge(edge, depth, anc_state, edges, ncells, t_sample)
    # print(sprintf('Path %d -> %d: %d', parent, children, nrow(result1)))
    if (j == length(path) - 1) {
      result1 <- result1[c(1:(length(result1[, 1] - 1))),]
    } 
    if (j == 1) {
      result <- result1 
    } else {
      result <- rbind(result, result1)
    }
  }
  
  return(result)
}

parse_spatial_params <- function(params, num_genes) {
  num_max_neighbours <- 4
  effects <- NULL
  regulators <- list()
  
  if (is.null(params)) {
    return(NULL)
  }
  
  spatial_list <- params$params
  regulators <- sort(unique(spatial_list[, 2]))
  num_regulators <- length(regulators)
  effects <- matrix(0, num_genes, num_regulators * num_max_neighbours)
  # e.g. 2 regulators => c(0, 2, 4, 8)
  rg_idx_base <- 0:(num_max_neighbours - 1) * num_regulators
  
  for (irow in 1:nrow(spatial_list)) {
    tg <- spatial_list[irow, 1]
    rg <- which(regulators %in% spatial_list[irow, 2])
    weight <- spatial_list[irow, 3]
    effects[tg, rg_idx_base + rg] <- weight
  }
  
  list(
    regulators = regulators,
    num_regulators = num_regulators,
    num_max_neighbours = num_max_neighbours,
    effects = effects
  )
}

sp_corr <- function(results, which_rg = 1) {
  reg_len <- 1 # length(lig_params$regulator)
  last_nb_means <- double(num_cells * 2 * reg_len)
  last_no_nb_means <- double(num_cells * 2 * reg_len)
  last_tg_all <- double(num_cells * 500 * reg_len)
  nb_means <- double(num_cells * 500 * reg_len)
  no_nb_means <- double(num_cells * 500 * reg_len)
  tg_all <- double(num_cells * 500 * reg_len)
  n <- 0
  m <- 0
  #set.seed(66)
  
  for (icell in 1:499) {
    neighbours <- results$grid$get_neighbours(icell)
    neighbours <- neighbours[!is.na(neighbours)]
    nbs <- neighbours
    non_nbs <- setdiff(1:500, nbs)
    if (length(nbs) == 0) next
    
    # last
    ct <- results$counts
    for (i in which_rg:which_rg) {
      rg <- lig_params$regulator[i]
      tg <- lig_params$target[i]
      m <- m + 1
      last_nb_means[m] <- mean(ct[rg, nbs])
      last_no_nb_means[m] <- mean(ct[rg, non_nbs])
      last_tg_all[m] <- ct[tg, icell]
    }
    
    # history
    ct <- results$counts_hist
    for (step in icell:499) {
      nbs <- setdiff(neighbours, icell:500)
      non_nbs <- setdiff(1:step, nbs)
      if (length(nbs) == 0) next
      for (i in which_rg:which_rg) {
        rg <- lig_params$regulator[i]
        tg <- lig_params$target[i]
        n <- n + 1
        nb_means[n] <- mean(ct[nbs, step, rg])
        no_nb_means[n] <- mean(ct[non_nbs, step, rg])
        tg_all[n] <- ct[icell, step, tg]
      }
    }
  }
  
  n <- n - 1
  m <- m - 1
  invisible(list(
    cor = list(
      last_n = cor(last_nb_means[1:m], last_tg_all[1:m]),
      last_nonn = cor(last_no_nb_means[1:m], last_tg_all[1:m]),
      n = cor(nb_means[1:n], tg_all[1:n]),
      nonn = cor(no_nb_means[1:n], tg_all[1:n])
    ),
    n = n, m = m,
    last_nb_means = last_nb_means[1:m], last_no_nb_means = last_no_nb_means[1:m], last_tg_all = last_tg_all[1:m],
    nb_means = nb_means[1:n], no_nb_means = no_nb_means[1:n], tg_all = tg_all[1:n]
  ))
}


sp_corr2 <- function(results, which_rg = 1) {
  reg_len <- 1 # length(lig_params$regulator)
  last_nb_means <- double(num_cells * 2 * reg_len)
  last_no_nb_means <- double(num_cells * 2 * reg_len)
  last_tg_all <- double(num_cells * 500 * reg_len)
  nb_means <- double(num_cells * 500 * reg_len)
  no_nb_means <- double(num_cells * 500 * reg_len)
  tg_all <- double(num_cells * 500 * reg_len)
  n <- 0
  m <- 0
  #set.seed(66)
  
  for (icell in 1:499) {
    neighbours <- results$grid$get_neighbours(icell)
    neighbours <- neighbours[!is.na(neighbours)]
    nbs <- neighbours
    non_nbs <- sample(setdiff(1:500, nbs), 4, replace = F)
    if (length(nbs) == 0) next
    
    # last
    ct <- results$counts
    for (i in which_rg:which_rg) {
      rg <- lig_params$regulator[i]
      tg <- lig_params$target[i]
      m <- m + 1
      last_nb_means[m] <- mean(ct[rg, nbs])
      last_no_nb_means[m] <- mean(ct[rg, non_nbs])
      last_tg_all[m] <- ct[tg, icell]
    }
    
    # history
    ct <- results$counts_hist
    for (step in icell:499) {
      nbs <- setdiff(neighbours, (step+1):500)
      non_nbs <- setdiff(1:step, nbs)
      non_nbs <- sample(non_nbs, max(length(nbs), 4), replace = F)
      if (length(nbs) == 0) next
      for (i in which_rg:which_rg) {
        rg <- lig_params$regulator[i]
        tg <- lig_params$target[i]
        n <- n + 1
        nb_means[n] <- mean(ct[nbs, step, rg])
        no_nb_means[n] <- mean(ct[non_nbs, step, rg])
        tg_all[n] <- ct[icell, step, tg]
      }
    }
  }
  
  n <- n - 1
  m <- m - 1
  invisible(list(
    cor = list(
      last_n = cor(last_nb_means[1:m], last_tg_all[1:m]),
      last_nonn = cor(last_no_nb_means[1:m], last_tg_all[1:m]),
      n = cor(nb_means[1:n], tg_all[1:n]),
      nonn = cor(no_nb_means[1:n], tg_all[1:n])
    )
    # n = n, m = m,
    # last_nb_means = last_nb_means[1:m], last_no_nb_means = last_no_nb_means[1:m], last_tg_all = last_tg_all[1:m],
    # nb_means = nb_means[1:n], no_nb_means = no_nb_means[1:n], tg_all = tg_all[1:n]
  ))
}
