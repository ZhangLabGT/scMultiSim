#' Simulate true scRNA and scATAC counts from the parameters
#'
#' @md
#' @param options See scMultiSim_help().
#' @param return_summarized_exp Whether to return a SummarizedExperiment object.
#'
#' @return scMultiSim returns an environment with the following fields:
#'
#' - `counts`: Gene-by-cell scRNA-seq counts.
#' - `atac_counts`: Region-by-cell scATAC-seq counts.
#' - `region_to_gene`: Region-by-gene 0-1 marix indicating the corresponding relationship between chtomatin regions and genes.
#' - `atacseq_data`: The "clean" scATAC-seq counts without added intrinsic noise.
#' - `cell_meta`: A dataframe containing cell type labels and pseudotime information.
#' - `cif`: The CIF used during the simulation.
#' - `giv`: The GIV used during the simulation.
#' - `kinetic_params`: The kinetic parameters used during the simulation.
#' - `.grn`: The GRN used during the simulation.
#' - `.grn$regulators`: The list of TFs used by all gene-by-TF matrices.
#' - `.grn$geff`: Gene-by-TF matrix representing the GRN used during the simulation.
#' - `.n`: Other metadata, e.g. `.n$cells` is the number of cells.
#'
#' If `do.velocity` is enabled, it has these additional fields:
#'
#'   - `unspliced_counts`: Gene-by-cell unspliced RNA counts.
#' - `velocity`: Gene-by-cell RNA velocity ground truth.
#' - `cell_time`: The pseudotime at which the cell counts were generated.
#'
#' If dynamic GRN is enabled, it has these additional fields:
#'
#'   - `cell_specific_grn`: A list of length `n_cells`. Each element is a gene-by-TF matrix, indicating the cell's GRN.
#'
#' If cell-cell interaction is enabled, it has these additional fields:
#'
#' - `grid`: The grid object used during the simulation.
#'   - `grid$get_neighbours(i)`: Get the neighbour cells of cell `i`.
#' - `cci_locs`: A dataframe containing the X and Y coordinates of each cell.
#' - `cci_cell_type_param`: A dataframe containing the CCI network ground truth: all ligand-receptor pairs between each pair of cell types.
#' - `cci_cell_types`: For continuous cell population, the sub-divided cell types along the trajectory used when simulating CCI.
#'
#' If it is a debug session (`debug = TRUE`), a `sim` field is available,
#' which is an environment contains all internal states and data structures.
#'
#' @export
#' @examples
#' data(GRN_params_100, envir = environment())
#' sim_true_counts(list(
#'   rand.seed = 0,
#'   GRN = GRN_params_100,
#'   num.cells = 100,
#'   num.cifs = 50,
#'   tree = Phyla5()
#' ))
sim_true_counts <- function(options, return_summarized_exp = FALSE) {
  # ==== options ===============================================================

  options <- .check_opt(options)
  phyla <- OP("tree")
  do_velocity <- OP("do.velocity")
  spatial_params <- options$cci
  is_discrete <- OP("discrete.cif")

  # debug?
  is_debug <- isTRUE(options$debug)

  # ==== initialization ========================================================

  cores <- OP("threads")
  if (cores != 1) {
    stop("Multithreading is not supported yet.")
  }

  # create session
  sim <- new.env()
  attr(sim, "name") <- "scMultiSim Session"
  sim$start_time <- Sys.time()

  # seeds
  # set.seed(OP("rand.seed"))
  seed <- sample(seq_len(1e5), size = 9)

  # get the GRN info and the numbers
  sim$do_spatial <- is.list(spatial_params)
  sim$is_dyn_grn <- is.list(OP("dynamic.GRN"))

  .grn_params <- OP("GRN")
  .sp_params <- if (sim$do_spatial) spatial_params$params else NULL
  c(.grn_params, .rn_sp) %<-% .renameGenes(sim, .grn_params, .sp_params)
  if (sim$do_spatial) {
    spatial_params$params <- .rn_sp
  }

  GRN <- .normalizeGRNParams(.grn_params)
  N <- .getNumbers(GRN, options)
  if (!is.null(GRN)) {
    GRN$geff <- .geneEffectsByRegulator(seed[1], GRN, N)
  }
  if (sim$is_dyn_grn) {
    dyngrn_opts <- .getDynGRNOpts(options$dynamic.GRN)
    GRN <- .CreateDynGRN(GRN, dyngrn_opts)
  }

  # velocity and spatial
  if (sim$do_spatial) {
    cat("CCI simulation is enabled.\n")

    .parseSpatialParams(spatial_params, N$gene, phyla, is_discrete) %->% c(
      sim$sp_params,
      sim$sp_regulators,
      sim$sp_targets,
      N$sp_regulators,
      N$max_nbs,
      sim$sp_effect,
      sim$sp_ctype_param,
      sim$cell_type_map,
      N$step_size,
      N$sp_del_lr_pair,
      same_type_prob,
      grid_size,
      sim$sp_layout,
      sim$sp_layout_param,
      sim$sp_sc_gt,
      sim$sp_static_steps,
      sim$sp_radius
    )

    sim$grid <- CreateSpatialGrid(
      N$cell, N$max_nbs,
      .grid.size = grid_size, .same.type.prob = same_type_prob,
      .method = sim$sp_layout, .method.param = sim$sp_layout_param,
      .nb.radius = sim$sp_radius
    )

    if (is_discrete) {
      N$max_layer <- N$cell
    } else {
      c(paths, total_ncell) %<-% .getPaths(N, options)
      sim$paths <- paths
      N$max_layer <- total_ncell
      sim$cell_path <- sample(seq_along(paths), OP("num.cells"), replace = TRUE)
    }
  }

  N$reg_cif <- if (is.null(GRN)) {
    0
  } else if (sim$do_spatial) {
    N$regulator + N$sp_regulators * N$max_nbs
  } else {
    N$regulator
  }

  # ==== simulation ============================================================

  sim$GRN <- GRN
  sim$N <- N
  sim$options <- options

  # 1.1 CIF
  if (!sim$do_spatial) {
    sim$CIF_all <- if (is_discrete) {
      sim$curr_cif <- "rna"
      .discreteCIF(seed[2], N, options, sim)
    } else {
      .continuousCIF(seed[2], N, options)
    }
  }

  # 1.2 RIV
  sim$RIV <- .regionIdentityVectors(seed[3], GRN, N, options)

  # 1.3 GIV
  sim$GIV <- .geneIdentifyVectors(seed[4], sim, options)

  # region-to-gene matrix
  sim$region_to_gene <- .regionToGeneMatrix(seed[5], N, options)

  # 1.4 ATAC-seq & CIF (for spatial)
  if (sim$do_spatial) {
    # ==== spatial ====
    CIF_atac_all <- if (is_discrete) {
      # discrete CIF: N$cell == number of layers
      .discreteCIF(seed[6], N, options, sim)
    } else {
      .continuousCIF(seed[6], N, options, ncell_key = "max_layer")
    }
    # get edge length
    atac_neutral <- CIF_atac_all$neutral[seq(N$max_layer), ]
    if (!is_discrete) {
      sim$path_len <- .getPathLen(atac_neutral, paths, N)
    }
    # atac
    sim$CIF_atac <- CIF_atac_all$cif$s
    .atacSeq(seed[7], sim)
    # CIF
    cat("Get CIF...")
    if (is_discrete) {
      cif <- .discreteCIFSpatial(seed[2], N, options, sim)
    } else {
      cif <- .continuousCIF(
        seed[2], N, options,
        is_spatial = TRUE,
        spatial_params = list(total_ncell, sim$paths, sim$cell_path, sim$path_len)
      )
      # add cell.type.idx
      if (!is.null(sim$cell_type_map)) {
        for (i in seq_along(cif$meta_by_path)) {
          cif$meta_by_path[[i]] <- cbind(
            cif$meta_by_path[[i]],
            data.frame(cell.type.idx = sim$cell_type_map[cif$meta_by_path[[i]]$cell.type])
          )
        }
      }
    }
    sim$CIF_spatial <- cif
  } else {
    # ==== not spatial ====
    sim$CIF_atac <- if (is_discrete) {
      sim$curr_cif <- "atac"
      .discreteCIF(seed[6], N, options, sim)$cif$s
    } else {
      .continuousCIF(seed[6], N, options)$cif$s
    }
    .atacSeq(seed[7], sim)
  }

  # 1.5 Params
  if (sim$do_spatial) {
    cat("Get params...")
    sim$params_spatial <- lapply(seq(N$cell), \(i) .getParams(
      seed[8] + i, sim,
      sp_cell_i = i, sp_path_i = sim$cell_path[i]
    ))
    cat("Done\n")
  } else {
    sim$params <- .getParams(seed[8], sim)
  }

  # check result
  if (is_debug) {
    .print_param_summary(sim)
  }

  # 1.6 RNA-seq
  if (sim$do_spatial) {
    .rnaSeqSpatial(seed[9], sim)
  } else {
    .rnaSeq(seed[9], sim)
  }

  .print_time(sim)

  # Results
  .getResult(sim, do_velocity, is_debug, return_summarized_exp)
}


# return the result of the simulation
.getResult <- function(sim, do_velocity, is_debug, return_summarized_exp) {
  cell_meta <- if (sim$do_spatial) {
    sim$meta_spatial
  } else {
    cbind(
      cell_id = paste0("cell", seq(sim$N$cell)),
      sim$CIF_all$meta
    )
  }
  rownames(cell_meta) <- paste0("cell", seq(sim$N$cell))

  if (is.null(sim$GRN)) {
    grn_params <- NULL
  } else {
    grn_params <- sim$GRN$params
    colnames(grn_params) <- c("target", "regulator", "effect")
    grn_params$regulator <- paste0("gene", grn_params$regulator)
    grn_params$target <- paste0("gene", grn_params$target)
  }

  # name other genes
  L <- length(sim$gene_name_map)
  N_other <- sim$N$gene - L
  sim$gene_name_map <- c(
    sim$gene_name_map,
    setNames(seq(N_other) + L, paste0("gene", seq(N_other) + L))
  )

  counts <- t(sim$counts_s)
  rownames(counts) <- names(sim$gene_name_map)
  colnames(counts) <- paste0("cell", seq(sim$N$cell))

  result <- list(
    counts = counts,
    cif = sim$CIF_all$cif,
    giv = sim$GIV,
    cell_meta = cell_meta,
    kinetic_params = sim$params,
    atacseq_data = if (sim$do_spatial) t(sim$sp_atac) else t(sim$atac_data),
    region_to_gene = sim$region_to_gene,
    num_genes = sim$N$gene,
    grn_params = grn_params,
    hge_scale = sim$hge_scale,
    .options = sim$options,
    .grn = sim$GRN,
    .n = sim$N
  )

  result$atac_counts <- .atacIntrNoise(result$atacseq_data)

  if (do_velocity) {
    result <- c(result, list(
      unspliced_counts = t(sim$counts_u),
      velocity = t(sim$velocity),
      cell_time = sim$cell_time
    ))
  }

  if (sim$is_dyn_grn) {
    result <- c(result, list(
      cell_specific_grn = lapply(
        seq(sim$N$cell), \(i) sim$GRN$history[[sim$dyngrn_ver_map[i]]]
      )
    ))
  }

  if (sim$do_spatial) {
    if (!is.null(sim$sp_ctype_param)) {
      ctype_param <- data.frame()
      for (i in seq(nrow(sim$sp_params))) {
        tg <- sim$sp_params[i, 1]
        rg <- sim$sp_params[i, 2]
        i_rg <- which(sim$sp_regulators %in% rg)
        for (ct1 in sim$cell_type_map) {
          for (ct2 in sim$cell_type_map) {
            if (sim$sp_ctype_param[ct1, ct2, i_rg] > 0) {
              ctype_param <- rbind(ctype_param, data.frame(
                ligand = rg, receptor = tg,
                effect = sim$sp_params[i, 3],
                ct1 = ct1, ct2 = ct2
              ))
            }
          }
        }
      }
      ctype_param$ligand <- names(sim$gene_name_map)[ctype_param$ligand]
      ctype_param$receptor <- names(sim$gene_name_map)[ctype_param$receptor]
      result$cif <- sim$CIF_spatial
    } else {
      ctype_param <- NULL
    }

    cci_locs <- do.call(rbind, sim$grid$locs)
    colnames(cci_locs) <- c("x", "y")
    rownames(cci_locs) <- paste0("cell", seq(sim$N$cell))

    result <- c(result, list(
      grid = sim$grid,
      cci_locs = cci_locs,
      cci_cell_type_param = ctype_param,
      cci_cell_types = sim$cell_type_map
    ))
    
    if (sim$sp_sc_gt) {
      result <- c(result, list(
        cci_gt = sim$cci_single_cell
      ))
    }
  }

  if (is_debug) {
    result <- c(result, list(
      sim = sim
    ))
  }

  if (return_summarized_exp) {
    return(.summarizeExp(result))
  }

  resenv <- new.env()
  attr(resenv, "name") <- "scMultiSim Result"

  for (n in names(result)) {
    resenv[[n]] <- result[[n]]
  }
  resenv
}


# normalize the gene names:
# - if the gene names are integers, convert them to strings
# - if there are additional unnamed genes (user specified num.genes), give them names
.renameGenes <- function(sim, grn_params, sp_params) {
  name_map <- integer()
  renamed_grn <- NULL
  renamed_sp <- NULL
  grn_genes <- NULL
  if (is.data.frame(grn_params)) {
    grn_genes <- sort(unique(c(grn_params[, 1], grn_params[, 2])))
    name_map <- setNames(seq_along(grn_genes), grn_genes)
    renamed_grn <- data.frame(
      target = name_map[grn_params[, 1]],
      regulator = name_map[grn_params[, 2]],
      effect = grn_params[, 3]
    )
  }
  if (is.data.frame(sp_params)) {
    sp_genes <- sort(unique(c(sp_params[, 1], sp_params[, 2])))
    sp_genes_only <- setdiff(sp_genes, grn_genes)
    name_map <- c(
      name_map,
      setNames(seq_along(sp_genes_only) + length(name_map), sp_genes_only)
    )
    renamed_sp <- data.frame(
      receptor = name_map[as.character(sp_params[, 1])],
      ligand = name_map[as.character(sp_params[, 2])],
      effect = sp_params[, 3]
    )
  }

  sim$gene_name_map <- name_map
  list(renamed_grn, renamed_sp)
}


#' Rename the original gene IDs in the GRN table to integers.
#' @param params GRN parameters.
#' @return list
.normalizeGRNParams <- function(params) {
  if (!is.data.frame(params)) {
    return(NULL)
  }

  if (nrow(unique(params[, seq_len(2)])) != nrow(params)) {
    stop("Duplicated edges found in the GRN.")
  }

  orig_gene_ids <- unique(c(params[, 1], params[, 2])) %>%
    sort() %>%
    as.character()
  new_gene_ids <- seq_along(orig_gene_ids)
  name_map <- setNames(new_gene_ids, orig_gene_ids)

  # rename gene ids
  for (j in c(1, 2)) {
    for (i in seq_len(nrow(params))) {
      params[i, j] <- name_map[[as.character(params[i, j])]]
    }
  }

  params[, 1] <- as.numeric(params[, 1])
  params[, 2] <- as.numeric(params[, 2])

  # get new target and regulator list
  targets <- sort(unique(params[, 1]))
  regulators <- sort(unique(params[, 2]))

  # return
  list(
    params = params,
    name_map = name_map,
    targets = targets,
    regulators = regulators,
    n_tgt = length(targets),
    n_reg = length(regulators)
  )
}


# construct the sim$N object, which contains the numbers of genes, CIFs, etc.
.getNumbers <- function(GRN, options) {
  N <- list()
  N$cell <- OP("num.cells")

  # gene
  if (is.null(GRN)) {
    N$grn.gene <- 0
    N$gene <- OP("num.genes")
    N$no.grn.gene <- N$gene
  } else {
    N$grn.gene <- length(GRN$name_map)
    N$gene <- if (is.numeric(options$num.genes) && options$num.genes >= N$grn.gene) {
      options$num.genes
    } else {
      ceiling((N$grn.gene + N$grn.gene * OP("unregulated.gene.ratio")) / 10) * 10
    }
    N$non.grn.gene <- N$gene - N$grn.gene
    N$regulator <- GRN$n_reg
    N$target <- GRN$n_tgt
  }

  # cif
  N$cif <- OP("num.cifs")
  n_diff_cif <- ceiling(N$cif * OP("diff.cif.fraction"))
  is_vary <- switch(OP("vary"),
    "all" = c(TRUE, TRUE, TRUE),
    "kon" = c(TRUE, FALSE, FALSE),
    "koff" = c(FALSE, TRUE, FALSE),
    "s" = c(FALSE, FALSE, TRUE),
    "except_kon" = c(FALSE, TRUE, TRUE),
    "except_koff" = c(TRUE, FALSE, TRUE),
    "except_s" = c(TRUE, TRUE, FALSE)
  )
  N$diff.cif <- vapply(
    is_vary, function(x) ifelse(x, n_diff_cif, 0),
    double(1)
  )
  N$nd.cif <- N$cif - N$diff.cif

  # regions
  N$region <- length(OP("region.distrib")) * N$gene

  # data: param density
  data(param_realdata.zeisel.imputed, envir = environment())
  match_params[, seq_len(3)] <- log(base = 10, match_params[, seq_len(3)])
  N$params_den <- lapply(seq_len(3), function(i) {
    density(match_params[, i], n = 2000)
  })

  # return
  N
}


# generate CIF for discrete cell population
.discreteCIF <- function(seed, N, options, sim) {
  phyla <- OP("tree")
  cif_center <- OP("cif.center")
  cif_sigma <- OP("cif.sigma")
  user_popsize <- OP("discrete.pop.size")
  min_popsize <- OP("discrete.min.pop.size")
  i_minpop <- OP("discrete.min.pop.index")

  npop <- length(phyla$tip.label)
  if (!is.null(sim$ncells_pop)) {
    ncells_pop <- sim$ncells_pop
  } else if (npop == 1) {
    ncells_pop <- N$cell
  } else if (is.integer(user_popsize)) {
    if (length(user_popsize) != npop) {
      stop("The number of discrete.pop.size must be equal to the total number of cell types.")
    }
    if (sum(user_popsize) != N$cell) {
      stop("The sum of discrete.pop.size must be equal to the total number of cells.")
    }
    ncells_pop <- user_popsize
  } else {
    ncells_pop <- rep(min_popsize, npop)
    if (N$cell < min_popsize * npop) {
      stop(sprintf(
        "The size of the smallest population (%g * %g) is too big for the total number of cells (%g)",
        min_popsize, npop, N$cell
      ))
    }

    larger_pops <- setdiff(seq(npop), i_minpop)
    ncells_pop[larger_pops] <- floor((N$cell - min_popsize) / length(larger_pops))
    leftover <- N$cell - sum(ncells_pop)
    if (leftover > 0) {
      temp <- sample(larger_pops, leftover, replace = FALSE)
      ncells_pop[temp] <- ncells_pop[temp] + 1
    }
  }

  if (is.null(sim$ncells_pop)) {
    sim$ncells_pop <- ncells_pop
  }

  vcv_evf_mean <- vcv.phylo(phyla, corr = TRUE)
  param_name <- c("kon", "koff", "s")

  evfs <- lapply(seq_len(3), function(iparam) {
    n_nd_cif <- N$nd.cif[iparam]
    n_diff_cif <- N$diff.cif[iparam]

    # ========== nd_cif ==========
    if (n_nd_cif > 0) {
      pop_evf_nonDE <- lapply(seq(npop), function(ipop) {
        evf <- vapply(seq(n_nd_cif), function(ievf) {
          rnorm(ncells_pop[ipop], cif_center, cif_sigma)
        }, double(ncells_pop[ipop]))
        return(evf)
      })
      pop_evf_nonDE <- do.call(rbind, pop_evf_nonDE)
      colnames(pop_evf_nonDE) <- rep("nonDE", n_nd_cif)
    } else {
      pop_evf_nonDE <- NULL
    }

    # if (iparam == 3) {
    #   if (sim$curr_cif == "rna") {
    #     pop_evf_nonDE[, seq_len(10)] <- ex_rna_cif
    #     message("rna cif")
    #   } else if (sim$curr_cif == "atac") {
    #     pop_evf_nonDE[, seq_len(10)] <- ex_atac_cif
    #     message("atac_cif")
    #   } else {
    #     stop("? cif")
    #   }
    # }

    # ========== de_cif ==========
    if (n_diff_cif > 0) {
      pop_evf_mean_DE <- MASS::mvrnorm(n_diff_cif, rep(cif_center, npop), vcv_evf_mean)
      pop_evf_DE <- lapply(seq(npop), function(ipop) {
        evf <- vapply(seq(n_diff_cif), function(ievf) {
          rnorm(ncells_pop[ipop], pop_evf_mean_DE[ievf, ipop], cif_sigma)
        }, double(ncells_pop[ipop]))
        return(evf)
      })
      pop_evf_DE <- do.call(rbind, pop_evf_DE)
      colnames(pop_evf_DE) <- rep("DE", n_diff_cif)
    } else {
      pop_evf_DE <- NULL
    }

    cif <- cbind(pop_evf_nonDE, pop_evf_DE)
    colnames(cif) <- sprintf(
      "%s_%s_cif%d", param_name[iparam], colnames(cif),
      seq(n_nd_cif + n_diff_cif)
    )

    # ========== generate reg_cif for k_on, k_off ===========
    if (iparam <= 2 && N$reg_cif > 0) {
      reg_cif <- lapply(
        seq(N$reg_cif),
        \(.) rnorm(N$cell, cif_center, cif_sigma)
      ) %>% do.call(cbind, .)
      colnames(reg_cif) <- paste(param_name, "reg", seq(N$reg_cif), sep = "_")
      cif <- cbind(cif, reg_cif)
    }
    return(cif)
  })

  names(evfs) <- param_name
  meta <- data.frame(pop = do.call(c, lapply(seq(npop), function(i) {
    rep(i, ncells_pop[i])
  })))

  list(cif = evfs, meta = meta)
}


#' Generates cifs for cells sampled along the trajectory of cell development
#'
#' @param seed random seed
#' @param N the number list
#' @param options the option list
#' @param ncell_key the key for the number of cells in N
#' @param is_spatial return a list of cifs for spatial
#' @param spatial_params the spatial parameters
#' @param .plot save the CIF plot
#' @param .plot.name plot name
#' @return a list containing the cif and meta data
.continuousCIF <- function(seed, N, options, ncell_key = "cell", is_spatial = FALSE, spatial_params = NULL,
                           .plot = FALSE, .plot.name = "cont_cif.pdf") {
  # set.seed(seed)

  ncells <- N[[ncell_key]]
  phyla <- OP("tree")
  cif_center <- OP("cif.center")
  cif_sigma <- OP("cif.sigma")
  use_impulse <- OP("use.impulse")
  tree_info <- .tree_info(phyla)
  neutral <- SampleSubtree(
    tree_info$root, 0, OP("cif.center"),
    tree_info$edges,
    if (is_spatial) N$max_layer else ncells,
    N$step_size,
    neutral = NA
  )


  # ==== generate cif: cell x n_cif ============================================

  cif <- .continuousCIFParam(
    is_spatial,
    ncells, N$nd.cif, N$diff.cif, N$reg_cif,
    cif_center, cif_sigma, N$step_size,
    neutral, phyla, tree_info,
    use_impulse,
    sp_params = spatial_params
  )

  # ===== metadata & output ====================================================

  if (is_spatial) {
    # metadata is already returned from .continuousCIFParamSpatial
    c(cif, list(neutral = neutral))
  } else {
    meta <- data.frame(
      pop = apply(neutral[, seq_len(2)], 1, \(X) paste0(X, collapse = "_")),
      depth = neutral[, 3]
    )[seq(ncells), ]

    list(cif = cif, meta = meta, neutral = neutral)
  }
}


# gene x regulator matrix for the GRN
.geneEffectsByRegulator <- function(seed, GRN, N) {
  # set.seed(seed)

  geff <- matrix(0L, nrow = N$grn.gene, ncol = N$regulator)
  for (r in seq_len(nrow(GRN$params))) {
    c(target, regulator, effect) %<-% GRN$params[r, ]
    regu_idx <- which(GRN$regulators %in% regulator)
    geff[target, regu_idx] <- effect
  }

  # add rows for non-GRN genes
  geff <- rbind(geff, matrix(0, nrow = N$non.grn.gene, ncol = N$regulator))
  colnames(geff) <- GRN$regulators
  rownames(geff) <- c(
    vapply(seq_len(N$grn.gene), \(i) names(GRN$name_map)[[i]], ""),
    seq_len(N$non.grn.gene) + N$grn.gene
  )
  geff
}


# generate GIV or RIV, return a size x cif matrix
.identityVectors <- function(size, n_cif, prob, mean, sd) {
  lapply(seq(size), function(i) {
    nonzero <- sample(c(0, 1),
      size = n_cif,
      prob = c(1 - prob, prob),
      replace = TRUE
    )
    nonzero[nonzero != 0] <- rnorm(sum(nonzero), mean, sd)
    nonzero
  }) %>% do.call(rbind, .)
}


# return a list of kon, koff (gene x cif+regu), s (gene x cif)
.geneIdentifyVectors <- function(seed, sim, options) {
  # set.seed(seed)

  GRN <- sim$GRN
  N <- sim$N

  # calculate initial gene effect values with original Symsim approach
  param_names <- c("kon", "koff", "s")
  giv <- lapply(seq_len(3), function(i) {
    .identityVectors(N$gene, N$cif,
      prob = OP("giv.prob"),
      mean = OP("giv.mean"),
      sd = OP("giv.sd")
    )
  })
  names(giv) <- param_names

  if (!is.null(GRN)) {
    # ==== k_on, k_off: all values are zero for regulators
    giv$kon <- cbind(giv$kon, matrix(0, N$gene, N$reg_cif))
    giv$koff <- cbind(giv$koff, matrix(0, N$gene, N$reg_cif))

    # ==== s: add the effect of regulator gene effects to gene_effects
    if (sim$do_spatial) {
      n_reg <- N$regulator + N$sp_regulators
      regu_list <- c(GRN$regulators, sim$sp_regulators)
      tgt_list <- c(GRN$targets, sim$sp_targets)
    } else {
      n_reg <- N$regulator
      regu_list <- GRN$regulators
      tgt_list <- GRN$targets
    }
    non_grn_gene <- setdiff(seq(N$gene), c(regu_list, tgt_list))
    giv$s[-non_grn_gene, ] <- 0
    regu_diff_cif <- matrix(0, n_reg, N$diff.cif[3])
    # regulator rows: gene effect of 2 added to two random differential cifs
    # for each master regulator gene
    indices <- replicate(n_reg, sample(seq(N$diff.cif[3]), 2, replace = FALSE)) %>% as.vector()
    regu_diff_cif[cbind(rep(seq(n_reg), each = 2), indices)] <- 2
    giv$s[regu_list, ] <- cbind(matrix(0, n_reg, N$nd.cif[3]), regu_diff_cif)
    # target rows: for every target gene, it should use the same gene_effect vector
    # as its regulators. if a gene has multiple regulators, its gene effects will
    # be the combination of that of the regulators.
    # (1 * reg) * (reg x cif)
    regu_counts <- rowSums(GRN$geff != 0)
    stopifnot(sum(regu_counts > 0) == GRN$n_tgt)
    grn_eff <- na.omit(GRN$geff %*% giv$s[GRN$regulators, ] /
      regu_counts /
      2)
    # ---
    rg_row <- which(rowSums(giv$s[GRN$targets, ]) > 0)
    grn_eff[rg_row, ] <- grn_eff[rg_row, ] * 0.5 + giv$s[GRN$targets[rg_row], ] * 0.5
    grn_target <- giv$s[GRN$targets, ] <- grn_eff
    # ---
    if (sim$do_spatial) {
      sp_eff <- sim$sp_effect[sim$sp_targets, seq(N$sp_regulators)]
      sp_target <- na.omit(sp_eff %*% giv$s[sim$sp_regulators, ] / 2)
      sp_factor <- mean(grn_target[grn_target > 0]) / mean(sp_target[sp_target > 0])
      if (is.nan(sp_factor)) {
        sp_factor <- 0
      }
      sp_eff <- (sp_target * sp_factor)
      # ---
      rg_row <- which(rowSums(giv$s[sim$sp_targets, ]) > 0)
      sp_eff[rg_row, ] <- sp_eff[rg_row, ] * 0.5 + giv$s[sim$sp_targets[rg_row], ] * 0.5
      giv$s[sim$sp_targets, ] <- sp_eff
      # ---
    }
    # remove small values
    giv$s[tgt_list, ][abs(giv$s[tgt_list, ]) < 0.2] <- 0
    # increase the scale of gene-effects for the new genes which do not have regulators
    giv$s[non_grn_gene, ] <- giv$s[non_grn_gene, ] * 3
  }

  # return
  giv
}


# return region x cif matrix
.regionIdentityVectors <- function(seed, GRN, N, options) {
  # set.seed(seed)
  .identityVectors(N$region, N$cif,
    prob = OP("riv.prob"),
    mean = OP("riv.mean"),
    sd = OP("riv.sd")
  )
}


# return region x gene matrix
.regionToGeneMatrix <- function(seed, N, options) {
  # set.seed(seed)

  res <- matrix(0, N$region, N$gene)
  # gene is regulated by 0, 1, or 2 regions
  regu_by <- sample(c(0, 1, 2), size = N$gene, replace = TRUE, prob = OP("region.distrib"))
  regu_by_1 <- which(regu_by == 1)
  regu_by_2 <- which(regu_by == 2)
  # for genes regulated by 1 region, select a random region
  res[cbind(sample(seq(N$region), length(regu_by_1), replace = TRUE), regu_by_1)] <- 1
  # for genes regulated by 2 regions, select 2 consecutive random regions
  region_idx <- sample(seq(N$region - 1), length(regu_by_2), replace = TRUE) %>% c(., . + 1)
  res[cbind(region_idx, rep(regu_by_2, 2))] <- 1

  # return
  res
}


.summarizeExp <- function (result) {
  se <- SummarizedExperiment(
    assays = list(
      counts = result$counts
    ),
    colData = result$cell_meta
  )
  result$atac_counts <- NULL
  result$cell_meta <- NULL

  metadata(se) <- result

  se
}
