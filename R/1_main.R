sim_true_counts <- function(options) {
  # ==== options ===============================================================

  options <- .check_opt(options)
  phyla <- OP(tree)
  do_velocity <- OP(do.velocity)
  spatial_params <- options$cci
  is_discrete <- OP(discrete.cif)

  # debug?
  is_debug <- isTRUE(options$debug)

  # ==== initialization ========================================================

  cores <- OP(threads)
  if (cores == 1) {
    registerDoSEQ()
  } else {
    library(doParallel)
    if (cores == 0) cores <- detectCores()
    registerDoParallel(cores = cores)
    cat(sprintf("Using %i threads for multithreading.\n", cores))
  }

  # create session
  sim <- new.env()
  attr(sim, "name") <- "scMultiSim Session"
  sim$start_time <- Sys.time()

  # seeds
  set.seed(OP(rand.seed))
  seed <- sample(1:1e5, size = 9)

  # get the GRN info and the numbers
  sim$is_dyn_grn <- is.list(OP(dynamic.GRN))
  GRN <- .normalize_GRN_params(OP(GRN))
  N <- .get_numbers(GRN, options)
  if (!is.null(GRN)) {
    GRN$geff <- .gene_effects_by_regulator(seed[1], GRN, N)
  }
  if (sim$is_dyn_grn) {
    dyngrn_opts <- .get_dyngrn_opts(options$dynamic.GRN)
    GRN <- .CreateDynGRN(GRN, dyngrn_opts)
  }

  # velocity and spatial
  sim$do_spatial <- is.list(spatial_params)
  if (sim$do_spatial) {
    cat("CCI simulation is enabled.\n")

    .parse_spatial_params(spatial_params, N$gene, phyla, is_discrete) %->% c(
      sim$sp_params,
      sim$sp_regulators,
      sim$sp_targets,
      N$sp_regulators,
      N$max_nbs,
      sim$sp_effect,
      sim$sp_ctype_param,
      sim$cell_type_map,
      N$step_size
    )
    sim$grid <- CreateSpatialGrid(N$cell, N$max_nbs)
    c(paths, total_ncell) %<-% .get_paths(N, options)
    sim$paths <- paths
    N$max_layer <- total_ncell
    sim$cell_path <- sample(1:length(paths), OP(num.cells), replace = T)
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
      .discrete_cif(seed[2], N, options)
    } else {
      .continuous_cif(seed[2], N, options)
    }
  }

  # 1.2 RIV
  sim$RIV <- .region_identity_vectors(seed[3], GRN, N, options)

  # 1.3 GIV
  sim$GIV <- .gene_identify_vectors(seed[4], sim, options)

  # region-to-gene matrix
  sim$region_to_gene <- .region_to_gene_matrix(seed[5], N, options)

  # 1.4 ATAC-seq & CIF (for spatial)
  if (sim$do_spatial) {
    CIF_atac_all <- .continuous_cif(seed[6], N, options, ncell_key = "max_layer")
    # get edge length
    atac_neutral <- CIF_atac_all$neutral[1:N$max_layer, ]
    sim$path_len <- .get_path_len(atac_neutral, paths, N)
    # atac
    sim$CIF_atac <- CIF_atac_all$cif$s
    .atac_seq(seed[7], sim)
    # CIF
    cat("Get CIF...")
    cif <- .continuous_cif(
      seed[2], N, options,
      is_spatial = T,
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
    sim$CIF_spatial <- cif
  } else {
    sim$CIF_atac <- if (is_discrete) {
      .discrete_cif(seed[6], N, options)$cif$s
    } else {
      .continuous_cif(seed[6], N, options)$cif$s
    }
    .atac_seq(seed[7], sim)
  }

  # 1.5 Params
  if (sim$do_spatial) {
    cat("Get params...")
    sim$params_spatial <- lapply(1:N$cell, \(i) .get_params(
      seed[8] + i, sim,
      sp_cell_i = i, sp_path_i = sim$cell_path[i]
    ))
    cat("Done\n")
  } else {
    sim$params <- .get_params(seed[8], sim)
  }
  
  # check result
  if (is_debug) {
    .print_param_summary(sim)
  }

  # 1.6 RNA-seq
  if (sim$do_spatial) {
    .rna_seq.spatial(seed[9], sim)
  } else {
    .rna_seq(seed[9], sim)
  }

  .print_time(sim)

  # Results
  .get_result(sim, do_velocity, is_debug)
}


.get_result <- function(sim, do_velocity, is_debug) {
  cell_meta <- if (sim$do_spatial) {
    sim$meta_spatial
  } else {
    cbind(
      cell_id = paste0("cell", 1:sim$N$cell),
      sim$CIF_all$meta
    )
  }
  rownames(cell_meta) <- paste0("cell", 1:sim$N$cell)
  
  if (is.null(sim$GRN)) {
    grn_params <- NULL
  } else {
    grn_params <- sim$GRN$params
    colnames(grn_params) <- c("target", "regulator", "effect")
    grn_params$regulator <- paste0("gene", grn_params$regulator)
    grn_params$target <- paste0("gene", grn_params$target)
  }
  
  counts <- t(sim$counts_s)
  rownames(counts) <- paste0("gene", 1:sim$N$gene)
  colnames(counts) <- paste0("cell", 1:sim$N$cell)

  result <- list(
    counts = counts,
    cif = sim$CIF_all$cif,
    giv = sim$GIV,
    cell_meta = cell_meta,
    kinetic_params = sim$params,
    atacseq_data = t(sim$atac_data),
    region_to_gene = sim$region_to_gene,
    num_genes = sim$N$gene,
    grn_params = grn_params,
    hge_scale = sim$hge_scale,
    .options = sim$options,
    .grn = sim$GRN,
    .n = sim$N
  )
  

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
        1:sim$N$cell, \(i) sim$GRN$history[[sim$dyngrn_ver_map[i]]]
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
      ctype_param$ligand <- paste0("gene", ctype_param$ligand)
      ctype_param$receptor <- paste0("gene", ctype_param$receptor)
    } else {
      ctype_param <- NULL
    }
    
    cci_locs <- do.call(rbind, sim$grid$locs)
    colnames(cci_locs) <- c("x", "y")
    rownames(cci_locs) <- paste0("cell", 1:sim$N$cell)
    
    result <- c(result, list(
      grid = sim$grid,
      cci_locs = cci_locs,
      cci_cell_type_param = ctype_param,
      cci_cell_types = sim$cell_type_map
    ))
  }

  if (is_debug) {
    result <- c(result, list(
      sim = sim
    ))
  }

  resenv <- new.env()
  attr(resenv, "name") <- "scMultiSim Result"
  
  for (n in names(result)) {
    resenv[[n]] <- result[[n]]
  }
  resenv
}


#' Rename the original gene IDs in the GRN table to integers.
#' @return list
.normalize_GRN_params <- function(params) {
  if (!is.data.frame(params)) {
    return(NULL)
  }
  
  if (nrow(unique(params[,1:2])) != nrow(params)) {
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


#' Title
#'
#' @param GRN
#' @param options
#'
#' @return A list with the following keys
#' gene
#' grn.gene
#' non.grn.gene
#' regulator
#' target
#' cif
#' diff.cif
#' nd.cif
.get_numbers <- function(GRN, options) {
  N <- list()
  N$cell <- OP(num.cells)

  # gene
  if (is.null(GRN)) {
    N$grn.gene <- 0
    N$gene <- OP(num.genes)
    N$no.grn.gene <- N$gene
  } else {
    N$grn.gene <- length(GRN$name_map)
    N$gene <- if (is.numeric(options$num.genes) && options$num.genes >= N$grn.gene) {
      options$num.genes
    } else {
      ceiling((N$grn.gene + N$grn.gene * OP(unregulated.gene.ratio)) / 10) * 10
    }
    N$non.grn.gene <- N$gene - N$grn.gene
    N$regulator <- GRN$n_reg
    N$target <- GRN$n_tgt
  }

  # cif
  N$cif <- OP(num.cifs)
  n_diff_cif <- ceiling(N$cif * OP(diff.cif.fraction))
  n_non_cif <- N$cif - n_diff_cif
  is_vary <- switch(OP(vary),
    "all"         = c(T, T, T),
    "kon"         = c(T, F, F),
    "koff"        = c(F, T, F),
    "s"           = c(F, F, T),
    "except_kon"  = c(F, T, T),
    "except_koff" = c(T, F, T),
    "except_s"    = c(T, T, F)
  )
  N$diff.cif <- sapply(is_vary, function(x) ifelse(x, n_diff_cif, 0))
  N$nd.cif <- N$cif - N$diff.cif

  # regions
  N$region <- length(OP(region.distrib)) * N$gene

  # data: param density
  data(param_realdata.zeisel.imputed)
  match_params[, 1:3] <- log(base = 10, match_params[, 1:3])
  N$params_den <- lapply(1:3, function(i) {
    density(match_params[, i], n = 2000)
  })

  # return
  N
}


.discrete_cif <- function(seed, N, options) {
  phyla <- OP(tree)
  cif_center <- OP(cif.center)
  cif_sigma <- OP(cif.sigma)
  min_popsize <- OP(discrete.min.pop.size)
  i_minpop <- OP(discrete.min.pop.index)

  npop <- length(phyla$tip.label)
  if (npop == 1) {
    ncells_pop <- N$cell
  } else {
    ncells_pop <- rep(min_popsize, npop)
    if (N$cell < min_popsize * npop) {
      stop(sprintf(
        "The size of the smallest population (%g * %g) is too big for the total number of cells (%g)",
        min_popsize, npop, N$cell))
    }

    larger_pops <- setdiff(1:npop, i_minpop)
    ncells_pop[larger_pops] <- floor((N$cell - min_popsize) / length(larger_pops))
    leftover <- N$cell - sum(ncells_pop)
    if (leftover > 0) {
      temp <- sample(larger_pops, leftover, replace = F)
      ncells_pop[temp] <- ncells_pop[temp] + 1
    }
  }

  vcv_evf_mean <- vcv.phylo(phyla, cor = T)
  param_name <- c("kon", "koff", "s")

  evfs <- lapply(1:3, function(iparam) {
    n_nd_cif <- N$nd.cif[iparam]
    n_diff_cif <- N$diff.cif[iparam]
    if (n_nd_cif > 0) {
      pop_evf_nonDE <- lapply(c(1:npop), function(ipop) {
        evf <- sapply(c(1:(n_nd_cif)), function(ievf) {
          rnorm(ncells_pop[ipop], cif_center, cif_sigma)
        })
        return(evf)
      })
      pop_evf_nonDE <- do.call(rbind, pop_evf_nonDE)
      colnames(pop_evf_nonDE) <- rep("nonDE", n_nd_cif)
    } else {
      pop_evf_nonDE <- NULL
    }
    if (n_diff_cif > 0) {
      pop_evf_mean_DE <- mvrnorm(n_diff_cif, rep(cif_center, npop), vcv_evf_mean)
      pop_evf_DE <- lapply(c(1:npop), function(ipop) {
        evf <- sapply(c(1:n_diff_cif), function(ievf) {
          rnorm(ncells_pop[ipop], pop_evf_mean_DE[ievf, ipop], cif_sigma)
        })
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
      1:(n_nd_cif + n_diff_cif)
    )
    if (iparam <= 2 && N$reg_cif > 0) {
      reg_cif <- lapply(
        1:N$reg_cif,
        \(.) rnorm(N$cell, cif_center, cif_sigma)
      ) %>% do.call(cbind, .)
      colnames(reg_cif) <- paste(param_name, "reg", 1:N$reg_cif, sep = "_")
      cif <- cbind(cif, reg_cif)
    }
    return(cif)
  })

  names(evfs) <- param_name
  meta <- data.frame(pop = do.call(c, lapply(c(1:npop), function(i) {
    rep(i, ncells_pop[i])
  })))

  list(cif = evfs, meta = meta)
}


#' Generates cifs for cells sampled along the trajectory of cell development
#'
#' @param seed
#' @param options the option list
#' @param N the number list
#' @param is_spatial return a list of cifs for spatial
#' @param .plot
#' @param .plot.name
.continuous_cif <- function(seed, N, options, ncell_key = "cell", is_spatial = F, spatial_params = NULL,
                            .plot = F, .plot.name = "cont_cif.pdf") {
  set.seed(seed)

  ncells <- N[[ncell_key]]
  phyla <- OP(tree)
  cif_center <- OP(cif.center)
  cif_sigma <- OP(cif.sigma)
  use_impulse <- OP(use.impulse)
  tree_info <- .tree_info(phyla)
  neutral <- SampleSubtree(
    tree_info$root, 0, OP(cif.center),
    tree_info$edges,
    if (is_spatial) N$max_layer else ncells,
    N$step_size,
    neutral = NA
  )


  # ==== generate cif: cell x n_cif ============================================

  cif <- .continuous_cif_param(
    is_spatial,
    ncells, N$nd.cif, N$diff.cif, N$reg_cif,
    cif_center, cif_sigma, N$step_size,
    neutral, phyla, tree_info,
    use_impulse,
    sp_params = spatial_params
  )

  # ===== metadata & output ====================================================

  if (is_spatial) {
    c(cif, list(neutral = neutral))
  } else {
    meta <- data.frame(
      pop = apply(neutral[, 1:2], 1, \(X) paste0(X, collapse = "_")),
      depth = neutral[, 3]
    )[1:ncells, ]

    list(cif = cif, meta = meta, neutral = neutral)
  }
}


#' gene x regulator
.gene_effects_by_regulator <- function(seed, GRN, N) {
  set.seed(seed)

  geff <- matrix(0L, nrow = N$grn.gene, ncol = N$regulator)
  for (r in 1:nrow(GRN$params)) {
    c(target, regulator, effect) %<-% GRN$params[r, ]
    regu_idx <- which(GRN$regulators %in% regulator)
    geff[target, regu_idx] <- effect
  }

  # add rows for non-GRN genes
  geff <- rbind(geff, matrix(0, nrow = N$non.grn.gene, ncol = N$regulator))
  colnames(geff) <- GRN$regulators
  rownames(geff) <- c(
    sapply(seq_len(N$grn.gene), \(i) GRN$name_map[[i]]),
    seq_len(N$non.grn.gene) + N$grn.gene
  )
  geff
}


#' Title
#'
#' @param size
#' @param prob
#' @param mean
#' @param sd
#'
#' @return size x cif matrix
.identity_vectors <- function(size, n_cif, prob, mean, sd) {
  lapply(1:size, function(i) {
    nonzero <- sample(c(0, 1),
      size = n_cif,
      prob = c(1 - prob, prob),
      replace = T
    )
    nonzero[nonzero != 0] <- rnorm(sum(nonzero), mean, sd)
    nonzero
  }) %>% do.call(rbind, .)
}


#' Title
#'
#' @param seed
#' @param GRN
#' @param N
#' @param options
#'
#' @return a list of kon, koff (gene x cif+regu), s (gene x cif)
.gene_identify_vectors <- function(seed, sim, options) {
  set.seed(seed)

  GRN <- sim$GRN
  N <- sim$N

  # calculate initial gene effect values with original Symsim approach
  param_names <- c("kon", "koff", "s")
  giv <- lapply(1:3, function(i) {
    .identity_vectors(N$gene, N$cif,
      prob = OP(giv.prob),
      mean = OP(giv.mean),
      sd = OP(giv.sd)
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
    non_grn_gene <- setdiff(1:N$gene, c(regu_list, tgt_list))
    giv$s[-non_grn_gene, ] <- 0
    regu_diff_cif <- matrix(0, n_reg, N$diff.cif[3])
    # regulator rows: gene effect of 2 added to two random differential cifs
    # for each master regulator gene
    indices <- replicate(n_reg, sample(1:(N$diff.cif[3]), 2, replace = F)) %>% as.vector()
    regu_diff_cif[cbind(rep(1:n_reg, each = 2), indices)] <- 2
    giv$s[regu_list, ] <- cbind(matrix(0, n_reg, N$nd.cif[3]), regu_diff_cif)
    # target rows: for every target gene, it should use the same gene_effect vector
    # as its regulators. if a gene has multiple regulators, its gene effects will
    # be the combination of that of the regulators.
    # (1 * reg) * (reg x cif)
    regu_counts <- rowSums(GRN$geff != 0)
    stopifnot(sum(regu_counts > 0) == GRN$n_tgt)
    grn_eff <- na.omit(GRN$geff %*% giv$s[GRN$regulators, ] / regu_counts / 2)
    # ---
    rg_row <- which(rowSums(giv$s[GRN$targets, ]) > 0)
    grn_eff[rg_row, ] <- grn_eff[rg_row, ] * 0.5 + giv$s[GRN$targets[rg_row], ] * 0.5
    grn_target <- giv$s[GRN$targets, ] <- grn_eff
    # ---
    if (sim$do_spatial) {
      sp_eff <- sim$sp_effect[sim$sp_targets, 1:N$sp_regulators]
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



#' Title
#'
#' @param seed
#' @param GRN
#' @param N
#' @param options
#'
#' @return region x cif matrix
.region_identity_vectors <- function(seed, GRN, N, options) {
  set.seed(seed)
  .identity_vectors(N$region, N$cif,
    prob = OP(riv.prob),
    mean = OP(riv.mean),
    sd = OP(riv.sd)
  )
}


#' Title
#'
#' @param seed
#' @param N
#' @param options
#'
#' @return region x gene matrix
.region_to_gene_matrix <- function(seed, N, options) {
  set.seed(seed)

  res <- matrix(0, N$region, N$gene)
  # gene is regulated by 0, 1, or 2 regions
  regu_by <- sample(c(0, 1, 2), size = N$gene, replace = T, prob = OP(region.distrib))
  regu_by_1 <- which(regu_by == 1)
  regu_by_2 <- which(regu_by == 2)
  # for genes regulated by 1 region, select a random region
  res[cbind(sample(1:N$region, length(regu_by_1), replace = T), regu_by_1)] <- 1
  # for genes regulated by 2 regions, select 2 consecutive random regions
  region_idx <- sample(1:(N$region - 1), length(regu_by_2), replace = T) %>% c(., . + 1)
  res[cbind(region_idx, rep(regu_by_2, 2))] <- 1

  # return
  res
}
