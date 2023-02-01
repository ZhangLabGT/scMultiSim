# Parameters:
# ncells, n_nd_cif, n_diff_cif, n_reg_cif,
# cif_center, cif_sigma,
# neutral, phyla, tree_info,
# use_impulse


.continuous_cif_param <- function(is_spatial, ...) {
  if (is_spatial) {
    .continuous_cif_param.spatial(...) 
  } else {
    .continuous_cif_param.normal(...) 
  }
}


.continuous_cif_param.spatial <- function(
  ncells, N_nd.cif, N_diff.cif, n_reg_cif,
  cif_center, cif_sigma, step_size,
  neutral, phyla, tree_info,
  use_impulse, sp_params, ...
) {
  # paths: list of int vector, each path
  # cell_path: int vector, the path idx of each cell
  # path_len: int vector, the length of each path
  param_names <- c("kon", "koff", "s")
  
  sp_params %->% c(
    max_layers, paths, cell_path, path_len
  )
  
  # nd and reg cif
  cif <- foreach (i_cell = 1:ncells) %dopar% {
    i_path <- cell_path[i_cell]
    n_layers <- path_len[i_path]
    
    if (i_cell %% 100 == 0) cat(sprintf("%i..", i_cell))
    # for each cell, generate n_layer x n_cif
    cif_cell <- lapply(1:3, function(i) {
      param_name <- param_names[i]
      n_nd_cif <- N_nd.cif[i]
      n_diff_cif <- N_diff.cif[i]
      
      # nd cif
      nd_cif <- lapply(1:n_nd_cif, \(icif) rnorm(n_layers, cif_center, cif_sigma)) %>% do.call(cbind, .)
      colnames(nd_cif) <- paste(param_name, "nonDE", 1:n_nd_cif, sep = "_")
      
      # diff cif
      need_diff_cif <- n_diff_cif > 0    
      # for cell 1, output the diff_cif itself; for other cells, only output T or F
      diff_cif <- need_diff_cif
      if (need_diff_cif && i_cell == 1) {
        # diff cif is shared among all cell & layers; generate them lazily
        # make sure only generated once for kon, koff and s
        # n_layers x n_diff_cif
        # =============================================== COPY
        diff_cif <- if (use_impulse) {
          c(edges, root, tips, internal) %<-% tree_info
          # impulse model
          # pdf(file = .plot.name, width = 15, height = 5)
          tip <- rep(tips, ceiling(n_diff_cif / length(tips)))
          lapply(1:n_diff_cif, function(cif_i) {
            impulse <- Impulsecifpertip(phyla, edges, root, tips, internal, neutral, tip[cif_i], cif_sigma, cif_center, step_size)
            # if (.plot) { PlotRoot2Leave(impulse, tips, edges, root, internal) }
            re_order <- match(
              apply(neutral[, 1:3], 1, \(X) paste0(X, collapse = "_")),
              apply(impulse[, 1:3], 1, \(X) paste0(X, collapse = "_"))
            )
            return(impulse[re_order, ])
          })
          # dev.off()
        } else {
          # Gaussian sample
          lapply(1:n_diff_cif, function(icif) {
            # supply neutral to have the same t_sample values for all cells
            SampleSubtree(tree_info$root, 0, cif_center, tree_info$edges, ncells, step_size, neutral = neutral)[, 4]
          }) %>%
            do.call(cbind, .) %>%
            .[1:max_layers,]
        }
        colnames(diff_cif) <- paste(param_name, "DE", 1:n_diff_cif, sep = "_")
        # ================================================ COPY
        
        diff_cif
      }
      
      # reg cif
      reg_cif <- NULL
      if (i <= 2 && n_reg_cif > 0) {
        reg_cif <- lapply(
          1:n_reg_cif,
          \(.) rnorm(n_layers, cif_center, cif_sigma)
        ) %>% do.call(cbind, .)
        colnames(reg_cif) <- paste(param_name, "reg", 1:n_reg_cif, sep = "_")
      }
      
      # T if diff_cif is needed to be combined later
      list(nd = nd_cif, diff = diff_cif, reg = reg_cif)
    })
    
    setNames(cif_cell, param_names)
  }
  
  cat("Done\n")
  # gather diff_cif
  diff_cif_all <- list(NULL, NULL, NULL)
  for (i in 1:3) {
    d_cif <- cif[[1]][[i]]$diff
    if (!is.logical(d_cif)) {
      # if this param has diff cif, move it to diff_cif_all and replace it as F
      diff_cif_all[[i]] <- d_cif
      cif[[1]][[i]]$diff <- T
    }
  }
  
  # get the index on each path
  neutral <- neutral[1:max_layers,]
  layer_idx_by_path <- lapply(paths, function(path) {
    idx <- integer()
    for (i in 1:(length(path) - 1)) {
      a <- path[i]
      b <- path[i + 1]
      idx <- c(idx, which(neutral[,1] == a & neutral[,2] == b))
    }
    idx
  })
  
  # now process diff cif
  diff_cif_by_path <- lapply(diff_cif_all, function(d_cif) {
    lapply(seq_along(paths), function(i_path){
      if (is.null(d_cif)) return(NULL)
      d_cif[layer_idx_by_path[[i_path]],]
    })
  })
  names(diff_cif_by_path) <- param_names
  
  # cell types & meta
  cell_types <- character(length = nrow(neutral))
  for (i in 1:nrow(tree_info$edges)) {
    c(id, from, to, len) %<-% tree_info$edges[i,]
    n_steps <- len %/% step_size + ceiling(len %% step_size)
    pts <- which(neutral[,1] == from & neutral[,2] == to)
    n_pts <- length(pts)
    cell_types[pts] <- if (n_steps == 1) {
      paste(from, to, sep = "_")
    } else {
      type_id <- ceiling(1:n_pts * (n_steps / n_pts))
      paste(from, to, type_id, sep = "_")
    }
  }
  
  meta_by_path <- lapply(seq_along(paths), function(i_path){
    idx <- layer_idx_by_path[[i_path]]
    n <- neutral[idx,]
    data.frame(
      pop = apply(n[,1:2], 1, \(X) paste0(X, collapse = "_")),
      depth = n[,3],
      cell.type = cell_types[idx]
    )
  })
  
  for (d_cif in diff_cif_by_path) {
    for (i in seq_along(paths)) {
      if (is.null(d_cif[[i]])) next
      stopifnot(nrow(d_cif[[i]]) == path_len[i])
    }
  }
  
  list(
    cif = cif, diff_cif_by_path = diff_cif_by_path,
    meta_by_path = meta_by_path,
    layer_idx_by_path = layer_idx_by_path
  )
}


.continuous_cif_param.normal <- function(
  ncells, N_nd.cif, N_diff.cif, n_reg_cif,
  cif_center, cif_sigma, step_size,
  neutral, phyla, tree_info,
  use_impulse, ...
) {
  param_names <- c("kon", "koff", "s")
  
  cif <- lapply(1:3, function(i) {
    param_name <- param_names[i]
    n_nd_cif <- N_nd.cif[i]
    n_diff_cif <- N_diff.cif[i]
  
    # ========== de_cif ==========
    nd_cif <- lapply(1:n_nd_cif, \(icif) rnorm(ncells, cif_center, cif_sigma)) %>% do.call(cbind, .)
    colnames(nd_cif) <- paste(param_name, "nonDE", 1:n_nd_cif, sep = "_")
    cifs <- nd_cif
    
    # ========== nd_cif ==========
    if (n_diff_cif > 0) {
      # generate de_cif if there exist de_cifs for the parameter we are looking at
      diff_cif <- if (use_impulse) {
        c(edges, root, tips, internal) %<-% tree_info
        # impulse model
        # pdf(file = .plot.name, width = 15, height = 5)
        tip <- rep(tips, ceiling(n_diff_cif / length(tips)))
        lapply(1:n_diff_cif, function(cif_i) {
          impulse <- Impulsecifpertip(phyla, edges, root, tips, internal, neutral, tip[cif_i], cif_sigma, cif_center, step_size)
          # if (.plot) { PlotRoot2Leave(impulse, tips, edges, root, internal) }
          re_order <- match(
            apply(neutral[, 1:3], 1, \(X) paste0(X, collapse = "_")),
            apply(impulse[, 1:3], 1, \(X) paste0(X, collapse = "_"))
          )
          return(impulse[re_order, ])
        })
        # dev.off()
      } else {
        # Gaussian sample
        lapply(1:n_diff_cif, function(icif) {
          # supply neutral to have the same t_sample values for all cells
          SampleSubtree(tree_info$root, 0, cif_center, tree_info$edges, ncells, step_size, neutral = neutral)[, 4]
        }) %>%
          do.call(cbind, .) %>%
          .[1:ncells, ]
      }
      colnames(diff_cif) <- paste(param_name, "DE", 1:n_diff_cif, sep = "_")
      cifs <- cbind(nd_cif, diff_cif)
    }
    
    # ========== generate reg_cif for k_on, k_off ===========
    if (i <= 2 && n_reg_cif > 0) {
      reg_cif <- lapply(
        seq_len(n_reg_cif),
        \(.) rnorm(ncells, cif_center, cif_sigma)
      ) %>% do.call(cbind, .)
      colnames(reg_cif) <- paste(param_name, "reg", seq_len(n_reg_cif), sep = "_")
      cifs <- cbind(cifs, reg_cif)
    }
    
    return(cifs)
  })
  
  names(cif) <- param_names
  cif
}


.discrete_cif.spatial <- function(
    seed, N, options, ...
) {
  set.seed(seed)
  param_names <- c("kon", "koff", "s")
  
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
  
  # nd and reg cif
  cif <- foreach (i_cell = 1:N$cell) %dopar% {
    # === each cell ===
    n_layers <- N$cell
    
    # for each cell, generate n_layer x n_cif
    cif_cell <- lapply(1:3, function(i) {
      param_name <- param_names[i]
      n_nd_cif <- N$nd.cif[i]
      n_diff_cif <- N$diff.cif[i]
      need_diff_cif <- n_diff_cif > 0    
      
      # nd cif
      nd_cif <- lapply(1:n_nd_cif, \(icif) rnorm(n_layers, cif_center, cif_sigma)) %>% do.call(cbind, .)
      colnames(nd_cif) <- paste(param_name, "nonDE", 1:n_nd_cif, sep = "_")
      
      # reg cif
      reg_cif <- NULL
      if (i <= 2 && N$reg_cif > 0) {
        reg_cif <- lapply(
          1:N$reg_cif,
          \(.) rnorm(n_layers, cif_center, cif_sigma)
        ) %>% do.call(cbind, .)
        colnames(reg_cif) <- paste(param_name, "reg", 1:N$reg_cif, sep = "_")
      }
      
      list(nd = nd_cif, diff = need_diff_cif, reg = reg_cif)
    })
    
    setNames(cif_cell, param_names)
    # === end: each cell ===
  }
  
  
  # diff cif
  diff_cif <- lapply(1:3, function(i) {
    n_diff_cif <- N$diff.cif[i]
    need_diff_cif <- n_diff_cif > 0    
    if (need_diff_cif) {
      pop_diff_cif_mean <- mvrnorm(n_diff_cif, rep(cif_center, npop), vcv_evf_mean)
      dcif <- lapply(1:npop, function(ipop) {
        evf <- sapply(1:n_diff_cif, function(ievf) {
          rnorm(ncells_pop[ipop], pop_diff_cif_mean[ievf, ipop], cif_sigma)
        })
        return(evf)
      }) %>% do.call(rbind, .)
      colnames(dcif) <- rep("DE", n_diff_cif)
      dcif
    } else {
      NULL
    }
  })
  diff_cif <- setNames(diff_cif, param_names)
  
  pop <- do.call(c, lapply(c(1:npop), function(i) rep(i, ncells_pop[i])))
  
  meta <- data.frame(
    pop = pop, cell.type = pop, cell.type.idx = pop
  )
  
  list(
    cif = cif,
    meta = meta,
    diff_cif = diff_cif
  )
}
