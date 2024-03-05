# parse the spatial parameters
.parseSpatialParams <- function(params, num_genes, phyla, is_discrete) {
  max_nb <- params$max.neighbors
  if (!is.numeric(max_nb) || !(max_nb %in% seq_len(4))) max_nb <- 4
  effects <- NULL
  regulators <- list()

  if (!is.data.frame(lig_param <- params$params)) {
    stop("cci$params should be a data frame")
  }

  if (!all(lig_param[, seq_len(2)] > 0) && all(lig_param[, seq_len(2)] <= num_genes)) {
    stop("spatial genes should be integres smaller than num.genes")
  }

  spatial_list <- lig_param
  targets <- sort(unique(spatial_list[, 1]))
  regulators <- sort(unique(spatial_list[, 2]))
  stopifnot(length(intersect) > 0)

  num_regulators <- length(regulators)
  effects <- matrix(0, num_genes, num_regulators * max_nb)
  # e.g. 2 regulators => c(0, 2, 4, 8)
  rg_idx_base <- 0:(max_nb - 1) * num_regulators

  for (irow in seq(nrow(spatial_list))) {
    tg <- spatial_list[irow, 1]
    rg <- which(regulators %in% spatial_list[irow, 2])
    weight <- spatial_list[irow, 3]
    effects[tg, rg_idx_base + rg] <- weight
  }

  step_size <- if (is.numeric(params$step.size)) params$step.size else Inf

  # check cell.type.params
  num_lr <- params$cell.type.lr.pairs %||% 4:6
  ctp <- params$cell.type.interaction
  ctp0 <- cci_cell_type_params(phyla, num_regulators, num_lr, step_size,
                               rand = FALSE, discrete = is_discrete)

  cell_type_map <- NULL
  cell_type_factors <- if (is.list(ctp) && length(ctp) == 2) {
    if (any(dim(ctp0$params) != dim(ctp$params))) {
      stop("invalid cell.type.interaction. Please use cci_cell_type_params()")
    }
    cell_type_map <- ctp$type_map
    ctp$params
  } else if (is.character(ctp) && ctp == "on") {
    .cellTypeParamToMatrix(ctp0)
  } else if (is.character(ctp) && ctp == "random") {
    res <- cci_cell_type_params(phyla, num_regulators, num_lr, step_size,
                                rand = TRUE, discrete = is_discrete)
    cell_type_map <- res$type_map
    res$params
  } else {
    cat("cell.type.interaction is off, cells will have the same CCI level between all cell types.\n")
    NULL
  }

  same_type_prob <- params$same.type.prob %||% 0.8
  grid_size <- params$grid.size %||% NA
  del_lr_pair <- params$del.lr.pair %||% TRUE
  layout <- params$layout %||% "enhanced"
  nb_radius <- params$radius %||% 1

  grid_ex_params <- NULL
  if (layout %in% c("basic", "enhanced", "enhanced2", "layers")) {
    grid_method <- layout
  } else if (startsWith(layout, "islands")) {
    grid_method <- "islands"
    grid_ex_params <- substring(layout, 9) %>%
      strsplit(",") %>%
      unlist() %>%
      as.numeric()
    if (length(grid_ex_params) == 0 || any(is.na(grid_ex_params))) {
      stop("layout=island: please specify the island cell types, e.g. 'layout = island:1,2'")
    }
  } else if (layout == "enhanced+oracle" || layout == "basic+oracle") {
    grid_method <- substring(layout, 1, nchar(layout) - 7)
    grid_ex_params <- "oracle"
  } else {
    stop(sprintf("CCI grid layout '%s' is not supported.", layout))
  }

  if (nb_radius < 1) {
    stop("radius must be >= 1")
  }
  if (!(layout %in% c("layers", "islands")) && nb_radius > 1) {
    stop("radius > 1 only supports layers and islands layout")
  }
  
  if (!is.null(params$single.cell.gt) && params$single.cell.gt == TRUE) {
    sc_gt <- TRUE
    static_steps <- params$static.state.len %||% 50
  } else {
    sc_gt <- FALSE
    static_steps <- 10 
  }
  
  list(
    params = spatial_list,
    regulators = regulators,
    targets = targets,
    N_regulators = num_regulators,
    N_max_neighbours = max_nb,
    effects = effects,
    cell_type_factors = cell_type_factors,
    cell_type_map = cell_type_map,
    step_size = step_size,
    del_lr_pair = del_lr_pair,
    same_type_prob = same_type_prob,
    grid_size = grid_size,
    grid_method = grid_method,
    grid_ex_params = grid_ex_params,
    sc_gt = sc_gt,
    static_steps = static_steps,
    nb_radius = nb_radius
  )
}


# get all possible paths in the tree (from root to any tip)
.getPaths <- function(N, options) {
  ncell <- N$cell
  phyla <- OP("tree")
  c(edges, root, tips, internal) %<-% .tree_info(phyla)
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
    getPaths(root, numeric())
  }

  path_abs_len <- vapply(seq_along(paths), function(i) {
    path <- paths[[i]]
    len <- 0
    for (j in seq(length(path) - 1)) {
      parent <- path[[j]]
      child <- path[[j + 1]]
      len <- len + edges[edges[, 2] == parent & edges[, 3] == child, 4]
    }
    len
  }, numeric(1))
  total_ncell <- ceiling((ncell - 2) / max(path_abs_len) * sum(phyla$edge.length))
  if (total_ncell < ncell) {
    total_ncell <- ncell
  }

  list(
    paths = paths,
    total_ncell = total_ncell
  )
}


# get the length of each path (number of cells along the path)
.getPathLen <- function(atac_neutral, paths, N) {
  edge_len <- apply(unique(atac_neutral[, seq_len(2)]), 1, function(edge) {
    c(edge,
      sum(atac_neutral[, 1] == edge[1] & atac_neutral[, 2] == edge[2], na.rm = TRUE))
  }) %>% t()

  path_len <- vapply(paths, function(path) {
    len <- 0
    for (j in seq(length(path) - 1)) {
      parent <- path[j]
      child <- path[j + 1]
      idx <- atac_neutral[, 1] == parent & atac_neutral[, 2] == child
      len <- len + edge_len[edge_len[, 1] == parent & edge_len[, 2] == child, 3]
    }
    len
  }, numeric(1))

  path_len
}


.cellTypeParamToMatrix <- function(params) {
  states <- sort(unique(params[, 1]))
  n <- length(states)
  mtx <- matrix(0, n, n)
  rownames(mtx) <- states
  colnames(mtx) <- states
  for (i in seq(nrow(params))) {
    mtx[params[i, 1], params[i, 2]] <-
      mtx[params[i, 2], params[i, 1]] <- params[i, 3]
  }
  mtx
}


#' Generate cell-type level CCI parameters
#'
#' See the return value if you want to specify the cell-type level ground truth.
#'
#' @param tree Use the same value for `sim_true_counts()`.
#' @param total.lr Total number of LR pairs in the database. Use the same value for `sim_true_counts()`.
#' @param ctype.lr If `rand` is `TRUE`, how many LR pairs should be enabled between each cell type pair. Should be a range, e.g. 4:6.
#' @param step.size Use the same value for `sim_true_counts()`.
#' @param rand Whether fill the matrix randomly
#' @param discrete Whether the cell population is discrete. Use the same value for `sim_true_counts()`.
#'
#' @return A 3D matrix of (n_cell_type, n_cell_type, n_lr). The value at (i, j, k) is 1 if there exist CCI of LR-pair k between cell type i and cell type j.
#' @export
#' @examples cci_cell_type_params(Phyla3(), 100, 4:6, 0.5, TRUE, FALSE)
#'
cci_cell_type_params <- function(tree, total.lr, ctype.lr = 4:6, step.size = 1, rand = TRUE, discrete = FALSE) {
  .tree_info(tree) %->% c(edges, root, tips, internal)

  states <- if (discrete) {
    (tree$tip.label %||% (tips %>% as.character())) %>% sort()
  } else {
    # number of steps on each edge
    n_steps <- if (is.na(step.size)) {
      rep(1, nrow(edges))
    } else {
      edges[, 4] %/% step.size + ceiling(edges[, 4] %% step.size)
    }
    lapply(seq(nrow(edges)), \(i) {
      branch <- as.character(edges[i, 2:3])
      if (n_steps[i] == 1) {
        paste(branch[1], branch[2], sep = "_")
      } else {
        paste(branch[1], branch[2], seq(n_steps[i]), sep = "_")
      }
    }) %>% unlist()
  }

  n <- length(states)
  cell_type_map <- setNames(seq(n), states)
  res <- array(0, dim = c(n, n, total.lr))

  if (rand) {
    min_lr <- max(0, min(min(ctype.lr), total.lr - 2))
    max_lr <- min(max(ctype.lr), total.lr)

    for (i in seq(n)) {
      for (j in seq(i)) {
        if (i == j && n > 1) next
        # pick 4-6 LR pair
        n_pair <- sample(min_lr:max_lr, 1)
        pairs <- sample(seq(total.lr), n_pair)
        for (p in pairs) {
          res[j, i, p] <- res[i, j, p] <- 1
        }
      }
    }
  }

  list(
    params = res,
    type_map = cell_type_map
  )
}

# generate a clutter of cells by growing from the center
.gen_clutter <- function(n_cell, grid_size = NA, center = c(0, 0),
                         existing_loc = NULL, existing_grid = NULL) {
  .in_grid <- function(x, y) {
    if (is.null(existing_grid)) {
      if (is.na(grid_size)) return(TRUE)
      x >= 1 &&
        x <= grid_size &&
        y >= 1 &&
        y <= grid_size
    } else {
      any((existing_grid[, 1] == x) & (existing_grid[, 2] == y))
    }
  }

  nb_list <- list(c(-1, 0), c(1, 0), c(0, -1), c(0, 1))
  center <- matrix(center, nrow = 1)
  locs <- center[rep(1, n_cell),]
  if (!is.null(existing_loc)) {
    locs <- rbind(locs, existing_loc)
  }
  for (i in 2:n_cell) {
    done <- FALSE
    while (!done) {
      # randomly pick a cell
      idx <- sample(seq(i - 1), 1)
      # pick a neighbor
      for (j in sample(seq_len(4), 4)) {
        nb <- locs[idx,] + nb_list[[j]]
        if (.in_grid(nb[1], nb[2]) &&
          !any((locs[, 1] == nb[1]) & (locs[, 2] == nb[2]))) {
          locs[i,] <- nb
          done <- TRUE
          break
        }
      }
    }
  }
  locs
}


#' The class for spatial grids
#' @exportClass spatialGrid
#' @field method the method to generate the cell layout
#' @field grid_size the width and height of the grid
#' @field ncells the number of cells
#' @field grid the grid matrix
#' @field locs a list containing the locations of all cells
#' @field loc_order deprecated, don't use; the order of the locations
#' @field cell_types a map to save the cell type of each allocated cell
#' @field same_type_prob the probability of a new cell placed next to a cell with the same type
#' @field max_nbs the maximum number of neighbors for each cell
#' @field nb_map a list containing the neighbors for each cell
#' @field nb_adj adjacency matrix for neighbors
#' @field nb_radius the radius of neighbors
#' @field final_types the final cell types after the final time step
#' @field pre_allocated_pos the pre-allocated positions for each cell, if any
#' @field method_param additional parameters for the layout method
#' @return a spatialGrid object
.SpatialGrid <- setRefClass("spatialGrid", fields = c(
  "method", "grid_size", "ncells", "grid", "locs", "loc_order",
  "cell_types",
  "same_type_prob",
  "max_nbs", "nb_map", "nb_adj", "nb_radius",
  "final_types", "pre_allocated_pos", "method_param"
))

.SpatialGrid$methods(
  set_final_ctypes = function(ctypes) {
    final_types <<- ctypes
    if (method == "islands") {
      #========================================================================
      pre_allocated_pos <<- data.frame(x = rep(0, ncells), y = rep(0, ncells))
      ct_other <- setdiff(unique(final_types), method_param)
      # generate the islands first
      clutter_loc <- list(); n_islands <- length(method_param)
      ncells_island <- 0
      for (ct in method_param) {
        ncells_ct <- sum(final_types == ct)
        if (ncells_ct == 0) stop(sprintf("cell type %d is not found in cell types", ct))
        clutter <- .gen_clutter(ncells_ct)
        clutter_loc <- c(clutter_loc, list(clutter))
        ncells_island <- ncells_island + ncells_ct
      }
      # generate the outline
      grid_center <- c(round(grid_size / 2), round(grid_size / 2))
      outline <- .gen_clutter(ncells_island * 2, grid_size, grid_center)
      # put the islands in the outline
      done <- FALSE
      while (!done) {
        # sample center for the islands
        centers <- outline[sample(seq(nrow(outline)), n_islands),]
        if (n_islands == 1) {
          break
        }
        # if islands don't overlap
        done <- TRUE
        for (i in seq_len(n_islands - 1)) {
          for (j in (i + 1):n_islands) {
            loc1 <- centers[i,] + clutter_loc[[i]]
            loc2 <- centers[j,] + clutter_loc[[j]]
            if (anyDuplicated(rbind(loc1, loc2), MARGIN = 1) > 0) {
              done <- FALSE
              break
            }
          }
        }
      }
      clutter_loc2 <- lapply(seq(n_islands), function(i) {
        centers[rep(i, nrow(clutter_loc[[i]])),] + clutter_loc[[i]]
      })
      # for cells not in the islands
      clutter_loc_all <- do.call(rbind, clutter_loc2)
      ncells_island <- nrow(clutter_loc_all)
      ncells_other <- ncells - ncells_island
      other_loc <- matrix(NA, nrow = ncells_other, ncol = 2)
      all_loc <- rbind(clutter_loc_all, other_loc)
      nb_list <- list(c(-1, 0), c(1, 0), c(0, -1), c(0, 1))
      for (i in (ncells_island + 1):ncells) {
        done <- FALSE
        while (!done) {
          # randomly pick a cell
          idx <- sample(seq_len(i - 1), 1)
          # pick a neighbor
          for (j in sample(seq_len(4), 4)) {
            nb <- all_loc[idx,] + nb_list[[j]]
            if (.in_grid(nb[1], nb[2]) &&
              !any((all_loc[, 1] == nb[1]) & (all_loc[, 2] == nb[2]), na.rm = TRUE)) {
              all_loc[i,] <- nb
              done <- TRUE
              break
            }
          }
        }
      }
      # shift all_loc to make sure coordinates are positive
      all_loc <- all_loc - apply(all_loc, 2, min) + 1
      # randomize background cell types
      other_loc_ct <- lapply(
        ct_other, function(i) rep(i, sum(final_types == i))
      ) %>% unlist()
      rand_cells <- sample(seq_along(other_loc_ct), round(length(other_loc_ct) * 0.5))
      other_loc_ct[rand_cells] <- sample(
        other_loc_ct[rand_cells], length(rand_cells), replace = FALSE)
      # assign cell types corresponding to all_loc
      islands_loc_ct <- lapply(
        method_param, function(i) rep(i, sum(final_types == i))
      ) %>% unlist()
      all_loc_ct <- c(islands_loc_ct, other_loc_ct)
      # randomize all cell types
      rand_cells <- sample(seq_along(all_loc_ct), round(ncells * 0.1))
      all_loc_ct[rand_cells] <- sample(
        all_loc_ct[rand_cells], length(rand_cells), replace = FALSE)
      # assign the locations based on cell types
      pre_allocated_pos <<- all_loc
      for (i in unique(final_types)) {
        pre_allocated_pos[final_types == i] <<- all_loc[all_loc_ct == i,]
      }
      #========================================================================
    } else if (method == "layers") {
      #========================================================================
      grid_center <- c(round(grid_size / 2), round(grid_size / 2))
      all_locs <- .gen_clutter(ncells, grid_size, grid_center)
      # center is bottom-left
      left_ones <- which(all_locs[,1] == min(all_locs[,1]))
      new_center <- all_locs[left_ones[which.min(all_locs[left_ones, 2])],]
      new_locs <- .gen_clutter(ncells, grid_size, new_center, existing_grid = all_locs)
      rand_cells <- sample(seq_len(ncells), round(ncells * 0.05))
      new_locs[rand_cells,] <- new_locs[sample(rand_cells, length(rand_cells), replace = FALSE),]
      pre_allocated_pos <<- new_locs[order(final_types),]
    } else {
      return()
    }
    if (nb_radius > 1) {
      nb_adj <<- matrix(FALSE, ncells, ncells)
      for (i in seq_len(ncells)) {
        for (j in seq_len(i - 1)) {
          if (.in_radius(pre_allocated_pos[i,], pre_allocated_pos[j,])) {
            nb_adj[i, j] <<- nb_adj[j, i] <<- TRUE
          }
        }
      }
    }
  },
  find_nearby = function(icell, cell.type) {
    # get available locations
    nb_list <- list(c(-1, 0), c(1, 0), c(0, -1), c(0, 1))
    locs_all <- list()
    locs_same_type <- list()
    locs_diff_type <- list()
    for (i in seq_len(icell - 1)) {
      loc <- locs[[i]]
      ctype <- cell_types[[i]]
      for (nb in nb_list) {
        x <- loc[1] + nb[1]
        y <- loc[2] + nb[2]
        # if the location is available, append it to the list
        if (.in_grid(x, y) && is.na(grid[x, y])) {
          locs_all <- append(locs_all, list(c(x, y)))
          if (ctype == cell.type)
            locs_same_type <- append(locs_same_type, list(c(x, y)))
        }
      }
    }
    locs_all <- unique(locs_all)
    locs_same_type <- unique(locs_same_type)
    locs_diff_type <- locs_all[!(locs_all %in% locs_same_type)]
    # determine the new location
    pool <- if (is.null(same_type_prob) || length(locs_same_type) == 0) {
      # same_type_prob is not defined, or same type cell unavailable
      # sample a random location
      locs_all
    } else {
      if (runif(1) < same_type_prob || length(locs_diff_type) == 0) {
        locs_same_type
      } else {
        locs_diff_type
      }
    }
    loc <- sample(pool, 1)[[1]]
  },
  allocate = function(icell, cell.type = NULL) {
    nb_map[[icell]] <<- if (is.null(nb_adj)) {
      # if nb_adj is not defined, layout is generated at the simulation time
      # don't support radius > 1, so max possible neighbors = 4
      sample(seq_len(4), max_nbs, replace = FALSE)
    } else {
      sample(which(nb_adj[icell,]), max_nbs, replace = FALSE)
    }
    loc <- if (method == "accumulated") {
      if (icell == 1) {
        # allocate the first cell at the center
        p <- floor(grid_size / 2)
        c(p, p)
      } else {
        find_nearby(icell, cell.type)
      }
    } else if (method == "enhanced") {
      if (icell == 1) {
        # allocate the first cell at the center
        p <- floor(grid_size / 2)
        c(p, p)
      } else if (icell <= 5) {
        repeat {
          l <- floor(0.3 * grid_size)
          h <- floor(0.7 * grid_size)
          x <- sample(l:h, 1)
          y <- sample(l:h, 1)
          if (is.na(grid[x, y])) break
        }
        c(x, y)
      } else {
        find_nearby(icell, cell.type)
      }
    } else if (method == "enhanced2") {
      if (icell <= 8) {
        repeat {
          l <- floor(0.15 * grid_size)
          h <- floor(0.85 * grid_size)
          x <- sample(l:h, 1)
          y <- sample(l:h, 1)
          if (is.na(grid[x, y])) break
        }
        c(x, y)
      } else {
        find_nearby(icell, cell.type)
      }
    } else if (method == "islands" || method == "layers") {
      pre_allocated_pos[icell,]
    } else {
      loc <- loc_order[icell]
      x <- ceiling(loc / grid_size)
      y <- loc %% grid_size
      y <- ifelse(y == 0, grid_size, y)
      c(x, y)
    }
    # assign the location
    locs[[icell]] <<- loc
    cell_types[[icell]] <<- cell.type
    grid[loc[1], loc[2]] <<- icell
  },
  get_neighbours = function(icell, omit.NA = TRUE) {
    # message(sprintf("nb_radius = %d", nb_radius))
    loc <- locs[[icell]]
    nbs <- nb_map[[icell]]
    x <- loc[1]; y <- loc[2]
    res <- if (nb_radius == 1) {
      c(grid_val(x - 1, y), grid_val(x + 1, y),
        grid_val(x, y - 1), grid_val(x, y + 1))[nbs]
    } else {
      if (is.null(nb_adj)) {
        stop("radius > 1 only supports layers and islands layout")
      }
      nbs
    }
    if (omit.NA) {
      na.omit(res)
    } else {
      res
    }
  },
  grid_val = function(x, y) {
    if (.in_grid(x, y)) grid[x, y] else NA
  },
  .in_grid = function(x, y) {
    x >= 1 &&
      x <= grid_size &&
      y >= 1 &&
      y <= grid_size
  },
  .in_radius = function (p1, p2) {
    sqrt((p1[1] - p2[1])^2 + (p1[2] - p2[2])^2) <= nb_radius
  }
)

CreateSpatialGrid <- function(ncells, max_nbs, .grid.size = NA, .same.type.prob = 0.8, .method = "enhanced", .method.param = NULL, .nb.radius = 1) {
  grid_size <- if (is.na(.grid.size)) ceiling(sqrt(ncells) * 3) else .grid.size
  grid <- matrix(NA, grid_size, grid_size)
  loc_order <- sample(seq(ncells))
  grid <- .SpatialGrid$new(
    method = .method,
    grid_size = grid_size, ncells = ncells, grid = grid,
    same_type_prob = .same.type.prob,
    locs = list(), loc_order = loc_order,
    cell_types = list(),
    max_nbs = max_nbs,
    nb_map = list(),
    nb_adj = NULL,
    nb_radius = .nb.radius,
    final_types = NULL,
    pre_allocated_pos = NULL,
    method_param = .method.param
  )
  grid
}

