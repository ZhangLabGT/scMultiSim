.parseSpatialParams <- function(params, num_genes, phyla, is_discrete) {
  max_nb <- params$max.neighbors
  if (!is.numeric(max_nb) || !(max_nb %in% 1:4)) max_nb <- 4
  effects <- NULL
  regulators <- list()

  if (!is.data.frame(lig_param <- params$params)) {
    stop("cci$params should be a data frame")
  }

  if (!all(lig_param[, 1:2] > 0) && all(lig_param[, 1:2] <= num_genes)) {
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

  for (irow in 1:nrow(spatial_list)) {
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
                               rand = F, discrete = is_discrete)

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
                                rand = T, discrete = is_discrete)
    cell_type_map <- res$type_map
    res$params
  } else {
    cat("cell.type.interaction is off, cells will have the same CCI level between all cell types.\n")
    NULL
  }

  same_type_prob = params$same.type.prob %||% 0.8
  grid_size = params$grid.size %||% NA
  del_lr_pair = params$del.lr.pair %||% TRUE
  layout = params$layout %||% "enhanced"

  if (layout %in% c("basic", "enhanced", "enhanced2")) {
    grid_method <- layout
  } else if (layout == "enhanced+oracle" || layout == "basic+oracle") {
    grid_method <- substring(layout, 1, nchar(layout) - 7)
  } else {
    stop(sprintf("CCI grid layout '%s' is not supported.", layout))
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
    grid_method = grid_method
  )
}


.getPaths <- function(N, options) {
  ncell <- N$cell
  # get all possible paths in the tree; save them in `paths`
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

  path_abs_len <- sapply(seq_along(paths), function(i) {
    path <- paths[[i]]
    len <- 0
    for (j in 1:(length(path) - 1)) {
      parent <- path[[j]]
      child <- path[[j + 1]]
      len <- len + edges[edges[, 2] == parent & edges[, 3] == child, 4]
    }
    len
  })
  total_ncell <- ceiling((ncell - 2) / max(path_abs_len) * sum(phyla$edge.length))
  if (total_ncell < ncell) {
    total_ncell <- ncell
  }

  list(
    paths = paths,
    total_ncell = total_ncell
  )
}


.getPathLen <- function(atac_neutral, paths, N) {
  edge_len <- apply(unique(atac_neutral[, 1:2]), 1, function(edge) {
    c(edge,
      sum(atac_neutral[, 1] == edge[1] & atac_neutral[, 2] == edge[2], na.rm = T))
  }) %>% t()

  path_len <- sapply(paths, function(path) {
    len <- 0
    for (j in 1:(length(path) - 1)) {
      parent <- path[j]
      child <- path[j + 1]
      idx <- atac_neutral[, 1] == parent & atac_neutral[, 2] == child
      len <- len + edge_len[edge_len[, 1] == parent & edge_len[, 2] == child, 3]
    }
    len
  })

  path_len
}


.cellTypeParamToMatrix <- function(params) {
  states <- sort(unique(params[, 1]))
  n <- length(states)
  mtx <- matrix(0, n, n)
  rownames(mtx) <- states
  colnames(mtx) <- states
  for (i in 1:nrow(params)) {
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
#'
cci_cell_type_params <- function(tree, total.lr, ctype.lr, step.size = 1, rand = T, discrete = F) {
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
    lapply(1:nrow(edges), \(i) {
      branch <- as.character(edges[i, 2:3])
      if (n_steps[i] == 1) {
        paste(branch[1], branch[2], sep = "_")
      } else {
        paste(branch[1], branch[2], 1:n_steps[i], sep = "_")
      }
    }) %>% unlist()
  }

  n <- length(states)
  cell_type_map <- setNames(1:n, states)
  res <- array(0, dim = c(n, n, total.lr))

  if (rand) {
    min_lr <- max(0, min(min(ctype.lr), total.lr - 2))
    max_lr <- min(max(ctype.lr), total.lr)

    for (i in 1:n) {
      for (j in 1:i) {
        if (i == j && n > 1) next
        # pick 4-6 LR pair
        n_pair <- sample(min_lr:max_lr, 1)
        pairs <- sample(1:total.lr, n_pair)
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


.SpatialGrid <- setRefClass("spatialGrid", fields = c(
  "method", "grid_size", "ncells", "grid", "locs", "loc_order",
  # a map to save the cell type of each allocated cell
  "cell_types",
  # the probability of a new cell placed next to a cell with the same type
  "same_type_prob",
  "max_nbs", "nb_map"
))

.SpatialGrid$methods(
  find_nearby = function(icell, cell.type) {
    # get available locations
    nb_list <- list(c(-1, 0), c(1, 0), c(0, -1), c(0, 1))
    locs_all = list()
    locs_same_type = list()
    locs_diff_type = list()
    for (i in 1:(icell - 1)) {
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
    nb_map[[icell]] <<- sample(1:4, max_nbs)
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
  get_neighbours = function(icell, omit.NA = T) {
    loc <- locs[[icell]]
    nbs <- nb_map[[icell]]
    x <- loc[1]; y <- loc[2]
    res <- c(
      grid_val(x - 1, y),
      grid_val(x + 1, y),
      grid_val(x, y - 1),
      grid_val(x, y + 1)
    )[nbs]
    if (omit.NA) na.omit(res) else res
    # return(c(
    #   grid_val(x - 1, y)
    # ))
  },
  grid_val = function(x, y) {
    if (.in_grid(x, y)) grid[x, y] else NA
  },
  .in_grid = function(x, y) {
    x >= 1 &&
      x <= grid_size &&
      y >= 1 &&
      y <= grid_size
  }
)

CreateSpatialGrid <- function(ncells, max_nbs, .grid.size = NA, .same.type.prob = 0.8, .method = "enhanced") {
  grid_size <- if (is.na(.grid.size)) ceiling(sqrt(ncells) * 3) else .grid.size
  grid <- matrix(NA, grid_size, grid_size)
  loc_order <- sample(1:ncells)
  .SpatialGrid$new(
    method = .method,
    grid_size = grid_size, ncells = ncells, grid = grid,
    same_type_prob = .same.type.prob,
    locs = list(), loc_order = loc_order,
    cell_types = list(),
    max_nbs = max_nbs,
    nb_map = list()
  )
}

