#' Launch the Shiny App to configure the simulation
#' @export
run_shiny <- function() {
  # appDir <- system.file("shiny-app", package = "scMultiSim")
  appDir <- "inst/shiny-app"
  shiny::runApp(appDir, port = 8888, launch.browser = T)
}

generateSpatialLoc <- function(opt) {
  phyla <- opt$tree
  step_size <- opt$step_size
  ncell <- opt$ncell
  is_discrete <- opt$is_discrete
  lr_num <- opt$lr_num
  ctype_lr <- opt$ctype_lr

  ctp <- cci_cell_type_params(phyla, lr_num, ctype_lr, step_size,
                              rand = TRUE, discrete = is_discrete)

  c(paths, max_layers) %<-% .getPaths(
    list(cell = ncell),
    list(tree = phyla)
  )
  cell_path <- sample(seq_along(paths), ncell, replace = TRUE)


  tree_info <- .tree_info(phyla)
  neutral <- SampleSubtree(
    tree_info$root, 0, 0,
    tree_info$edges,
    max_layers,
    step_size,
    neutral = NA
  )

  neutral <- neutral[1:max_layers,]
  layer_idx_by_path <- lapply(paths, function(path) {
    idx <- integer()
    for (i in 1:(length(path) - 1)) {
      a <- path[i]
      b <- path[i + 1]
      idx <- c(idx, which(neutral[, 1] == a & neutral[, 2] == b))
    }
    idx
  })

  cell_types <- character(length = nrow(neutral))
  for (i in 1:nrow(tree_info$edges)) {
    c(id, from, to, len) %<-% tree_info$edges[i,]
    n_steps <- len %/% step_size + ceiling(len %% step_size)
    pts <- which(neutral[, 1] == from & neutral[, 2] == to)
    n_pts <- length(pts)
    cell_types[pts] <- if (n_steps == 1) {
      paste(from, to, sep = "_")
    } else {
      type_id <- ceiling(1:n_pts * (n_steps / n_pts))
      paste(from, to, type_id, sep = "_")
    }
  }

  meta_by_path <- lapply(seq_along(paths), function(i_path) {
    idx <- layer_idx_by_path[[i_path]]
    n <- neutral[idx,]
    data.frame(
      pop = apply(n[, 1:2], 1, \(X) paste0(X, collapse = "_")),
      cell.type = cell_types[idx]
    )
  })

  if (!is.null(ctp$type_map)) {
    for (i in seq_along(meta_by_path)) {
      meta_by_path[[i]] <- cbind(
        meta_by_path[[i]],
        data.frame(cell.type.idx = ctp$type_map[meta_by_path[[i]]$cell.type])
      )
    }
  }

  final_ctype <- integer(length = ncell)
  for (i in seq_len(ncell)) {
    final_ctype[i] <- if (is_discrete) {
      meta[i, "cell.type.idx"]
    } else {
      path_i <- cell_path[i]
      layer <- min(ncell - i + 1, nrow(meta_by_path[[path_i]]))
      meta_by_path[[path_i]][layer, "cell.type.idx"]
    }
  }


  grid <- CreateSpatialGrid(
    ncells = ncell,
    max_nbs = opt$max_nbs,
    .grid.size = opt$grid.size,
    .same.type.prob = opt$same.type.prob,
    .method = opt$layout,
    .method.param = NULL,
    .nb.radius = 1
  )

  grid$set_final_ctypes(final_ctype)
  for (i in 1:ncell) {
    new_cell_type <- if (is_discrete) meta[i, "cell.type.idx"] else cell_path[i]
    grid$allocate(i, new_cell_type)
  }

  grid
}
