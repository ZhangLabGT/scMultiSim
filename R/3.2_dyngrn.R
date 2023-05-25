.dynGRN <- setRefClass("dynGRN", fields = c(
  "randseed", "randstate",
  "involved_genes",
  "geff", "params", "regulators", "targets", "n_tgt", "n_reg", "name_map",
  "del_edges", "gen_edges", "has_tf_edges",
  "max_steps", "remaining_steps", "remaining_cells",
  "cell_per_step", "n_edges", "n_changing_edges", "weight_mean", "weight_sd",
  "history", "version" 
))

.dynGRN$methods(
  restructure = function() {
    if (is.null(randstate)) {
      set.seed(randseed)
    } else {
      set.seed(randstate)
    }
    # set all deleted edges' weights to 0
    if (!is.null(del_edges)) {
      geff[del_edges[,1:2]] <<- 0
    }
    grn_region <- geff[involved_genes,]
    edges <- which(grn_region != 0, arr.ind = TRUE)
    nonedges <- which(grn_region == 0, arr.ind = TRUE)
    if (!has_tf_edges) {
      nonedges <- nonedges[-(nonedges[,1] %in% regulators),]
    }
    N_changed_edges <- if (n_changing_edges < 1) {
      as.integer(n_edges * n_changing_edges)
    } else {
      as.integer(n_changing_edges)
    }
    stopifnot(N_changed_edges > 0)
    # get new del_edges and gen_edges
    dedges <- edges[sample(nrow(edges), N_changed_edges),]
    del_edges <<- cbind(dedges, geff[dedges])
    gedges <- nonedges[sample(nrow(nonedges), N_changed_edges),]
    gen_edges <<- cbind(gedges, rnorm(N_changed_edges, mean = weight_mean, sd = weight_sd))
    stopifnot(all(geff[del_edges[,1:2]] != 0))
    stopifnot(all(geff[gen_edges[,1:2]] == 0))
    randstate <<- as.integer((sum(geff) / 1e-8) %% 1e9+7)
  }
)

.dynGRN$methods(
  update = function() {
    if (remaining_steps == 0) {
      # change grn structure
      restructure()
      remaining_steps <<- max_steps
    }
    if (remaining_cells == 0) {
      # update gradually
      s <- 1 / max_steps
      for (row in 1:nrow(del_edges)) {
        i <- del_edges[row, 1]
        j <- del_edges[row, 2]
        w <- del_edges[row, 3]
        geff[i, j] <<- geff[i, j] - w * s
        if (abs(geff[i, j]) <= 1e-5) {
          geff[i, j] <<- 0
        }
      }
      for (row in 1:nrow(gen_edges)) {
        i <- gen_edges[row, 1]
        j <- gen_edges[row, 2]
        w <- gen_edges[row, 3]
        geff[i, j] <<- geff[i, j] + w * s
      }
      remaining_steps <<- remaining_steps - 1
      remaining_cells <<- cell_per_step
    }
    remaining_cells <<- remaining_cells - 1
    # update history
    history[[version]] <<- geff
    # return version
    ver_ <- version
    version <<- version + 1
    ver_
  }
)

.CreateDynGRN <- function(grn, opts) {
  if (is.na(opts$involved.genes)) {
    opts$involved.genes <- sort(unique(c(grn$regulators, grn$targets)))
  }
  if (is.na(opts$weight.mean)) {
    opts$weight.mean <- round(mean(grn$params[,3]), digits = 2)
  }
  
  dyngrn <- .dynGRN$new(
    randseed = opts$seed,
    randstate = NULL,
    # opts
    involved_genes = opts$involved.genes,
    max_steps = opts$num.steps,
    n_changing_edges = opts$num.changing.edges,
    cell_per_step = opts$cell.per.step,
    weight_mean = opts$weight.mean,
    weight_sd = opts$weight.sd,
    has_tf_edges = opts$create.tf.edges,
    # grn
    geff = grn$geff,
    params = grn$params,
    regulators = grn$regulators,
    targets = grn$targets,
    n_tgt = grn$n_tgt,
    n_reg = grn$n_reg,
    name_map = grn$name_map,
    # other fields
    del_edges = NULL,
    gen_edges = NULL,
    remaining_steps = opts$num.steps,
    remaining_cells = opts$cell.per.step,
    n_edges = nrow(grn$params),
    history = list(),
    version = 1
  )
  dyngrn$restructure()
  return(dyngrn)
}


.getDynGRNOpts <- function(options) {
  opts <- .dynamic_grn_default_params()
  for (name in names(opts)) {
    val <- options[[name]]
    if (!is.null(val)) {
      opts[[name]] <- val
    }
  }
  opts
}
