# String concatenation
`%+%` <- function(a, b) paste0(a, b)


.defaultArgs <- function(args = NULL, ...) {
  defaults <- list2(...)
  if (is.null(args)) {
    args <- eval(substitute(list(...), env = parent.frame()))
  }
  for (name in names(args)) {
    defaults[[name]] <- args[[name]]
  }
  defaults
}


.regionToTFMatrix <- function(GRN, region_to_gene, .all.genes = FALSE) {
  res <- matrix(0, nrow = nrow(region_to_gene), ncol = GRN$n_reg)
  # GRN$geff: gene x tf
  geff <- GRN$geff > 0
  # region_to_gene: region x gene
  # for each region
  for (i in seq_len(nrow(region_to_gene))) {
    # get genes in this region
    genes <- which(region_to_gene[i, ] > 0)
    if (length(genes) == 0) {
      next
    }
    if (.all.genes) {
      # if a TF also regulates all these genes
      tfs <- which(colSums(geff[genes, , drop = F]) == length(genes))
    } else {
      tfs <- which(colSums(geff[genes, , drop = F]) > 0)
    }
    res[i, tfs] <- 1
  }
  res
}


#' sample from smoothed density function
#' @param nsample number of samples needed
#' @param den_fun density function estimated from density() from R default
#' @param reduce.mem use alternative implementation to reduce memory usage
#' @return a vector of samples
SampleDen <- function(nsample, den_fun, reduce.mem = FALSE) {
  probs <- den_fun$y / sum(den_fun$y)
  bw <- den_fun$x[2] - den_fun$x[1]
  probs_seq = seq_along(probs)
  mins <- den_fun$x[probs_seq] - 0.5 * bw
  maxs <- den_fun$x[probs_seq] + 0.5 * bw
  
  if (reduce.mem) {
    counts <- rmultinom(n = 1, size = nsample, prob = probs)
    total_samples <- sum(counts)
    samples <- runif(total_samples) *
      rep(maxs - mins, times = counts) +
      rep(mins, times = counts)
  } else {
    bin_id <- sample(size = nsample, x = probs_seq, prob = probs, replace = TRUE)
    counts <- tabulate(bin_id, nbins = length(probs))
    total_samples <- sum(counts)
    samples <- numeric(length = total_samples)
    cum_counts <- c(0, cumsum(counts))
    for (j in 1:length(counts)) {
      if (counts[j] > 0) {
        samples[(cum_counts[j] + 1):cum_counts[j + 1]] <- 
          runif(counts[j], min = mins[j], max = maxs[j])
      }
    }
  }
  
  return(samples)
}


#' Creating an example tree with 5 tips
#' @param plotting True for plotting the tree on console, False for no plot
#' @return a R phylo object
#' @export
#' @examples
#' Phyla5()
Phyla5 <- function(plotting = FALSE) {
  phyla <- rtree(2)
  phyla <- compute.brlen(phyla, 1)
  tip <- compute.brlen(phyla, 1)
  phyla <- bind.tree(phyla, tip, 1)
  phyla <- bind.tree(phyla, tip, 2)
  phyla <- bind.tree(phyla, tip, 2)
  phyla <- compute.brlen(phyla, c(1, 1, 1, 1, 1, 0.2, 0.2, 3))
  edges <- cbind(phyla$edge, phyla$edge.length)
  edges <- cbind(seq_along(edges[, 1]), edges)
  connections <- table(c(edges[, 2], edges[, 3]))
  root <- as.numeric(names(connections)[connections == 2])
  tips <- as.numeric(names(connections)[connections == 1])
  phyla$tip.label <- as.character(tips)
  if (plotting == TRUE) {
    plot(phyla, show.tip.label = FALSE, lwd = 2)
    tiplabels(cex = 2)
    nodelabels(cex = 2)
  }
  return(phyla)
}

#' Creating an example tree with 3 tips
#' @param plotting True for plotting the tree on console, False for no plot
#' @return a R phylo object
#' @export
#' @examples
#' Phyla3()
Phyla3 <- function(plotting = FALSE) {
  # par(mfrow=c(2,2))
  phyla <- rtree(2)
  phyla <- compute.brlen(phyla, 1)
  tip <- compute.brlen(phyla, 1)
  phyla <- bind.tree(phyla, tip, 1)
  phyla <- compute.brlen(phyla, c(1, 1, 1, 2))
  edges <- cbind(phyla$edge, phyla$edge.length)
  edges <- cbind(seq_along(edges[, 1]), edges)
  connections <- table(c(edges[, 2], edges[, 3]))
  root <- as.numeric(names(connections)[connections == 2])
  tips <- as.numeric(names(connections)[connections == 1])
  phyla$tip.label <- as.character(tips)

  if (plotting == TRUE) {
    plot(phyla, show.tip.label = FALSE, lwd = 2)
    tiplabels(cex = 2)
    nodelabels(cex = 2)
  }
  return(phyla)
}

#' Creating a linear example tree
#' @param len length of the tree
#' @return a R phylo object
#' @export
#' @examples
#' Phyla1(len = 1)
Phyla1 <- function(len = 1) {
  myTree <- ape::read.tree(text='(A);')
  myTree <- compute.brlen(myTree, len)
  myTree
}


# get root, internal nodes and tips from a tree.
.tree_info <- function(tree) {
  edges <- cbind(1:nrow(tree$edge), tree$edge, tree$edge.length)
  colnames(edges) <- c("id", "from", "to", "len")
  parents <- unique(edges[, 2])
  children <- unique(edges[, 3])
  root <- setdiff(parents, children) %>% as.numeric()
  tips <- setdiff(children, parents) %>% as.numeric()
  internal <- union(parents, children) %>% as.numeric()

  list(edges = edges, root = root, tips = tips, internal = internal)
}

.print_param_summary <- function(sim) {
  cat(sprintf("intr noise: %g\n", sim$options$intrinsic.noise))
  
  N <- sim$N
  cat("======== Params Summary ========\n")
  cat(sprintf("Genes: %d (%d GRN + %d Non-GRN)\n", N$gene, N$grn.gene, N$non.grn.gene))
  cat(sprintf("CIF_%s: %d (%d nd + %d diff) + %d reg",
              c("kon", "koff", "s"), N$cif, N$nd.cif, N$diff.cif, N$reg_cif), sep = "\n")
  if (!is.null(sim$GRN)) {
    cat(sprintf("GRN: %d regulators, %d targets\n", sim$GRN$n_reg, sim$GRN$n_tgt))
  }
  if (sim$do_spatial) {
    cat(sprintf("Spatial: %d regulators\n", length(sim$sp_regulators)))
  }
  
  cat("Params:\n")
  cat("  CIF ")
  if (sim$do_spatial) {
    cat("(NA)\n")
  } else {
    .print_matrix_dim(sim$CIF_all$cif$kon, "kon", newline = FALSE)
    .print_matrix_dim(sim$CIF_all$cif$koff, "koff", newline = FALSE)
    .print_matrix_dim(sim$CIF_all$cif$s, "s")
  }
  
  cat("  GIV ")
  .print_matrix_dim(sim$GIV$kon, "kon", newline = FALSE)
  .print_matrix_dim(sim$GIV$koff, "koff", newline = FALSE)
  .print_matrix_dim(sim$GIV$s, "s")
  
  cat("  Params ")
  if (sim$do_spatial) {
    .print_matrix_dim(sim$params_spatial[[1]]$kon, "kon", newline = FALSE)
    .print_matrix_dim(sim$params_spatial[[1]]$koff, "koff")
  } else {
    .print_matrix_dim(sim$params$kon, "kon", newline = FALSE)
    .print_matrix_dim(sim$params$koff, "koff")
  }
  
  .print_matrix_dim(sim$CIF_atac, "  CIF_atac")
  .print_matrix_dim(sim$region_to_gene, "  Region2Gene")
  .print_matrix_dim(sim$atac_data, "  ATAC")
  
  cat("================================\n")
}

.print_matrix_dim <- function(mtx, name = NULL, newline = TRUE) {
  if (is.null(name)) {
    cat(sprintf("%dx%d", nrow(mtx), ncol(mtx)))
  } else {
    cat(sprintf("%s: %dx%d  ", name, nrow(mtx), ncol(mtx)))
  }
  if (newline) {
    cat("\n")
  }
}

.print_time <- function(sim) {
  cat(sprintf("Time spent: %.2f mins\n", as.numeric(Sys.time() - sim$start_time, units = "mins")))
}

.print_gene_in_grn <- function(sim) {
  rg <- sim$GRN$regulators
  tg <- sim$GRN$targets
  
  if (sim$do_spatial) {
    
  } else {
    
  }
}


#' Simulate a small example dataset with 200 cells and the 100-gene GRN
#' @param velocity whether to simulate RNA velocity
#' @return the simulation result
#' @export
#' @examples
#' sim_example_200_cells()
sim_example_200_cells <- function(velocity = FALSE) {
  data(GRN_params_100, envir = environment())
  options <- list(
    rand.seed = 0,
    GRN = GRN_params_100,
    num.cells = 200,
    num.cifs = 20,
    cif.sigma = 0.5,
    tree = Phyla3(),
    diff.cif.fraction = 0.8,
    do.velocity = velocity
  )
  sim_true_counts(options)
}


#' Simulate a small example dataset with 200 cells and the 100-gene GRN, with CCI enabled
#' @return the simulation result
#' @export
#' @examples
#' sim_example_200_cells_spatial()
sim_example_200_cells_spatial <- function() {
  data(GRN_params_100, envir = environment())
  lig_params <- data.frame(
    target    = c(101, 102),
    regulator = c(103, 104),
    effect    = c(5.2, 5.9)
  )
  options <- list2(
    rand.seed = 0,
    GRN = GRN_params_100,
    num.genes = 110,
    num.cells = 200,
    num.cifs = 50,
    tree = Phyla3(),
    intrinsic.noise = 0.5,
    cci = list(
      params = lig_params,
      max.neighbors = 4,
      cell.type.interaction = "random",
      step.size = 0.5
    )
  )
  sim_true_counts(options)
}

atac_dens_nonzero <- function(data) {
  x <- data[data > 0]
  density(x = log2(x + 1), adjust = 1, n = 999)
}