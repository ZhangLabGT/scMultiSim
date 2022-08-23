#' String concatenation
`%+%` <- function(a, b) paste0(a, b)


default.args <- function(args = NULL, ...) {
  defaults <- list2(...)
  if (is.null(args)) {
    args <- eval(substitute(list(...), env = parent.frame()))
  }
  for (name in names(args)) {
    defaults[[name]] <- args[[name]]
  }
  defaults
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


Phyla1 <- function() {
  myTree <- ape::read.tree(text='(A);')
  myTree <- compute.brlen(myTree, 1)
  myTree
}


#' get root, internal nodes and tips from a tree.
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
    .print_matrix_dim(sim$CIF_all$cif$kon, "kon", newline = F)
    .print_matrix_dim(sim$CIF_all$cif$koff, "koff", newline = F)
    .print_matrix_dim(sim$CIF_all$cif$s, "s")
  }
  
  cat("  GIV ")
  .print_matrix_dim(sim$GIV$kon, "kon", newline = F)
  .print_matrix_dim(sim$GIV$koff, "koff", newline = F)
  .print_matrix_dim(sim$GIV$s, "s")
  
  cat("  Params ")
  if (sim$do_spatial) {
    .print_matrix_dim(sim$params_spatial[[1]]$kon, "kon", newline = F)
    .print_matrix_dim(sim$params_spatial[[1]]$koff, "koff")
  } else {
    .print_matrix_dim(sim$params$kon, "kon", newline = F)
    .print_matrix_dim(sim$params$koff, "koff")
  }
  
  .print_matrix_dim(sim$CIF_atac, "  CIF_atac")
  .print_matrix_dim(sim$region_to_gene, "  Region2Gene")
  .print_matrix_dim(sim$atac_data, "  ATAC")
  
  cat("================================\n")
}

.print_matrix_dim <- function(mtx, name = NULL, newline = T) {
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