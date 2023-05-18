## ----"setup", include=FALSE---------------------------------------------------
require("knitr")
opts_chunk$set(fig.width=4, fig.height=3)

## ----install-packages, include=F, message=F, warning=F, eval=T----------------
(function() {
  installed <- installed.packages()[,"Package"]
  install <- function(list, fn) {
    pkg <- setdiff(list, installed)
    if (length(pkg)) fn(pkg, dependencies=TRUE)
  }
  
  r_packages <- c(
    "devtools", "dplyr", "ggplot2", "Rtsne",
    "reshape", "ape", "phytools", "repr", "KernelKnn",
    "gridExtra", "parallel", 'foreach', 'phytools', "doParallel",
    "zeallot", "gtools", "gplots", "stringi", "roxygen2", "usethis",
    "purrr"
  )
  install(r_packages, install.packages)
  
  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
    BiocManager::install(version = "3.15")
  }
  bioc_packages <- c('Biobase','IRanges','S4Vectors','SummarizedExperiment')
  install(bioc_packages, BiocManager::install)
})()

## ----load-package, quietly=TRUE, message=FALSE, warning=FALSE-----------------
library("devtools")
devtools::load_all(".")
library("scMultiSim")

## ----scmultisim-help, echo = T, results = "hide"------------------------------
scmultisim_help("options")

## ----plot-tree, fig.width = 8, fig.height = 4---------------------------------
par(mfrow=c(1,2))
Phyla5(plotting = T)
Phyla3(plotting = T)

# It's not possible to plot Phyla1() because it only contains 1 branch connecting two nodes.
Phyla1()

## ----load-grn-----------------------------------------------------------------
data(GRN_params_100)
GRN_params <- GRN_params_100
head(GRN_params)

## ----define-options-----------------------------------------------------------
options <- list(
  rand.seed = 0,
  GRN = GRN_params,
  num.cells = 1000,
  num.cifs = 50,
  cif.sigma = 0.5,
  tree = Phyla5(),
  diff.cif.fraction = 0.8,
  do.velocity = T
)

## ----run-simulation-----------------------------------------------------------
results <- sim_true_counts(options)
names(results)

## ----plot-counts, fig.width = 4, fig.height = 3.5, out.width = "60%"----------
plot_tsne(log2(results$counts + 1),
         results$cell_meta$pop,
         legend = 'pop', plot.name = 'True RNA Counts Tsne')
plot_tsne(log2(results$atacseq_data + 1),
         results$cell_meta$pop,
         legend = 'pop', plot.name = 'True ATAC-seq Tsne')

## ----plot-velocity, fig.width = 4, fig.height = 3.5, out.width = "60%"--------
plot_rna_velocity(results, arrow.length = 2)

## ----add-expr-noise, fig.width = 4, fig.height = 3.5, out.width = "60%"-------
# adds `counts_obs` to `results`
add_expr_noise(results)
# adds `counts_with_batches` to `results`
divide_batches(results, nbatch = 2)

plot_tsne(log2(results$counts_with_batches + 1),
         results$cell_meta$pop,
         legend = 'pop', plot.name = 'RNA Counts Tsne with Batches')

## ----simulate-discrete, fig.width = 4, fig.height = 3.5, out.width = "60%"----
options <- list2(
  rand.seed = 0,
  GRN = GRN_params,
  num.cells = 1000,
  num.cifs = 50,
  tree = Phyla5(),
  diff.cif.fraction = 0.8,
  discrete.cif = T
)

results <- sim_true_counts(options)

plot_tsne(log2(results$counts + 1),
         results$cell_meta$pop,
         legend = 'pop', plot.name = 'True RNA Counts Tsne')

## ----adjust-diff-cif-fraction, fig.width = 4, fig.height = 3.5, out.width = "60%"----
options <- list2 (
  rand.seed = 0,
  GRN = GRN_params,
  num.cells = 1000,
  num.cifs = 50,
  tree = Phyla5(),
  diff.cif.fraction = 0.8,
)

results <- sim_true_counts(
        options %>% purrr::list_modify(diff.cif.fraction = 0.2))
plot_tsne(log2(results$counts + 1),
         results$cell_meta$pop,
         legend = 'pop', plot.name = 'RNA Counts (diff.cif.fraction = 0.2)')

results <- sim_true_counts(
        options %>% purrr::list_modify(diff.cif.fraction = 0.8))
plot_tsne(log2(results$counts + 1),
         results$cell_meta$pop,
         legend = 'pop', plot.name = 'RNA Counts (diff.cif.fraction = 0.8)')

## ----adjust-cif-sigma, fig.width = 4, fig.height = 3.5, out.width = "60%"-----
options <- list2 (
  rand.seed = 0,
  GRN = GRN_params,
  num.cells = 1000,
  num.cifs = 50,
  tree = Phyla5(),
  diff.cif.fraction = 0.8,
  cif.sigma = 0.5
)

results <- sim_true_counts(
        options %>% purrr::list_modify(cif.sigma = 0.1))
plot_tsne(log2(results$counts + 1),
         results$cell_meta$pop,
         legend = 'pop', plot.name = 'RNA Counts (cif.sigma = 0.1)')

results <- sim_true_counts(
        options %>% purrr::list_modify(cif.sigma = 1.0))
plot_tsne(log2(results$counts + 1),
         results$cell_meta$pop,
         legend = 'pop', plot.name = 'RNA Counts (cif.sigma = 1.0)')

## ----adjust-intrinsic-noise, fig.width = 4, fig.height = 3.5, out.width = "60%"----
options <- list2 (
  rand.seed = 0,
  GRN = GRN_params,
  num.cells = 1000,
  num.cifs = 50,
  tree = Phyla5(),
  diff.cif.fraction = 0.8,
  intrinsic.noise = 1
)

results <- sim_true_counts(
        options %>% purrr::list_modify(intrinsic.noise = 0.5))
plot_tsne(log2(results$counts + 1),
         results$cell_meta$pop,
         legend = 'pop', plot.name = 'RNA Counts (intrinsic.noise = 0.5)')

results <- sim_true_counts(
        options %>% purrr::list_modify(intrinsic.noise = 1))
plot_tsne(log2(results$counts + 1),
         results$cell_meta$pop,
         legend = 'pop', plot.name = 'RNA Counts (intrinsic.noise = 1)')

## ----help-dynamic-grn---------------------------------------------------------
scmultisim_help("dynamic.GRN")

## ----define-options-dynamic-grn-----------------------------------------------
options_ <- list2(
  rand.seed = 0,
  GRN = GRN_params,
  num.cells = 1000,
  num.cifs = 50,
  tree = Phyla1(),
  diff.cif.fraction = 0.8,
  do.velocity = F,
  dynamic.GRN = list(
    cell.per.step = 3,
    num.changing.edges = 5,
    weight.mean = 0,
    weight.sd = 4
  )
)

results <- sim_true_counts(options_)

## ----show-cell-specific-grn---------------------------------------------------
# GRN for cell 1 (first 10 rows)
results$cell_specific_grn[[1]][1:10,]

## ----check-cell-specific-grn--------------------------------------------------
print(all(results$cell_specific_grn[[1]] == results$cell_specific_grn[[2]]))
print(all(results$cell_specific_grn[[2]] == results$cell_specific_grn[[3]]))
print(all(results$cell_specific_grn[[3]] == results$cell_specific_grn[[4]]))

## ----session-info-------------------------------------------------------------
sessionInfo()

