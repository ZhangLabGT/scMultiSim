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

## ----help-cci-----------------------------------------------------------------
library(scMultiSim)
scmultisim_help("cci")

## ----cci-network--------------------------------------------------------------
lig_params <- data.frame(
  target    = c(101, 102),
  regulator = c(103, 104),
  effect    = c(5.2, 5.9)
)

## -----------------------------------------------------------------------------
options_ <- list2(
  rand.seed = 0,
  GRN = GRN_params_100,
  num.genes = 200,
  num.cells = 500,
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

results <- sim_true_counts(options_)

## ----plot-cell-loc, fig.width=4.5, fig.height=4-------------------------------
plot_cell_loc(results)

## ----print-cell-loc-----------------------------------------------------------
head(results$cci_locs)

## ----session-info-------------------------------------------------------------
sessionInfo()

