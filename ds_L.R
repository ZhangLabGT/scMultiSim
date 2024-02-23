devtools::load_all(".")
library(BiocParallel)
library(HDF5Array)

print_size <- function(x) {
  str <- ""
  for (k in names(x)) {
    str <- paste(str, sprintf("%s:\t%s\n", k, format(object.size(x[[k]]), units = "MiB")))
  }
  message(str)
}

BPPARAM <- BiocParallel::MulticoreParam(workers=8)
setAutoBlockSize(1e8)
setAutoBPPARAM(BPPARAM)

data(GRN_params_100)

res0 <- sim_true_counts(list2(
  num.cells = 500,
  num.genes = 20000,
  num.cifs = 20,
  GRN = GRN_params_100,
  # GRN = NA,
  speed.up = T,
  optimize.mem = T,
  cif.sigma = 0.5,
  tree = Phyla3(),
  diff.cif.fraction = 0.8,
  discrete.cif = F,
  do.velocity = F,
))

