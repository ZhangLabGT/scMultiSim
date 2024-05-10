.default <- \(...) list(FALSE, as.character(enexprs(...)), ...)
.required <- list(TRUE)

.should.be.logical <- list(
  is.logical,
  "The value should be a logical."
)

.should.be.int <- list(
  \(x) x %% 1 == 0,
  "The value should be a numeric."
)

.should.be.int.between <- function(a, b) list(
  \(x) x %% 1 == 0 && x >= a && x <= b,
  sprintf("The value should be an integer between %g and %g.", a, b)
)

.should.be.num <- list(
  is.numeric,
  "The value should be a numeric."
)

.should.be.num.between <- function(a, b) list(
  \(x) is.numeric(x) && x >= a && x <= b,
  sprintf("The value should be a numeric between %g and %g.", a, b)
)

.choose_from <- function(...) {
  opts <- list(...)
  list(
    \(x) x %in% opts,
    sprintf("The value should be one of [%s].", do.call(paste, c(opts, sep = ", ")))
  )
}


# ==============================================================================
# OPTIONS: each option should be a list(default, checker, description)
# ==============================================================================

.opt_list <- function() list(
  "GENERAL",
  rand.seed                                                              = list(
    .default(0),
    .should.be.int,
    "scMultiSim should produce the same result if all other parameters are the same."
  ),
  threads                                                                = list(
    .default(1),
    .should.be.int.between(1, 4096),
    "Set to larger than 1 to use multithreading for some part of the simulation."
  ),
  speed.up                                                               = list(
    .default(FALSE),
    .should.be.logical,
    "Use experimental speed and memory optimization."
  ),
  # ========================== Gene ============================================
  "GENE",
  GRN                                                                    = list(
    .default(NULL),
    list(
      \(x) (length(x) == 1 && is.na(x)) || (is.data.frame(x) && ncol(x) >= 3 && is.numeric(x[[3]])),
      "It should be a data frame with 3 columns (target, regulator, effect). Supply NA to disable the GRN effect."
    ),
    "The GRN network."
  ),
  grn.effect                                                             = list(
    .default(1),
    .should.be.num.between(0, Inf),
    "Overall strength of the GRN effect on the expression. Different from the effect column in the GRN data frame, which is the relative effect of each TF-target pair."
  ),
  num.genes                                                              = list(
    .default(NULL),
    .should.be.int.between(1, Inf),
    "Number of genes if GRN is disabled."
  ),
  unregulated.gene.ratio                                                 = list(
    .default(0.1),
    .should.be.num.between(0, 1),
    "Ratio of unreulated to regulated genes. Extra unregulated genes will be simulated in addition to the genes in GRN."
  ),
  giv.mean                                                               = list(
    .default(0),
    .should.be.num.between(-Inf, Inf),
    "Mean of the Gene Identity Vectors."
  ),
  giv.prob                                                               = list(
    .default(0.3),
    .should.be.num.between(0, 1),
    "Probability of non-zero values in the Gene Identity Vectors."
  ),
  giv.sd                                                                 = list(
    .default(1),
    .should.be.num.between(0, Inf),
    "Stddev of the Gene Identity Vectors."
  ),
  hge.range                                                              = list(
    .default(1),
    .should.be.num.between(1, Inf),
    "Only choose highly expressed genes after this range."
  ),
  hge.prop                                                               = list(
    .default(0),
    .should.be.num.between(0, 1),
    "Propotion of highly expressed genes."
  ),
  hge.mean                                                               = list(
    .default(5),
    .should.be.num.between(1, Inf),
    "Scale of highly expressed genes."
  ),
  hge.sd                                                                 = list(
    .default(1),
    .should.be.num.between(0, Inf),
    "Variation of highly expressed genes."
  ),
  hge.max.var                                                            = list(
    .default(500),
    .should.be.num.between(0, Inf),
    "Genes with higher variation will not be selected as highly expressed genes."
  ),
  dynamic.GRN                                                            = list(
    .default(NA),
    NULL,
    "Specification of the dynamic GRN. See scmultisim_help(\"dynamic.GRN\") for details."
  ),
  # ========================== Cell ============================================
  "CELL",
  num.cells                                                              = list(
    .default(1000),
    .should.be.int.between(0, Inf),
    "Total number of cells from all populations."
  ),
  tree                                                                   = list(
    .default(Phyla5()),
    NULL,
    "A tree defining relationship between populations."
  ),
  discrete.cif                                                           = list(
    .default(FALSE),
    .should.be.logical,
    "Whether the cell population is discrete."
  ),
  discrete.pop.size                                                      = list(
    .default(NA),
    list(
      \(x) (length(x) == 1 && is.na(x)) || all(is.integer(x)),
      "the value should be an integer vector"
    ),
    "Specify the cell numbers in each population."
  ),
  discrete.min.pop.size                                                  = list(
    .default(70),
    .should.be.int,
    "Size of the smallest discrete cell population."
  ),
  discrete.min.pop.index                                                 = list(
    .default(1),
    .should.be.int.between(0, Inf),
    "Index of the smallest discrete cell population."
  ),
  num.cifs                                                               = list(
    .default(50),
    .should.be.int,
    "Number of Cell Identity Factors for each kinetic parameter."
  ),
  diff.cif.fraction                                                      = list(
    .default(0.9),
    .should.be.num.between(0, 1),
    "Fraction of CIFs which are differential factors between cell types."
  ),
  cif.center                                                             = list(
    .default(1),
    .should.be.num,
    "Mean of the CIF values."
  ),
  cif.sigma                                                              = list(
    .default(0.1),
    .should.be.num.between(0, Inf),
    "Stddev of the CIF values."
  ),
  use.impulse                                                            = list(
    .default(FALSE),
    .should.be.logical,
    "Use the impulse model when generating the continuous CIF."
  ),
  # ========================== ATAC ============================================
  "SIMULATION - ATAC",
  atac.effect                                                            = list(
    .default(0.5),
    .should.be.num.between(0, 1),
    "The influence of chromatin accessability data on gene expression."
  ),
  region.distrib                                                         = list(
    .default(c(0.1, 0.5, 0.4)),
    list(
      \(x) x > 0 && length(x) == 3 && sum(x) == 1,
      "the value should be a vector with 3 elements sum to 1"
    ),
    "The probability that a gene is regulated by respectively 0, 1, 2 consecutive regions."
  ),
  atac.p_zero                                                            = list(
    .default(0.8),
    NULL,
    "The proportion of 0s we see in the ATAC-seq data."
  ),
  atac.density                                                           = list(
    .default(NA),
    list(
      \(x) class(x) == "density",
      "the value should be a density object."
    ),
    "Density of the non-zero ATAC-seq values. Use atac_dens_nonzero() to generate."
  ),
  riv.mean                                                               = list(
    .default(0),
    .should.be.num.between(0, Inf),
    "Mean of the Region Identity Vectors."
  ),
  riv.prob                                                               = list(
    .default(0.3),
    .should.be.num.between(0, 1),
    "Probability of non-zero values in the Region Identity Vectors."
  ),
  riv.sd                                                                 = list(
    .default(1),
    .should.be.num.between(0, Inf),
    "Stddev of the Region Identity Vectors."
  ),
  # ========================== Simulation ======================================
  "SIMULATION - RNA",
  vary                                                                   = list(
    .default("s"),
    .choose_from("all", "kon", "koff", "s", "except_kon", "except_koff", "except_s"),
    "Which kinetic parameters have differential CIFs."
  ),
  bimod                                                                  = list(
    .default(0),
    .should.be.num.between(0, 1),
    "Adjust the bimodality of gene expression, thus controlling intrinsic variation."
  ),
  scale.s                                                                = list(
    .default(1),
    NULL,
    "Scale of the s parameter. Use smaller value for cell types known to be small (like naive cells). When discrete.cif = T, it can be a vector specifying the scale.s for each cluster."
  ),
  intrinsic.noise                                                        = list(
    .default(1),
    .should.be.num.between(0, 1),
    "The weight assigned to the random sample from the Beta-Poisson distribution, where the weight of the Beta-Poisson mean value is given a weight of (1 - intrinsic.noise)."
  ),
  # ========================== Kinetic Model ===================================
  "SIMULATION - KINETIC MODEL",
  do.velocity                                                            = list(
    .default(FALSE),
    .should.be.logical,
    "Simulate using the whole kinetic model and generate RNA velocity data."
  ),
  beta                                                                   = list(
    .default(0.4),
    .should.be.num,
    "Splicing rate of each gene in the kinetic model."
  ),
  d                                                                      = list(
    .default(1),
    .should.be.num,
    "Degradation rate of each gene in the kinetic model."
  ),
  num.cycles                                                             = list(
    .default(3),
    .should.be.int.between(1, Inf),
    "For velocity mode, the number of cycles run before sampling the gene expression of a cell."),
  cycle.len                                                              = list(
    .default(1),
    .should.be.num.between(0, Inf),
    "For velocity mode, a factor multiplied by the expected time to transition from kon to koff and back to form the the length of a cycle."
  ),
  mod.cif.giv = list(
    .default(NA),
    list(
      is.function, "should be a function"
    ),
    "Modify the generated CIF and GIV. The function takes four arguments: the kinetic parameter index (1=kon, 2=koff, 3=s), the current CIF matrix, the GIV matrix, and the cell metadata dataframe. It should return a list of two elements: the modified CIF matrix and the modified GIV matrix."
  ),
  ext.cif.giv = list(
    .default(NA),
    list(
      is.function, "should be a function"
    ),
    "Add customized CIF and GIV. The function takes one argument, the kinetic parameter index (1=kon, 2=koff, 3=s). It should return a list of two elements: the extra CIF matrix (n_extra_cif x n_cells) and the GIV matrix (n_genes x n_extra_cif). Return NULL for no extra CIF and GIV."
  ),
  # ========================== Spatial =========================================
  "SIMULATION - SPATIAL",
  cci                                                                    = list(
    .default(NA),
    list(
      \(x) is.list(x) && is.data.frame(x[["params"]]),
      "Enables cell-cell interaction. See scmultisim_help(\"cci\") for details."
    ),
    "The regulation network for spatial cell-cell interaction."
  )
)

# utils: check if the option is valid
.check_opt <- function(options) {
  opt_list <- .opt_list()
  opt_list <- opt_list[!sapply(opt_list, is.character)]
  for (name in names(opt_list)) {
    c(val, checker, desc) %<-% opt_list[[name]]
    required <- val[[1]]
    user_val <- options[[name]]
    if (is.null(user_val)) {
      # if option not exist
      if (required) {
        abort(sprintf("ERROR: Option '%s' is required.\n%s", name, desc))
      } else {
        # assign default value
        options[[name]] <- val[[3]]
      }
    } else {
      # check the value    
      if (!is.null(checker)) {
        c(check, err_msg) %<-% checker
        if (!check(user_val)) {
          abort(sprintf("ERROR: Option '%s' is invalid.\n%s", name, err_msg))
        }
      }
    }
  }
  
  options
}


.split_long_string <- function(x) {
  if (!is.character(x)) return(NULL)
  ss <- strsplit(x, "(?<=.{72})", perl = TRUE)[[1]]
  do.call(paste, c(as.list(ss), sep = "\n\t"))
}


.print_opt <- function(name = NULL) {
  opt_list <- .opt_list()
  names <- names(opt_list)
  
  opts <- if (is.null(name)) {
    seq_along(names)
  } else {
    which(names %in% name)
  }
  
  if (is.null(opts) || length(opts) == 0) {
    stop(sprintf("Option %s doesn't exist.\n", name))
  }
  
  for (i in opts) {
    n <- names[i]
    opt <- opt_list[[i]]
    if (n == "") {
      sprintf("\n[%s]\n\n", opt) %>% cat()
    } else {
      c(val, checker, desc) %<-% opt
      if (val[[1]]) {
        sprintf("%s  (required)\n", n) %>% cat()
      } else {
        sprintf("%s  (default: %s)\n", n, val[[2]]) %>% cat()
      }
      sprintf("\t%s\n", .split_long_string(desc)) %>% cat()
      sprintf("\t%s\n", .split_long_string(checker[[2]])) %>% cat()
    }
  }
}


#' Get option from an object in the current environment
#'
#' @param ... the parameter name
#' @param .name get option from this object
#'
#' @return the parameter value
OP <- function(..., .name = 'options') {
  options <- get(.name, envir = caller_env())
  k <- as.character(expr(...))
  if (!(k %in% names(options))) {
    stop(sprintf("Option %s is required but not presented.", k))
  }
  options[[k]]
}


.dynamic_grn_default_params <- function(help = FALSE) {
  if (help) {
    cat("Dynamic GRN deletes and creates some edges in the GRN in each epoch.
One epoch contains multiple steps, and the change is done gradually in steps.
The specific GRN at each step will be used by one or more cells sequentially.
When an epoch is done, another epoch will start.

Available options for dynamic.GRN:
  - seed: the random seed
  - num.steps: number of steps in each epoch.
  - cell.per.step: how many cells share the GRN in the same step.
  - involved.genes: a new edge will only be created within these specified genes.
      The default value is NA, which will use all existing genes in the GRN.
  - num.changing.edges: if < 1, it means the portion of edges added/deleted in each epoch.
      if >= 1, it means the number of edges added/deleted in each epoch.
  - create.tf.edges: whether a new edge can connect two TFs in the GRN.
  - weight.mean: the mean value of the weight for a newly created edge.
      The default value is NA, meaning that it will use the mean value of the input GRN.
  - weight.sd: the standard deviation of the weight for a newly created edge.

See the returned list for the default values.
")
  }

  list(
    seed = 0,
    num.steps = 200,
    cell.per.step = 1,
    involved.genes = NA,
    num.changing.edges = 2,
    create.tf.edges = FALSE,
    weight.mean = NA,
    weight.sd = 1
  )  
}
