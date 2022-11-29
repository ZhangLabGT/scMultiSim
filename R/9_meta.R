.ver <- "1.0.0"

#' Show detailed documentations of scMultiSim's parameters
#'
#' @param topic Can be `options`, `dynamic.GRN`, or `cci`
#'
#' @export
#'
#' @examples
scmultisim_help <- function(topic = NULL) {
  if (is.null(topic)) {
    meta_help <- "Call scmultisim_help(topic) where topic can be \"options\" or an option name. Printing help for options by default.\n"
    sprintf(.split_long_string(meta_help)) %>% cat()
    topic <- "options"
  }
  
  if (topic == "options") {
    sprintf("scMultiSim  v%s\n", .ver) %>% cat()
    .print_opt()
    return()
  }
  
  if (topic == "dynamic.GRN") {
    .dynamic_grn_default_params(help = T)
    return()
  }
  
  if (topic == "cci") {
    .cci_help()
    return()
  }
  
  .print_opt(topic)
}


.cci_help <- function() {
  cat("
To enable simulating cell-cell interaction, the value should be a list including
the following names:

- params: (data.frame)
    The spatial effect between neighbor cells.
    It should be a data frame similar to the GRN parameter.
- step.size: (number, optional)
    If using continuous population, use this step size to further divide the
    cell types on the tree. For example, if the tree only has one branch 1 -> 2
    and the branch length is 1 while the step size is 0.34, there will be totally
    three cell types: 1_2_1, 1_2_2, 1_2_3.
- cell.type.interaction: (\"random\" or a matrix)
    The interaction level between different cell types.
    They act as factors multiplied to the ligand effect.
    Supply the string \"random\" to let scMultiSim generate these factors randomly.
    Otherwise, use cci_cell_type_params() to generate the template data structure.
- max.neighbors: (integer from 1 to 4, optional)
    Constraint the maxinum number of neighbors with CCI for each cell.
    The neighbors with CCI will be randomly sampled. 
      ") 
}
