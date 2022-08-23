.ver <- "1.0.0"

scmultisim_help <- function(topic = NULL) {
  if (is.null(topic)) {
    meta_help <- "Call scmultisim_help(topic) where topic can be in {\"options\"}. Printing help for options by default.\n"
    sprintf(.split_long_string(meta_help)) %>% cat()
    topic <- "options"
  }
  
  if (topic == "options") {
    sprintf("scMultiSim  v%s\n", .ver) %>% cat()
    .print_opt()
  }
  
  if (topic == "dynamicGRN") {
    .dynamic_grn_default_params(help = T)
  }
  
  if (topic == "CCI") {
    .cci_help()
  }
}
