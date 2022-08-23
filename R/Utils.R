phyla_nodes <- function(phyla) {
  edges <- cbind(phyla$edge, phyla$edge.length)
  edges <- cbind(1:length(edges[, 1]), edges)
  connections <- table(c(edges[, 2], edges[, 3]))
  if (length(edges[, 1]) == 1) {
    root <- phyla$edge[1]
    tips <- phyla$edge[2]
  } else {
    root <- as.numeric(names(connections)[connections == 2])
    tips <- as.numeric(names(connections)[connections == 1])
  }
  internal <- as.numeric(names(connections)[connections == 3])
  return(list(
    root = root, tips = tips, internal = internal, edges = edges
  ))
}

# .lst <- function(...) {
#   args <- enexprs(...)
#   res <- list2(...)
#   names_ <- names(args)
#   idx <- names_ == ""
#   names_[idx] <- as.character(args[idx])
#   names(res) <- names_
#   res
# }

.lst <- function(...) {
  args <- enexprs(...)
  env <- env(caller_env())
  res <- list()
  names_ <- names(args)
  for (i in seq_along(names_)) {
    name <- if (names_[i] == '') as.character(args[[i]]) else names_[i]
    res[[name]] <- eval(args[[i]], envir = env)
  }
  res
}