library(shiny)

getGlobalObjects <- function() {
  res <- list()
  i <- 1
  obj.names <- ls(envir=.GlobalEnv)
  for (n in obj.names) {
    obj <- get(n, envir=.GlobalEnv)
    if (is.data.frame(obj) || is.matrix(obj)) {
      res[[i]] <- list(name = n, type = "data.frame")
    } else if (is.vector(obj)) {
      res[[i]] <- list(name = n, type = "vector", length = length(obj))
    } else if (is.list(obj)) {
      res[[i]] <- list(name = n, type = "list", length = length(obj))
    } else {
      res[[i]] <- list(name = n, type = "other")
    }
    i <- i + 1
  }
  res
}

getOptDefaults <- function() {
  opt_list <- .opt_list()
  names <- names(opt_list)
  opts <- seq_along(names)
  res <- list()
  
  for (i in opts) {
    n <- names[i]
    opt <- opt_list[[i]]
    if (n == "") next;
    val <- opt[[1]]
    if (!val[[1]]) {
      res[[n]] <- val[[2]]
    }
  }

  res
}


# Define server logic for random distribution app ----
server <- function(input, output, session) {
  grn_info <- reactive({
    print(input$GRN)
    error <- NULL
    grn_df <- get(input$GRN)
    if (!is.data.frame(grn_df)) {
      return(list(error = "GRN must be a data frame"))
    }
    if (ncol(grn_df) != 3) {
      return(list(error = "GRN must have 3 columns"))
    }
    rg_genes <- unique(grn_df[,2])
    tg_genes <- unique(grn_df[,1])
    all_genes <- unique(c(rg_genes, tg_genes))
    list(
      data = grn_df,
      ngenes = length(all_genes),
      nrows = nrow(grn_df),
      nregulators = length(rg_genes),
      ntargets = length(tg_genes),
      error = FALSE
    )
  })

  tree_info <- reactive({
    tree_name <- input$tree
    if (is.null(tree_name)) {
      return(list(error = "No tree selected"))
    }
    tree <- switch(tree_name,
      phyla1 = Phyla1(),
      phyla3 = Phyla3(),
      phyla5 = Phyla5(),
      tryCatch(
        { eval(parse(text=tree_name), envir=.GlobalEnv) }, error = function(e) NULL
      )
    )
    if (is.null(tree)) {
      return(list(error = paste0("Error loading tree: ", tree_name, ", please check syntax")))
    }
    edges <- cbind(1:nrow(tree$edge), tree$edge, tree$edge.length)
    colnames(edges) <- c("id", "from", "to", "len")
    parents <- unique(edges[, 2])
    children <- unique(edges[, 3])
    root <- setdiff(parents, children) %>% as.numeric()
    tips <- setdiff(children, parents) %>% as.numeric()
    internal <- union(parents, children) %>% as.numeric()

    list(error = FALSE, tree = tree, edges = edges, root = root, tips = tips, internal = internal)
  })

  grid <- eventReactive(input$submit_spatial, {
    opt <- input$submit_spatial
    g <- generateSpatialLoc(list(
      layout = opt$layout,
      tree = tree_info()$tree,
      step_size = opt$stepSize,
      ncell = opt$ncell,
      is_discrete = F,
      lr_num = 0,
      ctype_lr = 0,
      grid.size = opt$gridSize,
      same.type.prob = opt$sameTypeProb,
      max_nbs = 4
    ))

    list(
      locs = g$locs,
      size = g$grid_size,
      final_types = g$final_types
    )
  })

  # Generate a plot of the data ----
  # Also uses the inputs to build the plot label. Note that the
  # dependencies on the inputs and the data reactive expression are
  # both tracked, and all expressions are called in the sequence
  # implied by the dependency graph.
  # output$plot <- renderPlot({
  #   dist <- input$dist
  #   n <- input$n

  #   hist(d(),
  #        main = paste("r", dist, "(", n, ")", sep = ""),
  #        col = "#007bc2", border = "white")
  # })

  # Generate a summary of the data ----
  # output$summary <- renderPrint({
  #   summary(d())
  # })

  # Generate an HTML table view of the head of the data ----
  output$grn_head <- renderTable({
    info <- grn_info()
    if (is.character(info$error)) {
      NULL
    } else if (is.null(info$data)) {
      "No data available"
    } else {
      head(data.frame(info$data))
    }
  }, html.table.attributes = 'class="table table-sm"')

  output$grn_summary <- renderText({
    info <- grn_info()
    if (is.character(info$error)) {
      paste0("Error: ", info$error)
    } else if (is.null(info$data)) {
      "No data available"
    } else {
      paste("GRN with", info$nrows, "edges and" , info$ngenes, "genes, incl.", info$nregulators, "regulators and", info$ntargets, "targets")
    }
  })

  output$tree_plot <- renderPlot({
    info <- tree_info()
    if (is.null(info)) {
      NULL
    } else {
      plot(info$tree, no.margin = TRUE)
      nodelabels()
    }
  })

  observe({
    session$sendCustomMessage(type = "RObjects", getGlobalObjects())
  })
  observe({
    session$sendCustomMessage(type = "GRNInfo", grn_info())
  })
  observe({
    session$sendCustomMessage(type = "TreeInfo", tree_info())
  })
  observe({
    session$sendCustomMessage(type = "Grid", grid())
  })
  observe({
    session$sendCustomMessage(type = "Defaults", getOptDefaults())
  })

  observe({
    print(input$generatedOptions)
    if (!is.null(input$generatedOptions)) {
      eval(parse(text=input$generatedOptions), envir=.GlobalEnv)
    }
  })

  observe({
    print(input$stopApp)
    if (is.character(input$stopApp) && input$stopApp == "YES") {
      stopApp()
    }
  })
}

a <- shinyApp(ui = htmlTemplate("www/index.html"), server)
