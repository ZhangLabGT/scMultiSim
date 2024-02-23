BPPARAM <- BiocParallel::MulticoreParam(workers=8)

to_ff <- function(x, use_ff, chunkdim = NULL) {
  if (!use_ff) {
    return(x)
  }
  writeHDF5Array(x, chunkdim = chunkdim)
}

to_ff_iv <- function(x, use_ff) {
  if (!use_ff) {
    return(x)
  }
  if (!is.list(x)) stop("x must be a list")
  if (!all(c("kon", "koff", "s") %in% names(x))) stop("x must have kon, koff, and s")
  list(
    kon = as(DelayedArray::DelayedArray(x$kon), "HDF5Array"),
    koff = as(DelayedArray::DelayedArray(x$koff), "HDF5Array"),
    s = as(DelayedArray::DelayedArray(x$s), "HDF5Array")
  )
}

to_ff_cif <- function(x, use_ff) {
  if (!use_ff) {
    return(x)
  }
  if (!is.list(x)) stop("x must be a list")
  list(
    cif = to_ff_iv(x$cif, use_ff),
    meta = x$meta,
    neutral = x$neutral
  )
}

create_hdf5array <- function(nrow, ncol, chunkdim = NULL, .storage.mode = "double") {
  filename <- tempfile(fileext = ".h5")
  h5createFile(filename)
  h5createDataset(filename, "data",
                  dim = c(nrow, ncol), level = 0,
                  chunk = if (is.null(chunkdim)) c(nrow, ncol) else chunkdim,
                  storage.mode = .storage.mode)
  HDF5Array(filename, "data")
}

copy_hdf5array <- function(x) {
  filename <- tempfile(fileext = ".h5")
  h5createFile(filename)
  fid1 <- H5Fopen(x@seed@filepath, flags = "H5F_ACC_RDONLY")
  fid2 <- H5Fopen(filename)
  H5Ocopy(fid1, x@seed@name, fid2, "data")
  H5Fclose(fid1)
  H5Fclose(fid2)
  HDF5Array(filename, "data")
}

update_hdf5array <- function(x, data, mask = NULL) {
  # blocked update
  filepath <- tempfile(fileext = ".h5")
  m_grid <- defaultAutoGrid(x)
  sink <- HDF5Array::HDF5RealizationSink(dim(x), filepath = filepath, name = "data")
  id <- ipcid()
  blockApply(x, function(chunk) {
    vp <- currentViewport()
    data_chunk <- if (is.null(mask)) {
      read_block(data, vp)
    } else {
      ifelse(read_block(data, vp), read_block(data, vp), chunk)
    }
    BiocParallel::ipclock(id)
    write_block(sink, vp, data_chunk)
    BiocParallel::ipcunlock(id)
    NULL
  }, grid = m_grid, BPPARAM = BPPARAM)
  HDF5Array(filepath, "data")
}

ff_apply <- function(x, FUN, ...) {
  time_ <- Sys.time()
  message("ff_apply: ", dim(x)[1], "x", dim(x)[2])
  filepath <- tempfile(fileext = ".h5")
  m_grid <- defaultAutoGrid(x)
  message("m_grid: ", dim(m_grid)[1], "x", dim(m_grid)[2])
  sink <- HDF5Array::HDF5RealizationSink(dim(x), filepath = filepath, name = "data")
  id <- ipcid()
  blockApply(x, function(chunk) {
    vp <- currentViewport()
    res <- FUN(chunk, vp, ...)
    BiocParallel::ipclock(id)
    write_block(sink, vp, res)
    BiocParallel::ipcunlock(id)
    NULL
  }, grid = m_grid, BPPARAM = BPPARAM)
  message("end: ", format(Sys.time() - time_))
  HDF5Array(filepath, "data")
}

ff_rowApply <- function(x, FUN) {
  n_row <- nrow(x)
  n_col <- ncol(x)
  message("ff_rowApply: ", n_row, "x", n_col)
  idx_list <- get_batches(n_row, n_col, batch_size = 1e9 / 8)
  arr <- create_hdf5array(n_col, n_row, chunkdim = c(n_col, length(idx_list[[1]])))

  # id <- ipcid()
  # bplapply(
  lapply(
    idx_list,
    function(idx, FUN, input_h5, input_name, output_h5, output_name, id) {
      x_chunk <- h5read(input_h5, input_name, index = list(idx, NULL))
      res <- apply(x_chunk, 1, FUN)
      # BiocParallel::ipclock(id)
      h5write(res, output_h5, name = output_name, index = list(NULL, idx))
      # BiocParallel::ipcunlock(id)
      NULL
    },
    FUN, x@seed@filepath, x@seed@name,
    arr@seed@filepath, arr@seed@name,
    # id, BPPARAM = BPPARAM)
  )

  arr
}

create_hdf5vector <- function(n, data = NULL, chunk = NULL, .storage.mode = "double") {
  filename <- tempfile(fileext = ".h5")
  h5createFile(filename)
  if (!is.null(data)) {
    n <- length(data)
  }
  h5createDataset(filename, "data",
                  dim = n, chunk = chunk, level = 0,
                  storage.mode = .storage.mode)
  if (!is.null(data)) {
    h5write(data, filename, "data")
  }
  HDF5Array(filename, "data")
}

get_batches <- function(dim1, dim2, batch_size = 1e7 / 8) {
  chunk_size <- ceiling(batch_size / dim2)
  n_chunks <- ceiling(dim1 / chunk_size)
  inc <- n_chunks / dim1
  s <- 0
  idx <- 1
  for (i in seq_len(dim1)) {
    if (s >= 1) {
      idx <- c(idx, i)
      s <- s - 1
    }
    s <- s + inc
  }
  lapply(seq_along(idx), function (i) {
    if (i == length(idx)) {
      idx[i]:dim1
    } else {
      idx[i]:(idx[i + 1] - 1)
    }
  })
}

min_nonzero <- function(x) {
  min_value <- Inf # Initialize with Inf, so any real number is smaller

  # Function to update min value with the smallest non-zero value in the current chunk
  update_min <- function(chunk) {
    non_zero_values <- chunk[chunk != 0] # Filter out zero values
    if (length(non_zero_values) > 0) { # If there are non-zero values
      current_min <- min(non_zero_values)
      if (current_min < min_value) {
        min_value <<- current_min # Update global min_value if current_min is smaller
      }
    }
  }

  blockApply(x, update_min)
  min_value
}

update_delayed <- function(sim, obj_name, data, row_idx)  {
  if (is(sim[[obj_name]], "DelayedArray")) {
    h5write(data,
            sim[[obj_name]]@seed@filepath,
            sim[[obj_name]]@seed@name,
            index = list(row_idx, NULL))
  } else {
    sim[[obj_name]][row_idx,] <- data
  }
}


ff_matmul <- function(A, B, t = FALSE, chunkdim = NULL, serial = FALSE) {
  stopifnot(is(A, "HDF5Array") && is(B, "HDF5Array"))

  dim_A <- dim(A)
  dim_B <- dim(B)

  n_row <- dim_A[1]
  n_col <- if (t) dim_B[1] else dim_B[2]
  n_mid <- dim_A[2]
  output_file <- tempfile(fileext = ".h5")
  message("ff_matmul: ", n_row, "x", n_col, "n_mid", n_mid)

  # Define chunk sizes
  if (!is.null(chunkdim)) {
    chunk_size_A <- chunkdim[1]
    chunk_size_B <- chunkdim[2]
  } else if (n_mid > 1000) {
    chunk_size_A <- ceiling(n_row / 10)
    chunk_size_B <- ceiling(n_col / 10)
  } else {
    chunk_size_A <- ceiling(n_row / 4)
    chunk_size_B <- ceiling(n_col / 4)
  }

  message("chunk_size_A: ", chunk_size_A, " chunk_size_B: ", chunk_size_B)
  time_ <- Sys.time()
  # Create an empty HDF5Array to store the result
  C <- HDF5Array::writeHDF5Array(
    array(0, dim = c(n_row, n_col)),
    filepath = output_file,
    name = "C",
    chunkdim = c(chunk_size_A, chunk_size_B)
  )
  message("writeHDF5Array: ", format(Sys.time() - time_))

  chunks <- list()
  # Loop through chunks of A and B
  for (i in seq(1, n_row, by = chunk_size_A)) {
    for (j in seq(1, n_col, by = chunk_size_B)) {
      # Calculate end indices for the chunks
      i_end <- min(i + chunk_size_A - 1, n_row)
      j_end <- min(j + chunk_size_B - 1, n_col)
      chunks <- c(chunks, list(c(i, i_end, j, j_end)))
    }
  }

  time_ <- Sys.time()
  if (serial) {
    lapply(chunks, function(chunk, output_file, t, id) {
      .process_chunk(chunk[1], chunk[2], chunk[3], chunk[4], t,
                     a_h5=A@seed@filepath, a_name = A@seed@name,
                     b_h5=B@seed@filepath, b_name = B@seed@name,
                     result_h5=output_file, id=id)
      NULL
    }, output_file, t, NULL)
  } else {
    id <- ipcid()
    tasks <- bplapply(chunks, function(chunk, output_file, t, id) {
      .process_chunk(chunk[1], chunk[2], chunk[3], chunk[4], t,
                     a_h5=A@seed@filepath, a_name = A@seed@name,
                     b_h5=B@seed@filepath, b_name = B@seed@name,
                     result_h5=output_file, id=id)
      NULL
    }, output_file, t, id, BPPARAM = BPPARAM)
  }
  message("bplapply: ", format(Sys.time() - time_))

  C
}

.process_chunk <- function(sA, eA, sB, eB, transposed,
                           a_h5, a_name, b_h5, b_name, result_h5, id) {
  log_file <- sprintf("worker_%s.log", Sys.getpid())
  # msg <- function(...) write(do.call(paste0, list(...)) , log_file, append = TRUE)
  msg <- function(...) NULL

  time_ <- Sys.time()
  # Load chunks into memory
  chunk_A <- h5read(a_h5, a_name, index = list(sA:eA, NULL))
  b_idx <- if (transposed) list(sB:eB, NULL) else list(NULL, sB:eB)
  chunk_B <- h5read(b_h5, b_name, index = b_idx)
  msg("h5read: ", format(Sys.time() - time_))
  time_ <- Sys.time()

  # Perform matrix multiplication on the chunks
  result_chunk <- chunk_A %*% (if (transposed) t(chunk_B) else chunk_B)

  msg("matmul: ", format(Sys.time() - time_))
  time_ <- Sys.time()
  # Save the result chunk to the HDF5 file
  if (!is.null(id)) BiocParallel::ipclock(id)
  msg("wait", format(Sys.time() - time_))
  time_ <- Sys.time()
  h5write(result_chunk, result_h5, "C", index = list(sA:eA, sB:eB))
  if (!is.null(id)) BiocParallel::ipcunlock(id)
  msg("h5write: ", format(Sys.time() - time_))
}

.mem_size <- function (n) {
  n * 8 / 1024^2
}
