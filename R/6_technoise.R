#' Add experimental noise to true counts
#'
#' @param results The scMultisim result object
#' @param ...
#' `randseed`: The random seed
#' `protocol`: `UMI` or `non-UMI`
#' `gene_len`:  A vector with lengths of all genes
#' `alpha_mean`, `alpha_sd`: rate of subsampling of transcripts during capture step
#' `depth_mean`, `depth_sd`: The sequencing depth
#'
#' @seealso
#' The underlying methods
#' \link{True2ObservedCounts} and \link{True2ObservedATAC}
#'
#' @return none
#' @export
#'
#' @examples
#' results <- sim_example(ncells = 10)
#' add_expr_noise(results)
add_expr_noise <- function(results, ...) {
  cat("Adding experimental noise...\n")
  start_time <- Sys.time()
  data(gene_len_pool, envir = environment())
  gene_len <- sample(gene_len_pool, results$num_genes, replace = FALSE)
  args <- list(...)
  if (length(args) > 0) {
    rna_args <- args[names(args)[!startsWith(names(args), "atac.")]]
    atac_args <- args[names(args)[startsWith(names(args), "atac.")]]
  } else {
    rna_args <- list(); atac_args <- list()
  }
  rna_args <- .defaultArgs(rna_args, randseed = 0,
                           protocol = "nonUMI", alpha_mean = 0.1, alpha_sd = 0.02,
                           gene_len = gene_len, depth_mean = 1e5, depth_sd = 3e3,
                           nPCR1 = 16, nPCR2 = 10)
  atac_args <- .defaultArgs(atac_args, atac.obs.prob = 0.3, atac.sd.frac = 0.5)
  rna_args$true_counts <- floor(results$counts)
  rna_args$meta_cell <- results$cell_meta
  results$counts_obs <- do.call(True2ObservedCounts, rna_args)

  atac_data <- if (!is.null(results$atac_counts)) {
    cat("Using atac_counts\n")
    results$atac_counts
  } else {
    stop()
    cat("Using atacseq_data\n")
    results$atacseq_data
  }
  results$atacseq_obs <- True2ObservedATAC(atac_data, randseed = args$randseed,
                                           observation_prob = atac_args$atac.obs.prob,
                                           sd_frac = atac_args$atac.sd.frac)
  message(sprintf("Time spent: %.2f mins\n",
                  as.numeric(Sys.time() - start_time, units = "mins")))
}


#' Divide batches for observed counts
#'
#' @param results The scMultisim result object, after running `addExprNoise()`
#' @param nbatch Number of batches
#' @param effect Batch effect size, default is 3
#' @param randseed Random seed
#'
#' @return none
#' @export
#'
#' @examples
#' results <- sim_example(ncells = 10)
#' add_expr_noise(results)
#' divide_batches(results)
divide_batches <- function(results, nbatch = 2, effect = 3, randseed = 0) {
  cat("Adding batch effects...\n")
  obs <- results$counts_obs
  if (is.list(obs)) {
    obs <- obs$counts
  }
  ngene <- nrow(obs)
  merged <- rbind(obs, results$atacseq_obs)
  if ("batch" %in% names(results$cell_meta)) {
    results$cell_meta <- results$cell_meta[, !(names(results$cell_meta) %in% "batch")]
  }
  b <- .divideBatchesImpl(
    counts = merged, meta_cell = results$cell_meta,
    nbatch = nbatch, batch_effect_size = effect, randseed = 0
  )
  results$counts_with_batches <- b$counts[seq(ngene),]
  results$atac_with_batches <- b$counts[-(seq(ngene)),]
  results$cell_meta <- b$cell_meta
}


#' Divide the observed counts into multiple batches by adding batch effect to each batch
#' @param counts gene cell matrix
#' @param meta_cell the meta information related to cells, will be combined with technical cell level information and returned
#' @param nbatch number of batches
#' @param batch_effect_size amount of batch effects. Larger values result in bigger differences between batches. Default is 1.
#' @param randseed random seed
#' @keywords internal
#' @return a list with two elements: counts and meta_cell
.divideBatchesImpl <- function(counts, meta_cell, nbatch, batch_effect_size = 1, randseed = 0) {
  # set.seed(randseed)
  ## add batch effects to observed counts
  # use different mean and same sd to create the multiplicative factor for different part (gene/region) in different batch
  ncells <- dim(counts)[2]; nparts <- dim(counts)[1]
  batchIDs <- sample(seq(nbatch), ncells, replace = TRUE)
  meta_cell2 <- data.frame(batch = batchIDs, stringsAsFactors = FALSE)
  meta_cell <- cbind(meta_cell, meta_cell2)

  mean_matrix <- matrix(0, nparts, nbatch)
  part_mean <- rnorm(nparts, 0, 1)
  temp <- lapply(seq(nparts), function(ipart) {
    return(runif(nbatch, min = part_mean[ipart] - batch_effect_size, max = part_mean[ipart] + batch_effect_size))
  })
  mean_matrix <- do.call(rbind, temp)

  batch_factor <- matrix(0, nparts, ncells)
  for (ipart in seq(nparts)) {
    for (icell in seq(ncells)) {
      batch_factor[ipart, icell] <- rnorm(n = 1, mean = mean_matrix[ipart, batchIDs[icell]], sd = 0.01)
    }
  }
  counts <- round(2^(log2(counts) + batch_factor))
  return(list(counts = counts, cell_meta = meta_cell))
}


#' This function simulates the amplification, library prep, and the sequencing processes.
#' @param true_counts_1cell the true transcript counts for one cell (one vector)
#' @param protocol a string, can be "nonUMI" or "UMI"
#' @param rate_2cap the capture efficiency for this cell
#' @param gene_len gene lengths for the genes/transcripts, sampled from real human transcript length
#' @param amp_bias amplification bias for each gene, a vector of length ngenes
#' @param rate_2PCR PCR efficiency, usually very high
#' @param nPCR1 the number of PCR cycles
#' @param nPCR2 the number of PCR cycles
#' @param LinearAmp if linear amplification is used for pre-amplification step, default is FALSE
#' @param LinearAmp_coef the coeficient of linear amplification, that is, how many times each molecule is amplified by
#' @param N_molecules_SEQ number of molecules sent for sequencing; sequencing depth
#' @keywords internal
#' @return read counts (if protocol="nonUMI") or UMI counts (if protocol="UMI)
.amplifyOneCell <- function(true_counts_1cell, protocol, rate_2cap, gene_len, amp_bias,
                            rate_2PCR, nPCR1, nPCR2, LinearAmp, LinearAmp_coef, N_molecules_SEQ) {
  ngenes <- length(gene_len)
  if (protocol == "nonUMI") {
    if (!exists("len2nfrag")) data(len2nfrag)
  } else if (protocol == "UMI") { } else  {
    stop("protocol input should be nonUMI or UMI")
  }
  inds <- vector("list", 2)
  # expand the original vector and apply capture efficiency
  # maintain a transcript index vector: which transcript the molecule belongs to
  expanded_res <- .expandToBinary(c(true_counts_1cell, 1))
  expanded_vec <- expanded_res[[1]]
  trans_idx <- expanded_res[[2]]
  inds[[1]] <- expanded_vec > 0
  expanded_vec <- expanded_vec[inds[[1]]]
  trans_idx <- trans_idx[inds[[1]]]

  rate_2cap_gene <- rate_2cap[trans_idx]
  captured_vec <- expanded_vec
  captured_vec[runif(length(captured_vec)) > rate_2cap_gene] <- 0
  if (sum(captured_vec[seq(length(captured_vec) - 1)]) < 1) { return(rep(0, ngenes)) }
  captured_vec[length(captured_vec)] <- 1

  inds[[2]] <- captured_vec > 0
  captured_vec <- captured_vec[inds[[2]]]
  trans_idx <- trans_idx[inds[[2]]]
  amp_rate <- c((rate_2PCR + amp_bias[trans_idx[seq(length(trans_idx) - 1)]]), 1)
  # pre-amplification:
  if (LinearAmp) {
    PCRed_vec <- captured_vec * LinearAmp_coef
  } else {
    temp <- runif(length(captured_vec)) < amp_rate
    temp <- temp * 2 + captured_vec - temp
    for (iPCR in 2:nPCR1) {
      eff <- runif(length(temp)) * amp_rate
      v1 <- temp * (1 - eff)
      round_down <- (v1 - floor(v1)) < runif(length(v1))
      v1[round_down] <- floor(v1[round_down])
      v1[!round_down] <- ceiling(v1[!round_down])
      temp <- v1 + 2 * (temp - v1)
    }
    PCRed_vec <- temp
  }

  if (protocol == "nonUMI") { # add fragmentation step here
    temp_vec <- PCRed_vec
    for (i in seq(2, 1, -1)) {
      temp_vec1 <- numeric(); temp_vec1[inds[[i]]] <- temp_vec;
      temp_vec <- temp_vec1; temp_vec[is.na(temp_vec)] <- 0
    }
    recovered_vec <- temp_vec[seq(length(temp_vec) - 1)]
    amp_mol_count <- numeric(ngenes);
    GI <- c(0, cumsum(true_counts_1cell));
    for (i in which(true_counts_1cell > 0)) {
      x <- recovered_vec[(GI[i] + 1):GI[i + 1]]
      amp_mol_count[i] <- sum(x)
    }

    # for every copy of each transcript, convert it into number of fragments
    frag_vec <- numeric(ngenes)
    for (igene in which(amp_mol_count > 0)) {
      frag_vec[igene] <- sum(sample(len2nfrag[as.character(gene_len[igene]),],
                                    amp_mol_count[igene], replace = TRUE)) }
    # another 8 rounds of amplification to the fragments (fragmentation bias gets amplified)
    for (iPCR in seq_len(2)) {
      frag_vec <- frag_vec + vapply(
        frag_vec,\(x) rbinom(n = 1, x, prob = rate_2PCR), numeric(1))
    }
    for (iPCR in 3:nPCR2) {
      frag_vec <- frag_vec + round(frag_vec * rate_2PCR)
    }
    SEQ_efficiency <- N_molecules_SEQ / sum(frag_vec)
    if (SEQ_efficiency >= 1) {
      read_count <- frag_vec
    } else {
      read_count <- vapply(
        frag_vec,
        \(Y) rbinom(n = 1, size = Y, prob = SEQ_efficiency), numeric(1))
    }
    return(read_count)
  } else if (protocol == "UMI") {

    prob_vec <- vapply(
      gene_len[trans_idx[seq(length(trans_idx) - 1)]], .getProb, numeric(1))
    # fragmentation:
    frag_vec <- vapply(
      seq(length(PCRed_vec) - 1),
      \(igene) rbinom(n = 1, size = PCRed_vec[igene], prob = prob_vec[igene]),
      numeric(1))

    # another 10 rounds of amplification to the fragments (fragmentation bias gets amplified)
    for (iPCR in seq_len(2)) {
      frag_vec <- frag_vec + vapply(
        frag_vec, \(x) rbinom(n = 1, x, prob = rate_2PCR), numeric(1))
    }

    frag_vec <- round(frag_vec * (1 + rate_2PCR)^(nPCR2 - 1))

    SEQ_efficiency <- N_molecules_SEQ / sum(frag_vec)
    if (SEQ_efficiency >= 1) {
      sequenced_vec <- frag_vec
    } else {
      sequenced_vec <- vapply(
        frag_vec, \(Y) rbinom(n = 1, size = Y, prob = SEQ_efficiency),
        numeric(1))
    }

    temp_vec <- c(sequenced_vec, 1)
    for (i in seq(2, 1, -1)) {
      temp_vec1 <- numeric(); temp_vec1[inds[[i]]] <- temp_vec;
      temp_vec <- temp_vec1; temp_vec[is.na(temp_vec)] <- 0
    }
    recovered_vec <- temp_vec[seq(length(temp_vec) - 1)]

    UMI_counts <- numeric(ngenes);
    GI <- c(0, cumsum(true_counts_1cell));
    for (i in which(true_counts_1cell > 0)) {
      x <- recovered_vec[(GI[i] + 1):GI[i + 1]];
      UMI_counts[i] <- sum(x > 0);
    }

    return(list(UMI_counts, sequenced_vec, sum(frag_vec > 0)))
  }
}


#' Simulate observed count matrix given technical biases and the true counts
#' @param true_counts gene cell matrix
#' @param meta_cell the meta information related to cells, will be combined with technical cell level information and returned
#' @param protocol a string, can be "nonUMI" or "UMI"
#' @param alpha_mean the mean of rate of subsampling of transcripts during capture step, default at 10 percent efficiency
#' @param alpha_sd the std of rate of subsampling of transcripts
#' @param alpha_gene_mean the per-gene scale factor of the alpha parameter, default at 1
#' @param alpha_gene_sd the standard deviation of the per-gene scale factor of the alpha parameter, default at 0
#' @param lenslope amount of length bias
#' @param nbins number of bins for gene length
#' @param gene_len a vector with lengths of all genes
#' @param amp_bias_limit range of amplification bias for each gene, a vector of length ngenes
#' @param rate_2PCR PCR efficiency, usually very high, default is 0.8
#' @param nPCR1 the number of PCR cycles in "pre-amplification" step, default is 16
#' @param nPCR2 the number of PCR cycles used after fragmentation.
#' @param LinearAmp if linear amplification is used for pre-amplification step, default is FALSE
#' @param LinearAmp_coef the coeficient of linear amplification, that is, how many times each molecule is amplified by
#' @param depth_mean mean of sequencing depth
#' @param depth_sd std of sequencing depth
#' @param randseed (should produce same result if nregions, nevf and randseed are all the same)
#' @return if UMI, a list with two elements, the first is the observed count matrix, the second is the metadata; if nonUMI, a matrix
#' @export
#' @examples
#' \donttest{
#' results <- sim_example(ncells = 10)
#' data(gene_len_pool)
#' gene_len <- sample(gene_len_pool, results$num_genes, replace = FALSE)
#' True2ObservedCounts(
#'   results$counts, results$cell_meta, protocol = "nonUMI", randseed = 1,
#'   alpha_mean = 0.1, alpha_sd = 0.05, gene_len = gene_len, depth_mean = 1e5, depth_sd = 3e3
#' )
#' }
True2ObservedCounts <- function(true_counts, meta_cell, protocol, randseed, alpha_mean = 0.1, alpha_sd = 0.002,
                                alpha_gene_mean = 1, alpha_gene_sd = 0,
                                gene_len, depth_mean, depth_sd, lenslope = 0.02, nbins = 20,
                                amp_bias_limit = c(-0.2, 0.2),
                                rate_2PCR = 0.8, nPCR1 = 16, nPCR2 = 10, LinearAmp = FALSE, LinearAmp_coef = 2000) {
  # set.seed(randseed)
  ngenes <- dim(true_counts)[1]; ncells <- dim(true_counts)[2]
  amp_bias <- .calAmpBias(lenslope, nbins, gene_len, amp_bias_limit)
  rate_2cap_lb <- 0.0005; depth_lb <- 200 # lower bound for capture efficiency and sequencing depth
  rate_2cap_vec <- .rnormTrunc(n = ncells, mean = alpha_mean, sd = alpha_sd, a = rate_2cap_lb, b = 1)
  rate_2cap_vec_gene <- .rnormTrunc(n = ngenes, mean = alpha_gene_mean, sd = alpha_gene_sd, a = 0, b = 3)
  rate_2cap <- rate_2cap_vec_gene %o% rate_2cap_vec
  depth_vec <- .rnormTrunc(n = ncells, mean = depth_mean, sd = depth_sd, a = depth_lb, b = Inf)
  observed_counts <- lapply(seq(ncells), function(icell) {
    if (icell %% 50 == 0) cat(sprintf("%d..", icell))
    .amplifyOneCell(true_counts_1cell = true_counts[, icell], protocol = protocol,
                    rate_2cap = c(rate_2cap[, icell], rate_2cap_vec[icell]),
                    gene_len = gene_len, amp_bias = amp_bias,
                    rate_2PCR = rate_2PCR, nPCR1 = nPCR1, nPCR2 = nPCR2, LinearAmp = LinearAmp,
                    LinearAmp_coef = LinearAmp_coef, N_molecules_SEQ = depth_vec[icell])
  })
  gc()

  meta_cell2 <- data.frame(alpha = rate_2cap_vec, depth = depth_vec, stringsAsFactors = FALSE)
  meta_cell <- cbind(meta_cell, meta_cell2)

  if (protocol == "UMI") {
    UMI_counts <- do.call(cbind, lapply(observed_counts, "[[", 1))
    nreads_perUMI <- lapply(observed_counts, "[[", 2)
    nUMI2seq <- vapply(observed_counts, "[[", numeric(1), 3)
    observed_counts <- UMI_counts
  } else
    observed_counts <- do.call(cbind, observed_counts)

  if (protocol == "UMI") { return(list(counts = observed_counts, cell_meta = meta_cell, nreads_perUMI = nreads_perUMI,
                                       nUMI2seq = nUMI2seq))
  } else
    return(observed_counts)
}


#' Simulate observed ATAC-seq matrix given technical noise and the true counts
#' @param atacseq_data true ATAC-seq data
#' @param observation_prob for each integer count of a particular region for a particular cell, the probability the count will be observed
#' @param sd_frac the fraction of ATAC-seq data value used as the standard deviation of added normally distrubted noise
#' @param randseed (should produce same result if nregions, nevf and randseed are all the same)
#' @return a matrix of observed ATAC-seq data
#' @export
#' @examples
#' results <- sim_example(ncells = 10)
#' True2ObservedATAC(results$atac_counts, randseed = 1)
True2ObservedATAC <- function(atacseq_data, randseed, observation_prob = 0.3, sd_frac = 0.1) {
  # set.seed(randseed)
  atacseq_data <- round(atacseq_data)
  atacseq_noisy <- atacseq_data
  for (icell in seq(ncol(atacseq_data))) {
    for (iregion in seq(nrow(atacseq_data))) {
      if (atacseq_data[iregion, icell] > 0) {
        atacseq_noisy[iregion, icell] <- rbinom(n = 1, size = atacseq_data[iregion, icell], prob = observation_prob)
        atacseq_noisy[iregion, icell] <- max(atacseq_noisy[iregion, icell] + rnorm(1, mean = 0, sd = atacseq_noisy[iregion, icell] * sd_frac), 0)
      }
    }
  }
  return(atacseq_noisy)
}


#' Simulate technical biases
#' @param lenslope amount of length bias. This value sould be less than 2*amp_bias_limit\[2\]/(nbins-1)
#' @param nbins number of bins for gene length
#' @param gene_len transcript length of each gene
#' @param amp_bias_limit range of amplification bias for each gene, a vector of length ngenes
#' @keywords internal
#' @return a vector
.calAmpBias <- function(lenslope, nbins, gene_len, amp_bias_limit) {
  ngenes <- length(gene_len)
  len_bias_bin <- (-(seq(nbins))) * lenslope
  len_bias_bin <- len_bias_bin - median(len_bias_bin)
  if (max(len_bias_bin) > amp_bias_limit[2]) {
    stop("The lenslope parameter is too large.")
  }
  max_rand_bias <- amp_bias_limit[2] - max(len_bias_bin)

  rand_bias <- rnorm(ngenes, mean = 0, sd = max_rand_bias)
  rand_bias[rand_bias > max_rand_bias] <- max_rand_bias
  rand_bias[rand_bias < -max_rand_bias] <- -max_rand_bias
  #rand_bias <- runif(ngenes, -max_rand_bias,  max_rand_bias)

  binsize <- floor(ngenes / nbins)
  genes_in_bins <- vector("list", nbins)
  bin4genes <- numeric(ngenes)
  for (ibin in seq(nbins - 1)) {
    genes_in_bins[[ibin]] <- order(gene_len)[((ibin - 1) * binsize + 1):(ibin * binsize)]
    bin4genes[genes_in_bins[[ibin]]] <- ibin
  }
  genes_in_bins[[nbins]] <- order(gene_len)[((nbins - 1) * binsize + 1):ngenes]
  bin4genes[genes_in_bins[[nbins]]] <- nbins

  len_bias <- numeric(ngenes); len_bias <- len_bias_bin[bin4genes]
  amp_bias <- rand_bias + len_bias
  return(amp_bias)
}


#' expand transcript counts to a vector of binaries of the same length of as the number of transcripts
#' @param true_counts_1cell number of transcript in one cell
#' @keywords internal
#' @return a list of two vectors, the first vector is a vector of 1s, the second vector is the index of transcripts
.expandToBinary <- function(true_counts_1cell) {
  names(true_counts_1cell) <- NULL
  expanded_vec <- rep(1, sum(true_counts_1cell))
  trans_idx <- lapply(which(true_counts_1cell > 0),
                      function(igene) rep(igene, true_counts_1cell[igene]))
  trans_idx <- unlist(trans_idx)
  return(list(expanded_vec, trans_idx))
}

#' sample from truncated normal distribution
#' @param n number of values to create
#' @param a the minimum value allowed
#' @param b the maximum value allowed
#' @param mean mean of the normal distribution
#' @param sd standard deviation of the normal distribution
#' @keywords internal
#' @return a vector of length n
.rnormTrunc <- function(n, mean, sd, a, b) {
  vec1 <- rnorm(n, mean = mean, sd = sd)
  beyond_idx <- which(vec1 < a | vec1 > b)
  if (length(beyond_idx) > 0) { # for each value < rate_2cap_lb
    substi_vec <- vapply(seq_along(beyond_idx), function(i) {
      while (TRUE) {
        temp <- rnorm(1, mean = mean, sd = sd)
        if (temp > a | temp > b) { break } }
      return(temp)
    }, numeric(1))
    vec1[beyond_idx] <- substi_vec
  }
  return(vec1)
}

.getProb <- function(glength) {
  if (glength >= 1000) { prob <- 0.7 } else {
    if (glength >= 100 & glength < 1000) { prob <- 0.78 }
    else if (glength < 100) { prob <- 0 }
  }
  return(prob)
}

#' Add outliers to the observed counts
#' @param res The scMultisim result object
#' @param prob The probability of adding outliers for each gene
#' @param factor The factor of the outliers
#' @param sd The standard deviation of the outliers
#' @param cell.num For a gene, the number of cells chosen to add outliers
#' @param max.var The maximum variance allowed
#' @export
#' @return none
add_outliers <- function (
  res, prob = 0.01, factor = 2, sd = 0.5, cell.num = 1, max.var = Inf
) {
  if (is.null(res$counts_obs)) {
    stop("No counts found in the result object")
  }
  ngenes <- nrow(res$counts_obs)
  ncells <- ncol(res$counts_obs)
  gene_range <- if (is.null(res$.grn)) {
    seq(ngenes)
  } else {
    (max(res$.grn$name_map) + 1):ngenes
  }
  gene_range <- setdiff(gene_range, which(rowVars(res$counts_obs) > max.var))
  chosen_genes <- sample(gene_range, floor(ngenes * prob))
  for (i in chosen_genes) {
    # chosen_cells <- sample(which(res$counts_obs[i,] > 0), cell.num)
    chosen_cells <- sample(seq(ncells), cell.num)
    q <- rnorm(1, factor, sd)
    message(sprintf("Gene %d, cells %s, factor %.2f", i, paste(chosen_cells, collapse = ", "), q))
    res$counts_obs[i, chosen_cells] <- res$counts_obs[i, chosen_cells] * q
  }
}
