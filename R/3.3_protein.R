.proteinSeq <- function(sim) {
  N <- sim$N
  options <- sim$options

  gene_to_prot <- OP("gene.to.prot")
  gene_to_prot <- if (is.null(gene_to_prot)) {
    seq_len(max(min(N$gene %/% 10, 1), 30))
  } else {
    gene_to_prot
  }
  sim$gene_to_prot <- gene_to_prot
  n_prot <- length(gene_to_prot)
  if (n_prot < 1) {
    return()
  }

  # "g1", "g4", "g3", ...
  prot_idx <- sim$gene_name_map[gene_to_prot]

  noise_sd <- OP("prot.noise")

  cif <- sim$CIF_spatial_s %||% sim$CIF_all$cif$s
  piv <- sim$GIV$s[prot_idx,]
  piv <- piv + matrix(rnorm(length(piv), 0, noise_sd), ncol = N$cif)

  x <- cif %*% t(piv)
  alpha_mean <- OP("prot.alpha.mean")
  alpha_sd <- OP("prot.alpha.sd")
  theta_mean <- OP("prot.theta.mean")
  theta_sd <- OP("prot.theta.sd")
  dist_alpha <- rnorm(n_prot, alpha_mean, alpha_sd)
  dist_theta <- .rnormTrunc(n_prot, theta_mean, theta_sd, 5, 400)

  counts <- lapply(seq_len(n_prot), function(i) {
    x_tg <- sort(rgamma(N$cell, shape = dist_alpha[i], scale = dist_theta[i]))
    x_tg[rank(x[, i])]
  }) %>% do.call(cbind, .)

  sim$counts_prot <- counts
}


#' Infer subsampling parameters from a real ADT count matrix.
#'
#' @param counts  Integer matrix, cells x proteins (raw ADT counts).
#' @param simulated  Integer matrix, cells x proteins (clean simulated counts).
#' @param target_zf  Target zero fraction to match.
#' @param verbose  Print fitting diagnostics.
#' @return A list used by protein_subsample().
#' @export
infer_protein_subsample_params <- function(
  counts, simulated = NULL, target_zf = NULL, verbose = TRUE
) {
  counts <- as.matrix(counts)
  stopifnot(is.numeric(counts), ncol(counts) >= 3)
  np <- ncol(counts)
  nc <- nrow(counts)

  # ---- per-protein summary stats --------------------------------
  pmeans  <- colMeans(counts)
  pzero   <- colMeans(counts == 0)
  pvars   <- apply(counts, 2, var)

  nz_means <- vapply(seq_len(np), function(i) {
    x <- counts[counts[, i] > 0, i]
    if (length(x) > 2) mean(x) else NA_real_
  }, numeric(1))

  # ---- Stage-1: capture efficiency ------------------------------
  # For proteins with enough non-zero signal, the ratio
  #   overall_mean / nonzero_mean ≈ (1 - zero_frac)
  # is driven by both capture loss and dropout.  We isolate the
  # capture component by restricting to high-abundance proteins
  # where dropout is minimal.
  stable <- !is.na(nz_means) & pzero < 0.5 & pmeans > quantile(pmeans, 0.25)
  if (sum(stable) < 3) stable <- !is.na(nz_means) & pmeans > 0

  # capture_rate ≈ median( overall_mean / nz_mean ) among stable proteins
  ratios <- pmeans[stable] / nz_means[stable]
  capture_rate <- median(ratios, na.rm = TRUE)
  capture_rate <- min(max(capture_rate, 0.05), 0.95)

  # ---- Stage-2: dropout logistic curve --------------------------
  # Across proteins, fit   P(zero) = sigmoid(shape * (mid - log1p(mean)))
  lmu <- log1p(pmeans)

  dropout_shape <- 2.0
  dropout_mid   <- median(lmu)

  # Primary: nonlinear least squares
  nls_ok <- FALSE
  tryCatch({
    fit <- nls(
      pzero ~ 1 / (1 + exp(-shape * (mid - lmu))),
      start   = list(shape = 2, mid = median(lmu)),
      control = nls.control(maxiter = 500, warnOnly = TRUE)
    )
    dropout_shape <- coef(fit)[["shape"]]
    dropout_mid   <- coef(fit)[["mid"]]
    nls_ok <- TRUE
  }, error = function(e) NULL)

  # Fallback: OLS on logit-transformed zero fractions
  if (!nls_ok) {
    safe_pz <- pmin(pmax(pzero, 0.01), 0.99)
    logit_z <- log(safe_pz / (1 - safe_pz))
    lfit <- lm(logit_z ~ lmu)
    dropout_shape <- max(-coef(lfit)[[2]], 0.1)
    dropout_mid   <- -coef(lfit)[[1]] / coef(lfit)[[2]]
    if (verbose) message("  NLS failed; used logistic-OLS fallback.")
  }

  # ---- per-protein NB dispersion (method of moments) ------------
  theta <- ifelse(
    pvars > pmeans & pmeans > 0,
    pmeans^2 / (pvars - pmeans),
    100   # near-Poisson default
  )
  theta <- pmin(theta, 1000)
  names(theta) <- colnames(counts)

  params <- structure(list(
    capture_rate      = as.numeric(capture_rate),
    dropout_shape     = as.numeric(dropout_shape),
    dropout_mid       = as.numeric(dropout_mid),
    protein_means     = pmeans,
    protein_zero_frac = pzero,
    protein_vars      = pvars,
    protein_dispersions = theta,
    n_proteins        = np,
    n_cells           = nc
  ), class = "adt_subsample_params")

  if (verbose) message("Calibrating dropout to match target zero fraction...")

  if (!is.null(simulated) && !is.null(target_zf)) {
    params <- .calibrateDropout(simulated, params, target_zf)
  }

  if (verbose) print(params)
  params
}


.printProtSubsampleParams <- function(x, ...) {
  cat("ADT Subsampling Parameters\n")
  cat(sprintf("  Proteins : %d  |  Cells : %d\n", x$n_proteins, x$n_cells))
  cat(sprintf("  Capture rate       : %.3f\n", x$capture_rate))
  cat(sprintf("  Dropout shape      : %.2f\n", x$dropout_shape))
  cat(sprintf("  Dropout midpoint   : %.2f  (log1p scale)\n", x$dropout_mid))
  cat(sprintf("  Median zero frac.  : %.3f\n", median(x$protein_zero_frac)))
  cat(sprintf("  Dispersion range   : [%.1f, %.1f]\n",
              min(x$protein_dispersions), max(x$protein_dispersions)))
  invisible(x)
}


.calibrateDropout <- function(sim_counts, params, target_zf,
                              tol = 0.01, max_iter = 20, seed = 789) {
  lo <- params$dropout_mid - 3

  hi <- params$dropout_mid + 3

  eval_zf <- function(mid_val) {
    p2 <- params
    p2$dropout_mid <- mid_val
    sub <- protein_subsample(sim_counts, p2, seed = seed)
    mean(sub == 0)
  }

  for (iter in seq_len(max_iter)) {
    mid_try <- (lo + hi) / 2
    zf <- eval_zf(mid_try)
    if (abs(zf - target_zf) < tol) break
    if (zf < target_zf) lo <- mid_try else hi <- mid_try
  }

  params_cal <- params
  params_cal$dropout_mid <- mid_try
  message(sprintf("  Calibrated dropout_mid: %.3f -> %.3f  (zf: %.3f -> %.3f)",
                  params$dropout_mid, mid_try, mean(sim_counts == 0), zf))
  params_cal
}


#' Apply two-stage subsampling to a clean simulated protein count matrix.
#'
#' Stage 1 – Binomial thinning:
#'   Each molecule is independently retained with probability `capture_rate`.
#' Stage 2 – Expression-dependent dropout:
#'   After thinning, cells with low remaining counts have an additional
#'   probability of being masked to zero, controlled by the logistic
#'   curve  P(drop | x) = sigmoid(shape * (mid - log1p(x))).
#'
#' @param sim_counts  Integer matrix, cells x proteins (clean simulated counts).
#' @param params      Output of infer_protein_subsample_params().
#' @param capture_rate  Override the inferred capture rate (NULL = use inferred).
#' @param seed        Random seed for reproducibility (NULL = no seed).
#' @return Integer matrix of the same dimensions with realistic zero inflation.
#' @export
protein_subsample <- function(sim_counts, params,
                              capture_rate = NULL, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)
  sim_counts <- as.matrix(sim_counts)
  n_cell <- nrow(sim_counts)
  n_prot <- ncol(sim_counts)

  p_cap <- if (!is.null(capture_rate)) capture_rate else params$capture_rate

  # --- Stage 1: Binomial thinning --------------------------------
  thinned <- matrix(
    rbinom(n    = n_cell * n_prot,
           size = as.integer(pmax(sim_counts, 0)),
           prob = p_cap),
    nrow = n_cell, ncol = n_prot,
    dimnames = dimnames(sim_counts)
  )

  # --- Stage 2: Expression-dependent dropout ---------------------
  log_expr <- log1p(thinned)
  p_drop   <- 1 / (1 + exp(-params$dropout_shape *
                             (params$dropout_mid - log_expr)))

  dropout_mask <- matrix(
    rbinom(n_cell * n_prot, size = 1L, prob = as.numeric(p_drop)),
    nrow = n_cell, ncol = n_prot
  )

  result <- thinned * (1L - dropout_mask)
  dimnames(result) <- dimnames(sim_counts)
  storage.mode(result) <- "integer"
  result
}

