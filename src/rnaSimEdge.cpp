#include <Rcpp.h>
#include <algorithm>
#include <numeric>
#include <vector>

using namespace Rcpp;

// ---------------------------------------------------------------------------
// sample_density: C++ reimplementation of SampleDen(reduce.mem = TRUE)
// ---------------------------------------------------------------------------
static std::vector<double> sample_density(int nsample,
                                          const NumericVector &den_x,
                                          const NumericVector &den_y) {
  int nbins = den_x.size();
  double bw = den_x[1] - den_x[0];

  // Normalize y to probabilities
  double sum_y = 0.0;
  for (int i = 0; i < nbins; i++) sum_y += den_y[i];
  std::vector<double> probs(nbins);
  for (int i = 0; i < nbins; i++) probs[i] = den_y[i] / sum_y;

  // Bin boundaries
  std::vector<double> mins(nbins), maxs(nbins);
  for (int i = 0; i < nbins; i++) {
    mins[i] = den_x[i] - 0.5 * bw;
    maxs[i] = den_x[i] + 0.5 * bw;
  }

  // Draw multinomial bin counts using R's API
  IntegerVector counts(nbins);
  R_CheckUserInterrupt();
  ::Rf_rmultinom(nsample, probs.data(), nbins, counts.begin());

  // Total samples (should equal nsample)
  int total = 0;
  for (int i = 0; i < nbins; i++) total += counts[i];

  // Draw uniform samples within each bin
  std::vector<double> samples(total);
  int idx = 0;
  for (int i = 0; i < nbins; i++) {
    int c = counts[i];
    for (int j = 0; j < c; j++) {
      samples[idx++] = R::runif(mins[i], maxs[i]);
    }
  }

  return samples;
}

// ---------------------------------------------------------------------------
// match_params_den: C++ reimplementation of .matchParamsDen()
// ---------------------------------------------------------------------------
static void match_params_den(double *x_out,           // output: length n_x
                             const double *x_in,       // input values, length n_x
                             int n_x,
                             std::vector<double> &prev_values,  // modified in place
                             int keep_prev_val,
                             const NumericVector &den_x,
                             const NumericVector &den_y) {
  // Concatenate x_in with prev_values
  int n_prev = (int)prev_values.size();
  int n_total = n_x + n_prev;
  std::vector<double> values(n_total);
  std::copy(x_in, x_in + n_x, values.begin());
  std::copy(prev_values.begin(), prev_values.end(), values.begin() + n_x);

  // Compute ranks matching R's rank() with ties.method = "average"
  std::vector<int> order(n_total);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](int a, int b) {
    return values[a] < values[b];
  });
  std::vector<double> ranks(n_total);
  {
    int i = 0;
    while (i < n_total) {
      int j = i;
      // Find the end of the group of tied values
      while (j < n_total - 1 && values[order[j + 1]] == values[order[j]]) {
        j++;
      }
      // Average rank for positions i..j (1-based: i+1..j+1)
      double avg_rank = (i + 1 + j + 1) / 2.0;
      for (int k = i; k <= j; k++) {
        ranks[order[k]] = avg_rank;
      }
      i = j + 1;
    }
  }

  // Sample from density — max(ranks) == n_total when no ties,
  // but with average ties it could be fractional. Use ceil to be safe.
  int max_rank = n_total; // n_total is always the correct max for average ranks
  std::vector<double> samples = sample_density(max_rank, den_x, den_y);

  // Sort the samples
  std::sort(samples.begin(), samples.end());

  // Map back via ranks for the x_in portion only
  // R truncates fractional indices: sorted[141.5] == sorted[141]
  for (int i = 0; i < n_x; i++) {
    x_out[i] = samples[(int)ranks[i] - 1];
  }

  // Update prev_values buffer
  if (n_prev < keep_prev_val) {
    prev_values.resize(n_total);
    std::copy(values.begin(), values.end(), prev_values.begin());
  }
}

// ---------------------------------------------------------------------------
// rnaSimEdgeCpp: main exported function
// Replaces the per-cell for-loop in .rnaSimEdge() for the common fast path:
//   - no velocity, no dynamic GRN, continuous CIF, no custom params_mpl_fn$s
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List rnaSimEdgeCpp(
    NumericMatrix s_base,       // n_cell x n_gene (full matrix, indexed by cell_idx)
    NumericMatrix kon,          // n_gene x n_cell
    NumericMatrix koff,         // n_gene x n_cell
    NumericMatrix geff,         // n_gene x n_reg (or 0x0 if no_grn)
    IntegerVector regulators,   // 1-based gene indices (empty if no_grn)
    NumericVector curr_cif_in,  // length n_reg (empty if no_grn)
    double grn_effect,
    NumericVector hge_scale,    // n_gene (or length 1 scalar)
    NumericVector scale_s,      // length 1 (scalar) or n_cluster (per-cluster)
    IntegerVector cell_pop,     // 1-based cluster index per cell (length n_cell; empty if scalar scale_s)
    double intr_noise,
    NumericVector den_x,        // density grid x
    NumericVector den_y,        // density grid y
    IntegerVector cell_idx,     // 1-based R indices of cells on this edge
    NumericVector prev_values_in, // matchParamsDen history buffer
    int keep_prev_val,          // max history size
    bool no_grn
) {
  int n_cells_on_edge = cell_idx.size();
  int n_gene = s_base.ncol();
  int n_reg = no_grn ? 0 : curr_cif_in.size();
  bool hge_scalar = (hge_scale.size() == 1);
  bool scale_s_scalar = (scale_s.size() == 1);

  // Copy curr_cif so we can mutate it
  std::vector<double> curr_cif(n_reg);
  for (int i = 0; i < n_reg; i++) curr_cif[i] = curr_cif_in[i];

  // Copy prev_values into a std::vector for in-place updates
  std::vector<double> prev_values(prev_values_in.begin(), prev_values_in.end());

  // Allocate output matrices
  NumericMatrix counts_s(n_cells_on_edge, n_gene);
  NumericMatrix params_s(n_gene, n_cells_on_edge);
  NumericMatrix cif_regu;
  if (!no_grn) {
    cif_regu = NumericMatrix(n_cells_on_edge, n_reg);
  }

  // Temp buffers
  std::vector<double> raw_s(n_gene);
  std::vector<double> s_cell(n_gene);

  for (int n = 0; n < n_cells_on_edge; n++) {
    int i_cell = cell_idx[n] - 1; // 0-based index

    // 1. Compute raw_s = s_base[i_cell,] + curr_cif * geff * grn_effect
    for (int g = 0; g < n_gene; g++) {
      raw_s[g] = s_base(i_cell, g);
    }
    if (!no_grn) {
      // R: curr_cif %*% t(geff) where geff is n_gene x n_reg
      // = (1 x n_reg) %*% (n_reg x n_gene)
      for (int g = 0; g < n_gene; g++) {
        double grn_term = 0.0;
        for (int r = 0; r < n_reg; r++) {
          grn_term += curr_cif[r] * geff(g, r);
        }
        raw_s[g] += grn_term * grn_effect;
      }
    }

    // 2. match_params_den
    match_params_den(s_cell.data(), raw_s.data(), n_gene,
                     prev_values, keep_prev_val, den_x, den_y);

    // 3. Scale: s_cell = 10^s_cell * scale_s_cell * hge_scale
    //    scale_s may be scalar or per-cluster (indexed by cell_pop)
    //    hge_scale may be scalar (length 1) or per-gene (length n_gene)
    double ss = scale_s_scalar ? scale_s[0] : scale_s[cell_pop[i_cell] - 1];
    for (int g = 0; g < n_gene; g++) {
      double hge = hge_scalar ? hge_scale[0] : hge_scale[g];
      s_cell[g] = std::pow(10.0, s_cell[g]) * ss * hge;
    }

    // 4. Beta-Poisson sampling
    //    Match R's .betaPoisson.vec RNG order: ALL betas first, THEN all poissons
    {
      std::vector<double> y_beta(n_gene);
      for (int g = 0; g < n_gene; g++) {
        y_beta[g] = R::rbeta(kon(g, i_cell), koff(g, i_cell));
      }
      for (int g = 0; g < n_gene; g++) {
        double k_on = kon(g, i_cell);
        double k_off = koff(g, i_cell);
        double yMean = k_on / (k_on + k_off);
        double xMean = yMean * s_cell[g];
        double x = R::rpois(y_beta[g] * s_cell[g]);
        counts_s(n, g) = intr_noise * x + (1.0 - intr_noise) * xMean;
      }
    }

    // 5. Store s_cell in params_s[, n]
    for (int g = 0; g < n_gene; g++) {
      params_s(g, n) = s_cell[g];
    }

    // 6. CIF update (if !no_grn): continuous case
    if (!no_grn) {
      // counts_regu = counts[regulators]
      // mean_counts = mean(counts)
      double sum_counts = 0.0;
      for (int g = 0; g < n_gene; g++) {
        sum_counts += counts_s(n, g);
      }
      double mean_counts = sum_counts / n_gene;

      for (int r = 0; r < n_reg; r++) {
        int reg_idx = regulators[r] - 1; // 0-based
        double count_r = counts_s(n, reg_idx);
        curr_cif[r] = count_r / (count_r + mean_counts);
      }

      // Store in cif_regu
      for (int r = 0; r < n_reg; r++) {
        cif_regu(n, r) = curr_cif[r];
      }
    }
  }

  // Build return list
  List result;
  result["counts_s"] = counts_s;
  result["params_s"] = params_s;
  if (!no_grn) {
    result["cif_regu"] = cif_regu;
  }
  // Convert prev_values back to NumericVector
  result["prev_values"] = NumericVector(prev_values.begin(), prev_values.end());

  return result;
}
