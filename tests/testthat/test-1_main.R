test_that("simulates true counts using Kinetic model", {
  data(GRN_params_100)

  options_ <- list(
    rand.seed = 0,
    GRN = GRN_params_100,
    num.cells = 100,
    num.cifs = 20,
    tree = Phyla5(),
    diff.cif.fraction = 0.8,
    do.velocity = T
  )

  selectedIndicies <- c(1:5, 1000:1005, 10000:10005)

  res <- sim_true_counts(options_)
  expect_equal(dim(res$counts), c(110, 100))
  expect_equal(
    res$counts[selectedIndicies],
    c(205, 26, 11, 0, 2166, 0, 4, 12, 5, 8, 2, 12, 19, 0, 215, 19, 0)
  )

  add_expr_noise(res)
  divide_batches(res, nbatch = 2)
  expect_equal(
    res$counts_obs[selectedIndicies],
    c(3062, 924, 308, 0, 32756, 0, 0, 284, 0, 251, 0, 557, 554, 0, 3702, 454, 0)
  )
  expect_equal(
    res$counts_with_batches[selectedIndicies],
    c(819, 5412, 461, 0, 251479, 0, 0, 174, 0, 234, 0, 78, 953, 0, 20046, 545, 0)
  )

  expect_equal(
    res$atac_counts[selectedIndicies],
    c(1.1468792, 6.4758649, 0.0000000, 0.2558878, 0.9397685, 2.9793797,
      0.0000000, 1.8782937, 2.1486754, 0.0000000, 1.4174428, 0.0000000,
      0.0000000, 3.1725567, 3.9448536, 0.0000000, 0.0000000)
  )

  expect_equal(
    res$velocity[selectedIndicies],
    c(-11.81783778, 0.07735364, -1.76212608, 0.00000000, -29.35071370,
      0.00000000, -0.77707389, -4.01203054, -0.24126909, 3.04964790,
      -0.21508495, 0.18758518, -0.23687444, 0.54647546, -27.75571875,
      -2.71532070, 0.56252536)
  )
})


test_that("simulates true counts using Beta-Poisson model", {
  data(GRN_params_100)

  options_ <- list(
    rand.seed = 0,
    GRN = GRN_params_100,
    num.cells = 100,
    num.cifs = 20,
    tree = Phyla5(),
    diff.cif.fraction = 0.8,
    do.velocity = F
  )

  res <- sim_true_counts(options_)
  expect_equal(dim(res$counts), c(110, 100))
  expect_equal(
    res$counts[c(1:5, 101:105, 10001:10005)],
    c(186, 14, 3, 0, 2902, 13, 1, 7, 7, 1, 14, 2, 184, 14, 3)
  )
})


