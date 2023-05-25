test_that("simulates data using Kinetic model", {
  data(GRN_params_100, envir = environment())

  set.seed(0)
  options_ <- list(
    GRN = GRN_params_100,
    num.cells = 100,
    num.cifs = 20,
    tree = Phyla5(),
    diff.cif.fraction = 0.8,
    do.velocity = T
  )

  res <- sim_true_counts(options_)

  selectedIndicies <- c(1:5, 1000:1005, 10000:10005)
  expect_equal(dim(res$counts), c(110, 100))
  expect_equal(
    res$counts[selectedIndicies],
    c(34, 5, 28, 21, 49, 0, 0, 10, 4, 90, 18, 18, 0, 12, 0, 37, 171)
  )

  add_expr_noise(res)
  divide_batches(res, nbatch = 2)

  expect_equal(
    res$counts_obs[selectedIndicies],
    c(585, 307, 141, 187, 309, 0, 0, 326, 0, 2692, 0, 401, 0, 22, 0, 187, 2291)
  )
  expect_equal(
    res$counts_with_batches[selectedIndicies],
    c(2331, 5031, 180, 1263, 131, 0, 0, 93, 0, 21462, 0, 467, 0, 357, 0,
      1020, 1495)
  )

  expect_equal(
    res$atac_counts[selectedIndicies],
    c(0.0000000, 0.0000000, 0.6887022, 0.7694152, 0.0000000, 2.5172858,
      0.3549220, 5.6001554, 3.5924741, 0.2475742, 0.0000000, 1.0374884,
      3.4453020, 4.0073263, 0.0000000, 0.0000000, 0.0000000)
  )

  expect_equal(
    res$velocity[selectedIndicies],
    c(4.3146894, 0.4199234, -2.1558627, 2.0012763, -23.2763955, 0.0000000,
      0.0000000, -2.0025398, 1.2857669, 2.8214759, 2.5935024, 1.3315967,
      0.0000000, -2.5724074, 0.2749523, 1.9946433, -2.0456734)
  )

  expect_no_error(plot_gene_module_cor_heatmap(res, save = FALSE))
  expect_no_error(gene_corr_regulator(res, 2))
  expect_no_error(plot_rna_velocity(res, perplexity = 20))
})


test_that("simulates data using Beta-Poisson model", {
  data(GRN_params_100, envir = environment())

  set.seed(0)
  options_ <- list(
    GRN = GRN_params_100,
    num.cells = 100,
    num.cifs = 20,
    tree = Phyla5(),
    diff.cif.fraction = 0.8,
    do.velocity = F
  )

  res <- sim_true_counts(options_)

  selectedIndicies <- c(1:5, 101:105, 10001:10005)
  expect_equal(dim(res$counts), c(110, 100))
  expect_equal(
    res$counts[selectedIndicies],
    c(120, 5, 18, 33, 88, 5, 4, 0, 8, 96, 0, 18, 0, 15, 146)
  )
})


test_that("simulates spatial data", {
  data(GRN_params_100, envir = environment())

  lig_params <- data.frame(
    target    = c(101, 102),
    regulator = c(103, 104),
    effect    = c(5.2, 5.9)
  )

  options_ <- list2(
    GRN = GRN_params_100,
    num.genes = 200,
    num.cells = 100,
    num.cifs = 20,
    tree = Phyla3(),
    intrinsic.noise = 0.5,
    cci = list(
      params = lig_params,
      max.neighbors = 4,
      cell.type.interaction = "random",
      step.size = 0.5
    )
  )

  set.seed(0)
  res <- sim_true_counts(options_)

  selectedIndicies <- c(1:5, 1000:1005, 10000:10005)
  expect_equal(
    res$counts[selectedIndicies],
    c(40.675564, 30.876988, 29.984167, 49.430348, 25.113605, 4.093944,
      45.194247, 29.063519, 47.389263, 42.516067, 43.014273, 7.110385,
      55.992341, 13.604489, 14.811897, 10.213004, 24.046141)
  )

  expect_no_error(plot_cell_loc(res))
  expect_no_error(gene_corr_cci(res))
})


test_that("simulates spatial data with discrete population and HGE", {
  data(GRN_params_100, envir = environment())

  lig_params <- data.frame(
    target    = c(101, 102),
    regulator = c(103, 104),
    effect    = c(5.2, 5.9)
  )

  options_ <- list2(
    GRN = GRN_params_100,
    num.genes = 200,
    num.cells = 100,
    num.cifs = 20,
    tree = Phyla3(),
    discrete.cif = T,
    discrete.min.pop.size = 20,
    intrinsic.noise = 0.5,
    hge.prop = 0.05,
    cci = list(
      params = lig_params,
      max.neighbors = 4,
      cell.type.interaction = "random",
      step.size = 0.5
    )
  )

  set.seed(0)
  res <- sim_true_counts(options_)

  selectedIndicies <- c(1:5, 1000:1005, 10000:10005)
  expect_equal(
    res$counts[selectedIndicies],
    c(109.0693303, 60.4151790, 91.9120934, 0.5816326, 177.6741585, 197.0663584,
      102.3145704, 65.9978484, 89.1613630, 3.1734446, 156.9202179, 29.8315553,
      92.3944947, 66.1421921, 105.4677530, 0.5729707, 110.5115346)
  )
})
