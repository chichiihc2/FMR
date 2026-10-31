# fmr() against two stored replications of the paper's simulation study (Scenario S2, n = 500,
# Case 2, s_alpha = 0.6; from simulations/results/main_S2 of the repository). The
# stored values were produced on R 4.4.2 (x86_64, reference BLAS); the tolerance allows for
# other BLAS / platforms.
test_that("fmr() reproduces the paper's simulation replications", {
  ref <- readRDS(test_path("fixtures", "reference_S2.rds"))
  for (i in ref$i) {
    r <- ref[ref$i == i, ]
    set.seed(i)
    d <- simulate_fmr_data(n = r$n, s_alpha = r$Dsignal, case = r$case)
    f <- fmr(d$y, d$x, d$Z, pos = d$pos, k = 11, k_max = 11, pve = 0.99, B = 199,
             seed = 5e6 + i, talpha = FALSE, compare = FALSE)
    expect_equal(f$lambda, r$lambda)
    expect_equal(f$beta, r$beta, tolerance = 1e-8)
    expect_equal(f$se_vhat, r$se, tolerance = 1e-8)
    expect_equal(f$c_hat, r$c_hat, tolerance = 1e-6)
    expect_equal(f$se_tau, r$se_tau01, tolerance = 1e-8)
  }
})

test_that("fmr() runs with covariates, BIC k and comparators", {
  set.seed(3)
  d <- simulate_fmr_data(n = 600, s_alpha = 0.3, case = 6)
  cov <- data.frame(age = rnorm(600), sex = factor(sample(1:2, 600, TRUE)))
  f <- fmr(d$y, d$x, d$Z, pos = 1e5 + 1e3 * d$pos, covariates = cov, B = 59)
  expect_s3_class(f, "fmr")
  expect_true(f$k >= 3 && f$k <= 10)
  expect_true(all(c("FMR, bootstrap-t", "FMR, Wald (tau = 0.1)", "OLS") %in% f$coef$method))
  expect_true(is.finite(f$talpha[["p"]]))
  expect_equal(nrow(f$alpha), 101)
  expect_output(print(f), "Functional Mendelian randomization")
  expect_equal(names(coef(f)), "beta")
})

test_that("the caller's random-number stream is restored", {
  set.seed(10); d <- simulate_fmr_data(n = 400)
  set.seed(99); a <- runif(1)
  set.seed(99); invisible(fmr(d$y, d$x, d$Z, pos = d$pos, B = 20)); b <- runif(1)
  expect_identical(a, b)
})
