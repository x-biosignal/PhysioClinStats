library(testthat)
library(PhysioClinStats)

test_that("split-conformal 90% intervals attain marginal coverage (88-92%)", {
  set.seed(1)
  alpha <- 0.1; nsim <- 1000; ntr <- 50; nca <- 50; covered <- 0L
  for (s in seq_len(nsim)) {
    tr <- data.frame(x = rnorm(ntr)); tr$y <- 2 * tr$x + rnorm(ntr)
    ca <- data.frame(x = rnorm(nca)); ca$y <- 2 * ca$x + rnorm(nca)
    fit <- lm(y ~ x, data = tr)
    xte <- rnorm(1); yte <- 2 * xte + rnorm(1)
    v <- PhysioCore::resultValue(
      conformalInterval(fit, ca, data.frame(x = xte), alpha = alpha))
    if (yte >= v$lower && yte <= v$upper) covered <- covered + 1L
  }
  cov <- covered / nsim
  expect_gte(cov, 0.88); expect_lte(cov, 0.93)     # distribution-free guarantee
})

test_that("jackknife+ intervals are valid and near-nominal", {
  set.seed(2)
  alpha <- 0.1; nsim <- 150; n <- 25; covered <- 0L
  for (s in seq_len(nsim)) {
    d <- data.frame(x = rnorm(n)); d$y <- 2 * d$x + rnorm(n)
    fit <- lm(y ~ x, data = d)
    xte <- rnorm(1); yte <- 2 * xte + rnorm(1)
    v <- PhysioCore::resultValue(
      conformalInterval(fit, d, data.frame(x = xte), alpha = alpha,
                        type = "jackknife_plus"))
    expect_true(v$lower <= v$upper)
    if (yte >= v$lower && yte <= v$upper) covered <- covered + 1L
  }
  expect_gte(covered / nsim, 0.84); expect_lte(covered / nsim, 0.97)
})

test_that("the analytic credible interval matches the conjugate closed form", {
  set.seed(3)
  y <- rnorm(30, 5, 2); pm <- 1; ps <- 3; sig <- 2; lev <- 0.95
  v <- PhysioCore::resultValue(
    bayesEstimate(y, prior_mean = pm, prior_sd = ps, sigma = sig, level = lev))
  pp <- 1 / ps^2 + length(y) / sig^2; pv <- 1 / pp
  pmean <- pv * (pm / ps^2 + length(y) * mean(y) / sig^2)
  z <- qnorm((1 + lev) / 2)
  expect_equal(v$posterior_mean, pmean, tolerance = 1e-6)
  expect_equal(v$posterior_sd, sqrt(pv), tolerance = 1e-6)
  expect_equal(v$ci_lower, pmean - z * sqrt(pv), tolerance = 1e-6)
  expect_equal(v$ci_upper, pmean + z * sqrt(pv), tolerance = 1e-6)
  # a flat prior reproduces the sampling posterior N(mean(y), sigma^2/n)
  vf <- PhysioCore::resultValue(bayesEstimate(y, sigma = sig))
  expect_equal(vf$posterior_mean, mean(y), tolerance = 1e-6)
  expect_equal(vf$posterior_sd, sig / sqrt(length(y)), tolerance = 1e-6)
})

test_that("results are valid AnalysisResults with seed provenance", {
  set.seed(4)
  d <- data.frame(x = rnorm(60)); d$y <- d$x + rnorm(60)
  fit <- lm(y ~ x, data = d[1:30, ])
  r <- conformalPredict(fit, d[31:60, ], data.frame(x = 0.5), seed = 42)
  expect_s4_class(r, "AnalysisResult")
  expect_equal(PhysioCore::resultType(r), "conformal")
  expect_true("seed" %in% names(r@provenance))
  expect_equal(r@provenance$seed, 42L)
  b <- bayesEstimate(rnorm(20), seed = 7)
  expect_equal(PhysioCore::resultType(b), "bayes")
  expect_equal(b@provenance$seed, 7L)
})

test_that("the Bayesian formula path degrades gracefully without Stan", {
  skip_if(requireNamespace("rstanarm", quietly = TRUE) ||
          requireNamespace("brms", quietly = TRUE),
          "Stan backend is installed")
  expect_error(bayesEstimate(y ~ 1, data = data.frame(y = rnorm(10))),
               "rstanarm or brms")
  # the analytic numeric path still runs
  expect_s4_class(bayesEstimate(rnorm(15)), "AnalysisResult")
})

test_that("credibleInterval: HDI is narrower than the equal-tailed interval", {
  set.seed(5)
  post <- rexp(5000)                                # skewed posterior
  hdi <- credibleInterval(post, 0.9, "hdi")
  qi <- credibleInterval(post, 0.9, "quantile")
  expect_lt(diff(hdi), diff(qi))
  # symmetric posterior: HDI ~ quantile
  sym <- rnorm(5000)
  expect_equal(diff(credibleInterval(sym, 0.9, "hdi")),
               diff(credibleInterval(sym, 0.9, "quantile")), tolerance = 0.15)
})

test_that("uncertainty functions validate their inputs", {
  d <- data.frame(x = rnorm(40)); d$y <- d$x + rnorm(40)
  fit <- lm(y ~ x, data = d)
  expect_error(conformalInterval(fit, d, data.frame(x = 0), alpha = 1.5),
               "alpha")
  expect_error(conformalInterval(fit, data.frame(x = 1), data.frame(x = 0)),
               "response")                          # calib lacks the response
  expect_error(conformalPredict(fit, d, data.frame(x = c(0, 1))), "one-row")
  expect_error(bayesEstimate(numeric(0)), "finite value")
  expect_error(credibleInterval(1, level = 0.9), "2 finite")
  expect_error(bayesEstimate(rnorm(5), level = 2), "level")
})

test_that("conformal scores a transformed response on the model's scale", {
  # log(y) ~ x: the interval must be built from log-scale residuals, and cover
  # the log-scale response at the nominal rate (not the raw-y mismatch)
  set.seed(9); alpha <- 0.1; nsim <- 500; n <- 40; covered <- 0L
  for (s in seq_len(nsim)) {
    tr <- data.frame(x = rnorm(n)); tr$y <- exp(0.5 * tr$x + rnorm(n, 0, 0.3))
    ca <- data.frame(x = rnorm(n)); ca$y <- exp(0.5 * ca$x + rnorm(n, 0, 0.3))
    fit <- lm(log(y) ~ x, data = tr)
    xte <- rnorm(1); yte <- exp(0.5 * xte + rnorm(1, 0, 0.3))
    v <- PhysioCore::resultValue(conformalInterval(fit, ca, data.frame(x = xte),
                                                   alpha = alpha))
    if (log(yte) >= v$lower && log(yte) <= v$upper) covered <- covered + 1L
  }
  expect_gte(covered / nsim, 0.86)                  # model-scale coverage holds
})

# --- regression tests for adversarial-review findings (WS8-03) -----------------

test_that("conformal rejects calibration rows with missing values", {
  set.seed(1)
  d <- data.frame(x = rnorm(60)); d$y <- 2 * d$x + rnorm(60)
  fit <- lm(y ~ x, data = d)
  ca <- data.frame(x = rnorm(50)); ca$y <- 2 * ca$x + rnorm(50)
  ca$x[c(5, 17, 33)] <- NA                          # missing predictors
  expect_error(conformalInterval(fit, ca, data.frame(x = 0)),
               "missing values")
  ca2 <- data.frame(x = rnorm(50)); ca2$y <- 2 * ca2$x + rnorm(50)
  ca2$y[c(2, 9)] <- NA                              # missing response
  expect_error(conformalInterval(fit, ca2, data.frame(x = 0)), "missing values")
})

test_that("jackknife+ reports an informative error when the model cannot refit", {
  set.seed(1)
  d <- data.frame(x = rnorm(40)); d$y <- d$x + rnorm(40)
  w <- runif(40)                                    # weights var not in calib
  fitw <- lm(y ~ x, data = d, weights = w)
  expect_error(
    conformalInterval(fitw, d, data.frame(x = 0), type = "jackknife_plus"),
    "could not refit")
})
