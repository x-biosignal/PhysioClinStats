library(testthat)
library(PhysioClinStats)

causal_rv <- function(x) PhysioCore::resultValue(x)

.mediation_models <- function() {
  data(jobs, package = "mediation")
  list(
    data = jobs,
    mediator = stats::lm(
      job_seek ~ treat + econ_hard + sex + age, data = jobs),
    outcome = stats::lm(
      depress2 ~ treat + job_seek + econ_hard + sex + age, data = jobs)
  )
}

test_that("causalMediation maps every backend effect and preserves RNG", {
  skip_if_not_installed("mediation")
  models <- .mediation_models()
  set.seed(80211)
  before <- .Random.seed
  result <- causalMediation(
    models$mediator, models$outcome,
    treat = "treat", mediator = "job_seek",
    sims = 100L, seed = 481L, sensitivity = FALSE)
  expect_identical(.Random.seed, before)

  set.seed(481)
  reference <- mediation::mediate(
    models$mediator, models$outcome,
    treat = "treat", mediator = "job_seek",
    control.value = 0, treat.value = 1,
    sims = 100L, boot = FALSE, boot.ci.type = "perc",
    conf.level = 0.95, robustSE = FALSE)
  table <- causal_rv(result)$effects
  expected <- c(
    reference$d0, reference$d1, reference$d.avg,
    reference$z0, reference$z1, reference$z.avg,
    reference$tau.coef,
    reference$n0, reference$n1, reference$n.avg)
  expect_equal(table$estimate, expected, tolerance = 1e-12)
  expect_equal(table$conf_low, unname(c(
    reference$d0.ci[1], reference$d1.ci[1], reference$d.avg.ci[1],
    reference$z0.ci[1], reference$z1.ci[1], reference$z.avg.ci[1],
    reference$tau.ci[1],
    reference$n0.ci[1], reference$n1.ci[1], reference$n.avg.ci[1])),
  tolerance = 1e-12)
  expect_equal(
    unname(result@estimate),
    c(reference$d.avg, reference$z.avg, reference$tau.coef),
    tolerance = 1e-12)
  expect_equal(PhysioCore::resultType(result), "causal_mediation")
})

test_that("causalMediation is reproducible and accepts AnalysisResult fits", {
  skip_if_not_installed("mediation")
  models <- .mediation_models()
  wrapped_m <- PhysioCore::AnalysisResult(
    "model", result = list(fit = models$mediator))
  wrapped_y <- PhysioCore::AnalysisResult(
    "model", result = list(fit = models$outcome))
  first <- causalMediation(
    wrapped_m, wrapped_y, "treat", "job_seek",
    sims = 75L, seed = 19L, sensitivity = FALSE)
  second <- causalMediation(
    wrapped_m, wrapped_y, "treat", "job_seek",
    sims = 75L, seed = 19L, sensitivity = FALSE)
  expect_identical(causal_rv(first)$effects, causal_rv(second)$effects)
  expect_equal(causal_rv(first)$sample_size, nrow(models$data))
})

test_that("causalMediation sensitivity matches medsens for linear models", {
  skip_if_not_installed("mediation")
  models <- .mediation_models()
  set.seed(117)
  rng_before <- .Random.seed
  result <- causalMediation(
    models$mediator, models$outcome, "treat", "job_seek",
    sims = 60L, seed = 31L, sensitivity = TRUE)
  expect_identical(.Random.seed, rng_before)
  sensitivity <- causal_rv(result)$sensitivity
  expect_equal(sensitivity$status, "computed")

  set.seed(31)
  direct <- mediation::mediate(
    models$mediator, models$outcome, treat = "treat",
    mediator = "job_seek", sims = 60L)
  direct_sensitivity <- mediation::medsens(direct, effect.type = "both")
  expect_equal(sensitivity$curves$rho, direct_sensitivity$rho,
               tolerance = 1e-12)
  expect_equal(sensitivity$curves$acme_control,
               as.numeric(direct_sensitivity$d0), tolerance = 1e-12)
  expect_equal(sensitivity$zero_crossing$acme,
               direct_sensitivity$err.cr.d, tolerance = 1e-12)
})

test_that("pinned mediation fixture reproduces backend effects and sensitivity", {
  skip_if_not_installed("mediation")
  fixture <- readRDS(test_path(
    "fixtures", "causal", "mediation-jobs.rds"))
  model_m <- stats::lm(
    job_seek ~ treat + econ_hard + sex + age, data = fixture$data)
  model_y <- stats::lm(
    depress2 ~ treat + job_seek + econ_hard + sex + age,
    data = fixture$data)
  result <- causalMediation(
    model_m, model_y, treat = "treat", mediator = "job_seek",
    sims = fixture$settings$sims, seed = fixture$settings$seed,
    sensitivity = TRUE)
  value <- causal_rv(result)
  expect_equal(
    value$effects$estimate,
    c(
      fixture$effects$d0, fixture$effects$d1, fixture$effects$d_avg,
      fixture$effects$z0, fixture$effects$z1, fixture$effects$z_avg,
      fixture$effects$total,
      value$backend$n0, value$backend$n1, value$backend$n.avg
    ),
    tolerance = 1e-12)
  expect_equal(value$sensitivity$curves$rho,
               fixture$sensitivity$rho, tolerance = 1e-12)
  # The medsens sensitivity curves (acme_control = d0, ade_control = z0) are
  # analytic functions of the two lm fits, computed through LAPACK matrix/QR
  # operations. They are RNG-invariant (identical across seeds) but subject to
  # cross-platform BLAS/LAPACK drift. Against a different BLAS than the one that
  # pinned this fixture (measured under OpenBLAS vs reference LAPACK) the curves
  # differ by up to 8.0e-13 absolute / 2.1e-10 relative -- above the 1e-12 used
  # by the RNG-reproduced quantities, but far below any meaningful magnitude
  # (|d0|,|z0| ~ 1e-2..1.5e-1). Assert ~8-significant-figure reproduction.
  expect_equal(value$sensitivity$curves$acme_control,
               as.numeric(fixture$sensitivity$d0), tolerance = 1e-8)
  expect_equal(value$sensitivity$curves$ade_control,
               as.numeric(fixture$sensitivity$z0), tolerance = 1e-8)
  expect_equal(value$sensitivity$zero_crossing$acme,
               fixture$sensitivity$err_cr_d, tolerance = 1e-12)
  expect_equal(value$sensitivity$zero_crossing$ade,
               fixture$sensitivity$err_cr_z, tolerance = 1e-12)
})

test_that("unsupported sensitivity is structured rather than fatal", {
  skip_if_not_installed("mediation")
  models <- .mediation_models()
  binary <- transform(models$data, depressed = as.integer(depress2 > 3))
  model_y <- stats::glm(
    depressed ~ treat + job_seek + econ_hard + sex + age,
    family = stats::binomial(), data = binary)
  model_m <- stats::lm(
    job_seek ~ treat + econ_hard + sex + age, data = binary)
  result <- causalMediation(
    model_m, model_y, "treat", "job_seek",
    sims = 50L, seed = 2L, sensitivity = TRUE)
  expect_equal(causal_rv(result)$sensitivity$status, "not_applicable")
  expect_match(causal_rv(result)$sensitivity$reason, "linear")
})

test_that("causalMediation rejects broken model and setting contracts", {
  skip_if_not_installed("mediation")
  models <- .mediation_models()
  reversed <- models$data[nrow(models$data):1L, ]
  row.names(reversed) <- sprintf("r%04d", seq_len(nrow(reversed)))
  model_y_reversed <- stats::lm(
    depress2 ~ treat + job_seek + econ_hard + sex + age,
    data = reversed)
  expect_error(causalMediation(
    models$mediator, model_y_reversed, "treat", "job_seek",
    sims = 20L, sensitivity = FALSE), "same observations")
  expect_error(causalMediation(
    models$mediator, models$outcome, "missing", "job_seek",
    sims = 20L), "both model")
  expect_error(causalMediation(
    models$mediator, models$outcome, "treat", "wrong",
    sims = 20L), "response")
  expect_error(causalMediation(
    models$mediator, models$outcome, "treat", "job_seek",
    control.value = 5, sims = 20L), "observed treatment")
  expect_error(causalMediation(
    models$mediator, models$outcome, "treat", "job_seek",
    covariates = list(nope = 1), sims = 20L), "absent")
  expect_error(causalMediation(
    models$mediator, models$outcome, "treat", "job_seek",
    sims = 0), "positive integer")
  expect_error(causalMediation(
    models$mediator, models$outcome, "treat", "job_seek",
    sims = 20L, conf.level = 1), "strictly between")
  expect_error(causalMediation(
    models$mediator, models$outcome, "treat", "job_seek",
    sims = 20L, covariates = list(age = Inf)), "finite")
  expect_error(causalMediation(
    models$mediator, models$outcome, "treat", "job_seek",
    sims = 20L, boot.ci.type = "wrong"), "perc.*bca")
  expect_error(causalMediation(
    models$mediator, models$outcome, "treat", "job_seek",
    sims = 20L, boot = TRUE, robustSE = TRUE),
  "not available")
  expect_error(causalMediation(
    PhysioCore::AnalysisResult("empty"), models$outcome,
    "treat", "job_seek"), "without result")
})
