library(testthat)
library(PhysioClinStats)

rv <- function(r) PhysioCore::resultValue(r)
est <- function(r) rv(r)$estimate

# canonical A/B example used throughout SingleCaseES / scan documentation
A <- c(20, 20, 26, 25, 22, 23)
B <- c(28, 25, 24, 27, 30, 30, 29, 28)

test_that("non-overlap effect sizes match SingleCaseES within 1e-6", {
  skip_if_not_installed("SingleCaseES")
  expect_equal(est(scedPND(A, B)), SingleCaseES::PND(A, B)$Est, tolerance = 1e-6)
  expect_equal(est(scedPEM(A, B)), SingleCaseES::PEM(A, B)$Est, tolerance = 1e-6)
  expect_equal(est(scedNAP(A, B)), SingleCaseES::NAP(A, B)$Est, tolerance = 1e-6)
  expect_equal(est(scedTau(A, B)), SingleCaseES::Tau(A, B)$Est, tolerance = 1e-6)
  expect_equal(est(scedTauU(A, B, method = "parker")),
               SingleCaseES::Tau_U(A, B)$Est, tolerance = 1e-6)
})

test_that("NAP standard error and Newcombe interval match SingleCaseES", {
  skip_if_not_installed("SingleCaseES")
  mine <- rv(scedNAP(A, B, confidence = 0.95))
  ref <- SingleCaseES::NAP(A, B, confidence = 0.95)
  expect_equal(mine$se, ref$SE, tolerance = 1e-6)
  expect_equal(mine$ci_lower, ref$CI_lower, tolerance = 1e-6)
  expect_equal(mine$ci_upper, ref$CI_upper, tolerance = 1e-6)
})

test_that("effect sizes match SingleCaseES on the bundled McKissick data", {
  skip_if_not_installed("SingleCaseES")
  utils::data("McKissick", package = "SingleCaseES")
  for (case in unique(McKissick$Case_pseudonym)) {
    d <- McKissick[McKissick$Case_pseudonym == case, ]
    a <- d$Outcome[d$Condition == "A"]; b <- d$Outcome[d$Condition == "B"]
    # McKissick is a reduction outcome
    expect_equal(est(scedNAP(a, b, improvement = "decrease")),
                 SingleCaseES::NAP(a, b, improvement = "decrease")$Est,
                 tolerance = 1e-6)
    expect_equal(est(scedTauU(a, b, improvement = "decrease")),
                 SingleCaseES::Tau_U(a, b, improvement = "decrease")$Est,
                 tolerance = 1e-6)
    expect_equal(est(scedPND(a, b, improvement = "decrease")),
                 SingleCaseES::PND(a, b, improvement = "decrease")$Est,
                 tolerance = 1e-6)
  }
})

test_that("Tau-U scan method matches scan::tau_u estimate, Z and p", {
  skip_if_not_installed("scan")
  sc <- scan::scdf(c(A, B), phase_design = c(A = length(A), B = length(B)))
  row <- scan::tau_u(sc)$table[[1]]["A vs. B - Trend A", ]
  mine <- rv(scedTauU(A, B, method = "scan"))
  expect_equal(mine$estimate, row$Tau, tolerance = 1e-4)
  expect_equal(mine$z, row$Z, tolerance = 1e-4)
  expect_equal(mine$p_value, row$p, tolerance = 1e-4)
})

test_that("scedTauU reports the same S under both denominator conventions", {
  a <- rv(scedTauU(A, B, method = "parker"))
  b <- rv(scedTauU(A, B, method = "scan"))
  expect_equal(a$S, b$S)                       # numerator is convention-free
  expect_equal(a$S, a$S_AB - a$S_trendA)
  expect_false(isTRUE(all.equal(a$estimate, b$estimate)))  # denominators differ
})

test_that("results are AnalysisResults carrying an analytic uncertainty", {
  r <- scedNAP(A, B)
  expect_s4_class(r, "AnalysisResult")
  expect_equal(PhysioCore::resultType(r), "sced_nap")
  u <- r@uncertainty
  expect_equal(u$type, "analytic")
  expect_equal(u$lower, rv(r)$ci_lower)
  expect_equal(u$level, 0.95)
})

test_that("the 2-SD band fires on a shift and stays silent on null data", {
  base <- c(10, 12, 11, 9, 10, 8)
  shifted <- scedTwoSDBand(base, c(15, 16, 17, 16, 18))
  expect_true(est(shifted))
  expect_equal(rv(shifted)$first_run_at, 1L)
  expect_equal(rv(shifted)$upper, mean(base) + 2 * sd(base))

  null <- scedTwoSDBand(base, c(10, 11, 9, 10, 12))
  expect_false(est(null))
  expect_true(is.na(rv(null)$first_run_at))

  # a therapeutic reduction registers under improvement = "decrease"
  dec <- scedTwoSDBand(base, c(2, 1, 3, 2), improvement = "decrease")
  expect_true(est(dec))
})

test_that("split-middle celeration matches the White & Haring hand computation", {
  # x = 1:8, mid-points (2.5, 3.5) and (6.5, 7.5) -> slope 1, intercept 1,
  # median residual 0 (no vertical shift needed).
  y <- c(2, 4, 3, 5, 6, 8, 7, 9)
  r <- rv(scedCelerationLine(y, method = "split_middle"))
  expect_equal(r$estimate, 1, tolerance = 1e-12)
  expect_equal(r$intercept, 1, tolerance = 1e-12)
  # the defining split-middle property: the median residual is zero
  expect_equal(stats::median(r$residuals), 0, tolerance = 1e-12)
})

test_that("split-middle applies the median-residual vertical shift", {
  # provisional line leaves a nonzero median residual; the adjustment zeroes it.
  y <- c(1, 1, 1, 5, 9, 9, 9)               # step up at the middle
  r <- rv(scedCelerationLine(y, method = "split_middle"))
  expect_equal(stats::median(r$residuals), 0, tolerance = 1e-12)
})

test_that("OLS celeration equals lm and honours an explicit time axis", {
  y <- c(2, 4, 3, 5, 6, 8, 7, 9); x <- (1:8) * 2
  r <- rv(scedCelerationLine(y, time = x, method = "ols"))
  expect_equal(r$estimate, unname(coef(stats::lm(y ~ x))[2]), tolerance = 1e-12)
  expect_equal(r$intercept, unname(coef(stats::lm(y ~ x))[1]), tolerance = 1e-12)
})

test_that("scedABAB separates four phases and computes both reversal contrasts", {
  df <- data.frame(
    value = c(10, 11, 9, 10, 18, 19, 20, 19, 10, 12, 11, 10, 21, 22, 20, 23),
    phase = rep(c("A1", "B1", "A2", "B2"), each = 4))
  res <- scedABAB(df)
  expect_s3_class(res, "sced_abab")
  expect_equal(names(res$contrasts), c("B1 vs A1", "B2 vs A2"))
  expect_equal(nrow(res$summary), 2L)
  expect_true(all(res$summary$band_flag))
  expect_true(all(res$summary$NAP > 0.9))
  # the B1-vs-A1 contrast really compares those two phases
  expect_equal(res$contrasts[["B1 vs A1"]]$baseline, "A1")
  expect_equal(res$contrasts[["B1 vs A1"]]$intervention, "B1")
})

test_that("scedABAB flags nothing on a null series", {
  df <- data.frame(
    value = c(10, 11, 9, 10, 10, 11, 9, 10, 10, 12, 11, 10, 10, 9, 11, 10),
    phase = rep(c("A1", "B1", "A2", "B2"), each = 4))
  res <- scedABAB(df)
  expect_false(any(res$summary$band_flag))
})

test_that("scedABAB reads a PhysioExperiment single-channel series", {
  skip_if_not_installed("SummarizedExperiment")
  vals <- c(10, 11, 9, 18, 19, 20, 10, 12, 21, 22, 20, 23)
  ph <- rep(c("A1", "B1", "A2", "B2"), each = 3)
  pe <- PhysioCore::PhysioExperiment(
    assays = S4Vectors::SimpleList(raw = matrix(vals, ncol = 1)),
    colData = S4Vectors::DataFrame(label = "ch1"),
    rowData = S4Vectors::DataFrame(phase = ph),
    samplingRate = 1)
  res <- scedABAB(pe, assay = "raw")
  expect_s3_class(res, "sced_abab")
  expect_equal(nrow(res$phase_data), length(vals))
  expect_equal(res$summary$NAP, scedABAB(
    data.frame(value = vals, phase = ph))$summary$NAP)
})

test_that("scedABAB errors informatively on missing phases", {
  df <- data.frame(value = 1:6, phase = rep(c("A1", "B1"), each = 3))
  expect_error(scedABAB(df), "missing")
})

test_that("plot.sced_abab returns a ggplot", {
  skip_if_not_installed("ggplot2")
  df <- data.frame(
    value = c(10, 11, 9, 10, 18, 19, 20, 19, 10, 12, 11, 10, 21, 22, 20, 23),
    phase = rep(c("A1", "B1", "A2", "B2"), each = 4))
  expect_s3_class(plot(scedABAB(df)), "ggplot")
})


# --- regression tests for adversarial-review findings ---

test_that("scedABAB orders observations by time (sequence-dependent stats)", {
  skip_if_not_installed("SingleCaseES")
  Aord <- c(20, 21, 22, 23, 24, 25); Bv <- c(30, 31, 32, 33)
  shuf <- c(3, 1, 5, 2, 6, 4)                      # scramble baseline rows...
  df <- rbind(
    data.frame(time = (1:6)[shuf], value = Aord[shuf], phase = "A1"),  # ...true order in time
    data.frame(time = 7:10,  value = Bv, phase = "B1"),
    data.frame(time = 11:16, value = Aord, phase = "A2"),
    data.frame(time = 17:20, value = Bv, phase = "B2"))
  res <- scedABAB(df)
  # the baseline-trend-corrected Tau-U must use time order, matching SingleCaseES
  expect_equal(res$summary$Tau_U[1],
               SingleCaseES::Tau_U(A_data = Aord, B_data = Bv)$Est,
               tolerance = 1e-9)
  # and that genuinely differs from the (wrong) row-order value
  expect_false(isTRUE(all.equal(
    PhysioCore::resultValue(scedTauU(Aord[shuf], Bv))$estimate,
    res$summary$Tau_U[1])))
})

test_that("factor value columns are read by label, not by integer code", {
  skip_if_not_installed("SingleCaseES")
  fa <- factor(c("20", "20", "26")); fb <- factor(c("28", "25", "30"))
  expect_equal(PhysioCore::resultValue(scedNAP(fa, fb))$estimate,
               SingleCaseES::NAP(c(20, 20, 26), c(28, 25, 30))$Est,
               tolerance = 1e-12)
  dff <- data.frame(
    value = factor(c("20", "20", "26", "28", "25", "30",
                     "22", "24", "31", "33", "20", "29")),
    phase = rep(c("A1", "B1", "A2", "B2"), each = 3))
  num <- data.frame(value = as.numeric(as.character(dff$value)), phase = dff$phase)
  expect_equal(scedABAB(dff)$summary$NAP, scedABAB(num)$summary$NAP)
  # genuinely non-numeric input errors instead of coercing to NA
  expect_error(scedNAP(factor(c("hi", "lo", "mid")), c(1, 2, 3)), "not numeric")
})
