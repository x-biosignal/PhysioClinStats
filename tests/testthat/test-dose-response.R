# Therapy dose-response model.

test_that("linear doseResponse recovers a planted dose slope", {
  set.seed(1)
  d <- data.frame(dose = rep(c(0, 5, 10, 20), each = 40))
  d$change <- 0.4 * d$dose + rnorm(nrow(d), 0, 1)
  dr <- doseResponse(d, outcome = "change", dose = "dose", form = "linear")
  expect_s4_class(dr, "AnalysisResult")
  expect_equal(PhysioCore::resultType(dr), "dose_response")
  co <- PhysioCore::resultValue(dr)$coefficients
  slope <- co$estimate[co$term == "dose"]
  expect_equal(slope, 0.4, tolerance = 0.1)
  expect_equal(unname(dr@estimate[["dose"]]), slope)
  # predicted curve spans the observed dose range
  crv <- PhysioCore::resultValue(dr)$curve
  expect_equal(range(crv$dose), c(0, 20))
  expect_equal(nrow(crv), 50L)
})

test_that("covariate/baseline adjustment is applied (ANCOVA)", {
  set.seed(2)
  d <- data.frame(dose = rep(c(0, 10, 20), each = 50),
                  base = rnorm(150, 50, 5))
  d$post <- 0.3 * d$dose + 0.8 * d$base + rnorm(150, 0, 1)
  dr <- doseResponse(d, outcome = "post", dose = "dose", baseline = "base")
  co <- PhysioCore::resultValue(dr)$coefficients
  expect_equal(co$estimate[co$term == "dose"], 0.3, tolerance = 0.12)
  expect_true("base" %in% co$term)                    # baseline entered as covariate
})

test_that("emax form recovers ED50 / Emax", {
  set.seed(3)
  d <- data.frame(dose = rep(c(0, 2, 5, 10, 20, 40), each = 20))
  d$y <- 2 + 8 * d$dose / (5 + d$dose) + rnorm(nrow(d), 0, 0.5)
  dr <- doseResponse(d, outcome = "y", dose = "dose", form = "emax")
  est <- dr@estimate
  expect_equal(unname(est[["Emax"]]), 8, tolerance = 1.5)
  expect_equal(unname(est[["ED50"]]), 5, tolerance = 2)
  expect_equal(unname(est[["E0"]]), 2, tolerance = 1)
})

test_that("log and spline forms run and produce a monotone-ish curve", {
  set.seed(4)
  d <- data.frame(dose = rep(c(1, 3, 6, 12, 24), each = 20))
  d$y <- 3 * log(d$dose) + rnorm(nrow(d), 0, 1)
  lg <- doseResponse(d, outcome = "y", dose = "dose", form = "log")
  expect_gt(PhysioCore::resultValue(lg)$curve$predicted[50],
            PhysioCore::resultValue(lg)$curve$predicted[1])   # rises with dose
  sp <- doseResponse(d, outcome = "y", dose = "dose", form = "spline", df = 3)
  expect_equal(PhysioCore::resultType(sp), "dose_response")
  expect_equal(nrow(PhysioCore::resultValue(sp)$coefficients), 4L) # intercept + 3 ns
})

test_that("subject random intercepts recover the slope (lme4)", {
  skip_if_not_installed("lme4")
  set.seed(5)
  ns <- 20
  subj <- rep(seq_len(ns), each = 4)
  ri <- rnorm(ns, 0, 3)[subj]
  dose <- rep(c(0, 5, 10, 20), times = ns)
  y <- 0.5 * dose + ri + rnorm(length(dose), 0, 1)
  d <- data.frame(subject = factor(subj), dose = dose, y = y)
  dr <- doseResponse(d, outcome = "y", dose = "dose", subject = "subject")
  co <- PhysioCore::resultValue(dr)$coefficients
  expect_equal(co$estimate[co$term == "dose"], 0.5, tolerance = 0.1)
})
