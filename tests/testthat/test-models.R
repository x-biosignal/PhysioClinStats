library(testthat)
library(PhysioClinStats)

rv <- function(r) PhysioCore::resultValue(r)

test_that("fitMixedModel fixed effects match lme4::lmer on sleepstudy", {
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")
  fm <- fitMixedModel(sleepstudy, Reaction ~ Days + (Days | Subject))
  ref <- lme4::fixef(lme4::lmer(Reaction ~ Days + (Days | Subject),
                                data = sleepstudy))
  fx <- rv(fm)$fixed
  expect_equal(fx$estimate, unname(ref), tolerance = 1e-6)
  expect_s4_class(fm, "AnalysisResult")
  expect_equal(PhysioCore::resultType(fm), "mixed_model")
  expect_equal(rv(fm)$sigma, sigma(lme4::lmer(Reaction ~ Days + (Days | Subject),
                                              data = sleepstudy)),
               tolerance = 1e-6)
  expect_true(nrow(rv(fm)$random) >= 2)
})

test_that("fitMixedModel offers Satterthwaite and Kenward-Roger df with valid p", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  data(sleepstudy, package = "lme4")
  s <- rv(fitMixedModel(sleepstudy, Reaction ~ Days + (Days | Subject),
                        df = "satterthwaite"))$fixed
  k <- rv(fitMixedModel(sleepstudy, Reaction ~ Days + (Days | Subject),
                        df = "kenward-roger"))$fixed
  expect_true(all(is.finite(s$df)) && all(is.finite(s$p_value)))
  expect_true(all(is.finite(k$df)) && all(is.finite(k$p_value)))
})

test_that("fitMixedModel accepts a separate random-effects formula", {
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")
  a <- rv(fitMixedModel(sleepstudy, Reaction ~ Days, random = ~ (Days | Subject)))
  b <- rv(fitMixedModel(sleepstudy, Reaction ~ Days + (Days | Subject)))
  expect_equal(a$fixed$estimate, b$fixed$estimate, tolerance = 1e-8)
  expect_error(fitMixedModel(sleepstudy, Reaction ~ Days), "random")
})

test_that("fitMMRM matches mmrm::mmrm on the FEV1 example", {
  skip_if_not_installed("mmrm")
  data(fev_data, package = "mmrm")
  fit <- fitMMRM(fev_data, "FEV1", "ARMCD", "AVISIT", "USUBJID",
                 covariates = c("RACE", "SEX"), df = "kenward-roger")
  ref <- mmrm::mmrm(
    FEV1 ~ RACE + SEX + ARMCD * AVISIT + us(AVISIT | USUBJID),
    data = fev_data, method = "Kenward-Roger")
  cr <- summary(ref)$coefficients
  cm <- rv(fit)$coefficients
  cm <- cm[match(rownames(cr), cm$term), ]
  expect_equal(cm$estimate, unname(cr[, "Estimate"]), tolerance = 1e-4)
  expect_equal(cm$std_error, unname(cr[, "Std. Error"]), tolerance = 1e-4)
  expect_equal(cm$df, unname(cr[, "df"]), tolerance = 0.1)
  expect_equal(rv(fit)$backend, "mmrm")
})

test_that("fitMMRM Satterthwaite df is selectable and gives valid p-values", {
  skip_if_not_installed("mmrm")
  data(fev_data, package = "mmrm")
  fit <- fitMMRM(fev_data, "FEV1", "ARMCD", "AVISIT", "USUBJID",
                 covariates = c("RACE", "SEX"), df = "satterthwaite")
  co <- rv(fit)$coefficients
  expect_true(all(is.finite(co$df)))
  expect_true(all(is.finite(co$p_value) & co$p_value >= 0 & co$p_value <= 1))
})

test_that("fitMMRM falls back to nlme::gls when mmrm is unavailable", {
  skip_if_not_installed("nlme")
  set.seed(1); n <- 24
  d <- expand.grid(subject = factor(seq_len(n)), time = factor(1:4))
  d$y <- 30 + 2 * (d$time == "2") + 3 * (d$time == "3") + 4 * (d$time == "4") +
    rep(rnorm(n, 0, 3), times = 4) + rnorm(nrow(d), 0, 1)
  d$arm <- factor(rep(c("A", "B"), length.out = n)[as.integer(d$subject)])
  # shadow requireNamespace so the mmrm branch is skipped
  fit <- local({
    requireNamespace <- function(pkg, ...) if (identical(pkg, "mmrm")) FALSE
      else base::requireNamespace(pkg, ...)
    environment(fitMMRM) <- environment()
    fitMMRM(d, "y", "arm", "time", "subject")
  })
  expect_equal(rv(fit)$backend, "gls")
  co <- rv(fit)$coefficients
  expect_true(all(is.finite(co$df)))
  # recovers the simulated visit-4 effect (~4)
  expect_equal(co$estimate[co$term == "time4"], 4, tolerance = 1)
})

test_that("estimatedMarginalMeans matches emmeans on the mmrm fit", {
  skip_if_not_installed("mmrm")
  skip_if_not_installed("emmeans")
  data(fev_data, package = "mmrm")
  fit <- fitMMRM(fev_data, "FEV1", "ARMCD", "AVISIT", "USUBJID",
                 covariates = c("RACE", "SEX"))
  ref <- mmrm::mmrm(
    FEV1 ~ RACE + SEX + ARMCD * AVISIT + us(AVISIT | USUBJID),
    data = fev_data, method = "Kenward-Roger")
  emm <- estimatedMarginalMeans(fit, ~ ARMCD | AVISIT)
  eref <- as.data.frame(summary(emmeans::emmeans(ref, ~ ARMCD | AVISIT)))
  mine <- rv(emm)$emmeans
  expect_equal(mine$emmean, eref$emmean, tolerance = 1e-6)
  expect_equal(mine$SE, eref$SE, tolerance = 1e-6)
})

test_that("pairwiseContrasts matches emmeans pairwise differences", {
  skip_if_not_installed("mmrm")
  skip_if_not_installed("emmeans")
  data(fev_data, package = "mmrm")
  fit <- fitMMRM(fev_data, "FEV1", "ARMCD", "AVISIT", "USUBJID",
                 covariates = c("RACE", "SEX"))
  ref <- mmrm::mmrm(
    FEV1 ~ RACE + SEX + ARMCD * AVISIT + us(AVISIT | USUBJID),
    data = fev_data, method = "Kenward-Roger")
  pc <- rv(pairwiseContrasts(fit, ~ ARMCD))$contrasts
  pref <- as.data.frame(summary(
    emmeans::contrast(emmeans::emmeans(ref, ~ ARMCD), method = "pairwise")))
  expect_equal(pc$estimate[1], pref$estimate[1], tolerance = 1e-6)
})

test_that("estimatedMarginalMeans can also compute contrasts in one call", {
  skip_if_not_installed("mmrm")
  skip_if_not_installed("emmeans")
  data(fev_data, package = "mmrm")
  fit <- fitMMRM(fev_data, "FEV1", "ARMCD", "AVISIT", "USUBJID",
                 covariates = c("RACE", "SEX"))
  r <- rv(estimatedMarginalMeans(fit, ~ ARMCD, contrasts = "pairwise"))
  expect_false(is.null(r$contrasts))
  expect_true(nrow(r$emmeans) >= 2)
})

test_that("peToLong flattens 2D and 3D assays with a value round-trip", {
  skip_if_not_installed("SummarizedExperiment")
  pe2 <- PhysioCore::PhysioExperiment(
    assays = S4Vectors::SimpleList(raw = matrix(1:12, nrow = 4)),
    colData = S4Vectors::DataFrame(label = c("c1", "c2", "c3")),
    samplingRate = 1)
  l2 <- peToLong(pe2)
  expect_equal(nrow(l2), 12L)
  expect_setequal(names(l2), c("subject", "time", "channel", "value"))

  arr <- array(1:24, dim = c(4, 3, 2))
  pe3 <- PhysioCore::PhysioExperiment(
    assays = S4Vectors::SimpleList(raw = arr),
    colData = S4Vectors::DataFrame(label = c("c1", "c2", "c3")),
    samplingRate = 1)
  l3 <- peToLong(pe3, subjects = c("s1", "s2"))
  expect_equal(nrow(l3), 24L)
  expect_equal(length(unique(l3$subject)), 2L)
  expect_equal(
    l3$value[l3$time == 1 & l3$channel == "c1" & l3$subject == "s2"],
    arr[1, 1, 2])
})

test_that("the model wrappers validate their inputs", {
  skip_if_not_installed("mmrm")
  data(fev_data, package = "mmrm")
  expect_error(fitMMRM(fev_data, "NOPE", "ARMCD", "AVISIT", "USUBJID"),
               "not in 'data'")
  expect_error(estimatedMarginalMeans(
    PhysioCore::AnalysisResult("x", result = list()), ~ a),
    "no fitted model")
})


# --- regression tests for adversarial-review findings (gls fallback) ---

.fit_mmrm_gls <- function(...) {
  local({
    requireNamespace <- function(pkg, ...) if (identical(pkg, "mmrm")) FALSE
      else base::requireNamespace(pkg, ...)
    environment(fitMMRM) <- environment()
    suppressMessages(fitMMRM(...))
  })
}

test_that("the gls fallback fits homogeneous cs/ar1, not the heterogeneous variant", {
  skip_if_not_installed("mmrm")
  skip_if_not_installed("nlme")
  data(fev_data, package = "mmrm")
  d <- fev_data[stats::complete.cases(
    fev_data[, c("FEV1", "RACE", "SEX", "ARMCD", "AVISIT", "USUBJID")]), ]
  d$AVISIT <- droplevels(d$AVISIT); d$USUBJID <- droplevels(d$USUBJID)
  for (cv in c("compound-symmetry", "ar1")) {
    covf <- if (cv == "compound-symmetry") "cs" else "ar1"
    ref <- mmrm::mmrm(stats::reformulate(
      c("RACE", "SEX", "ARMCD*AVISIT", sprintf("%s(AVISIT | USUBJID)", covf)),
      response = "FEV1"), data = d, method = "Satterthwaite")
    g <- rv(.fit_mmrm_gls(d, "FEV1", "ARMCD", "AVISIT", "USUBJID",
                          covariates = c("RACE", "SEX"), covariance = cv))$coefficients
    rc <- summary(ref)$coefficients
    cm <- g$estimate[match(rownames(rc), g$term)]
    # matches the HOMOGENEOUS structure to which mmrm cs/ar1 map
    expect_equal(cm, unname(rc[, "Estimate"]), tolerance = 1e-3)
  }
})

test_that("the gls fallback uses between-within df for between-subject effects", {
  skip_if_not_installed("nlme")
  set.seed(3); n <- 40
  d <- expand.grid(subject = factor(seq_len(n)), time = factor(1:4))
  d$arm <- factor(rep(c("A", "B"), length.out = n)[as.integer(d$subject)])
  d$y <- 10 + 3 * (d$arm == "B") + 2 * (d$time == "4") +
    rep(rnorm(n, 0, 2), times = 4) + rnorm(nrow(d), 0, 1)
  co <- rv(.fit_mmrm_gls(d, "y", "arm", "time", "subject"))$coefficients
  arm_df <- co$df[co$term == "armB"]                 # between-subject effect
  # ~ n_subjects, and far below the anti-conservative residual df (N - p)
  expect_lt(arm_df, nrow(d) - nrow(co))
  expect_lte(arm_df, n)
  # p-value is recomputed from that df
  t_arm <- co$statistic[co$term == "armB"]
  expect_equal(co$p_value[co$term == "armB"],
               2 * stats::pt(-abs(t_arm), df = arm_df), tolerance = 1e-10)
})
