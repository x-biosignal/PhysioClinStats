library(testthat)
library(PhysioClinStats)

rv <- function(r) PhysioCore::resultValue(r)

test_that("survivalFit KM estimates and median match survival::survfit on lung", {
  skip_if_not_installed("survival")
  km <- survivalFit(survival::lung, "time", "status")
  ref <- survival::survfit(survival::Surv(time, status) ~ 1,
                           data = survival::lung)
  expect_equal(km$fit$surv, ref$surv, tolerance = 1e-8)
  expect_equal(km$fit$time, ref$time)
  expect_equal(km$quantiles$median, 310, tolerance = 1e-8)
  expect_equal(km$quantiles$median_lower, 285, tolerance = 1e-8)
  expect_equal(km$quantiles$median_upper, 363, tolerance = 1e-8)
})

test_that("survivalFit stratifies and matches survfit medians per group", {
  skip_if_not_installed("survival")
  km <- survivalFit(survival::lung, "time", "status", group = "sex")
  expect_equal(sort(km$quantiles$median), c(270, 426), tolerance = 1e-8)
  ref <- survival::survfit(survival::Surv(time, status) ~ sex,
                           data = survival::lung)
  expect_equal(km$fit$surv, ref$surv, tolerance = 1e-8)
})

test_that("survivalFit matches survfit on the veteran dataset too", {
  skip_if_not_installed("survival")
  km <- survivalFit(survival::veteran, "time", "status")
  ref <- survival::survfit(survival::Surv(time, status) ~ 1,
                           data = survival::veteran)
  expect_equal(km$fit$surv, ref$surv, tolerance = 1e-8)
  expect_equal(unname(summary(ref)$table["median"]), km$quantiles$median,
               tolerance = 1e-8)
})

test_that("coxModel HR, SE and CI match survival::coxph on lung", {
  skip_if_not_installed("survival")
  cx <- coxModel(survival::lung, time = "time", event = "status",
                 covariates = c("age", "sex", "ph.ecog"))
  ref <- summary(survival::coxph(
    survival::Surv(time, status) ~ age + sex + ph.ecog, data = survival::lung))
  co <- rv(cx)$coefficients
  expect_equal(co$hr, unname(ref$coefficients[, "exp(coef)"]), tolerance = 1e-6)
  expect_equal(co$std_error, unname(ref$coefficients[, "se(coef)"]),
               tolerance = 1e-6)
  expect_equal(co$hr_lower, unname(ref$conf.int[, 3]), tolerance = 1e-6)
  expect_equal(co$hr_upper, unname(ref$conf.int[, 4]), tolerance = 1e-6)
  expect_s4_class(cx, "AnalysisResult")
  expect_equal(PhysioCore::resultType(cx), "cox_model")
  expect_false(rv(cx)$ph_violated)
})

test_that("coxModel surfaces the proportional-hazards violation on veteran", {
  skip_if_not_installed("survival")
  expect_warning(
    cx <- coxModel(survival::veteran, time = "time", event = "status",
                   covariates = c("trt", "karno", "celltype")),
    "Proportional-hazards")
  cx <- suppressWarnings(coxModel(
    survival::veteran, time = "time", event = "status",
    covariates = c("trt", "karno", "celltype")))
  expect_true(rv(cx)$ph_violated)
  ref <- survival::cox.zph(survival::coxph(
    survival::Surv(time, status) ~ trt + karno + celltype,
    data = survival::veteran))
  expect_equal(rv(cx)$ph_test["GLOBAL", "p"],
               unname(ref$table["GLOBAL", "p"]), tolerance = 1e-8)
})

test_that("coxModel accepts an explicit formula", {
  skip_if_not_installed("survival")
  cx <- coxModel(survival::lung,
                 formula = survival::Surv(time, status) ~ age + sex)
  ref <- summary(survival::coxph(survival::Surv(time, status) ~ age + sex,
                                 data = survival::lung))
  expect_equal(rv(cx)$coefficients$hr, unname(ref$coefficients[, "exp(coef)"]),
               tolerance = 1e-6)
  expect_error(coxModel(survival::lung), "formula|covariates")
})

test_that("timeToMilestone right-censors non-attainers and marks first crossing", {
  d <- data.frame(
    id = rep(c("a", "b", "c"), each = 4), time = rep(c(0, 4, 8, 12), 3),
    value = c(10, 30, 55, 60,       # a crosses at t = 8
              12, 20, 28, 35,       # b never reaches 50
              5, 52, 40, 48))       # c crosses at t = 4, then dips back
  m <- timeToMilestone(d, "id", "time", "value", threshold = 50)
  expect_equal(m$event[m$id == "a"], 1L)
  expect_equal(m$time[m$id == "a"], 8)
  expect_equal(m$event[m$id == "b"], 0L)          # censored
  expect_equal(m$time[m$id == "b"], 12)           # at last observation
  expect_equal(m$time[m$id == "c"], 4)            # FIRST crossing, not the dip
  # the Surv attribute is only attached when survival is installed
  if (requireNamespace("survival", quietly = TRUE)) {
    expect_s3_class(attr(m, "surv"), "Surv")
  }
})

test_that("timeToMilestone honours direction and censor_at", {
  d <- data.frame(id = c("x", "x", "y", "y"), time = c(1, 2, 1, 2),
                  value = c(100, 40, 90, 80))
  m <- timeToMilestone(d, "id", "time", "value", threshold = 50,
                       direction = "decrease")
  expect_equal(m$event[m$id == "x"], 1L)          # 40 <= 50
  expect_equal(m$event[m$id == "y"], 0L)          # never <= 50
  m2 <- timeToMilestone(d, "id", "time", "value", threshold = 50,
                        direction = "decrease", censor_at = 99)
  expect_equal(m2$time[m2$id == "y"], 99)
})

test_that("milestoneHazard fits a Cox HR for the grouped milestone endpoint", {
  skip_if_not_installed("survival")
  set.seed(2); n <- 60
  traj <- do.call(rbind, lapply(seq_len(n), function(i) {
    grp <- if (i <= n / 2) "trt" else "ctrl"
    rate <- if (grp == "trt") 8 else 4
    data.frame(id = i, group = grp, time = c(0, 4, 8, 12, 16),
               value = cumsum(c(0, stats::rexp(4, 1 / rate))))
  }))
  ms <- timeToMilestone(traj, "id", "time", "value", threshold = 20)
  g <- unique(traj[, c("id", "group")])
  ms$group <- g$group[match(ms$id, g$id)]
  mh <- suppressWarnings(milestoneHazard(ms, "group"))
  expect_equal(PhysioCore::resultType(mh), "milestone_hazard")
  co <- rv(mh)$coefficients
  expect_true(all(c("hr", "hr_lower", "hr_upper") %in% names(co)))
  # treatment attains the milestone faster -> HR > 1 for the treated group
  expect_gt(co$hr[1], 1)
})

test_that("survivalFit computes a numbers-at-risk table matching survfit", {
  skip_if_not_installed("survival")
  km <- survivalFit(survival::lung, "time", "status", group = "sex")
  expect_true(all(c("time", "group", "n_risk") %in% names(km$risk_table)))
  t1 <- km$risk_table$time[km$risk_table$time > 0][1]
  ref <- summary(survival::survfit(survival::Surv(time, status) ~ sex,
                                   data = survival::lung),
                 times = t1, extend = TRUE)
  expect_equal(sort(km$risk_table$n_risk[km$risk_table$time == t1]),
               sort(ref$n.risk))
})

test_that("plot.physio_km returns a Kaplan-Meier ggplot with optional risk table", {
  skip_if_not_installed("survival")
  skip_if_not_installed("ggplot2")
  km <- survivalFit(survival::lung, "time", "status", group = "sex")
  expect_s3_class(plot(km, risk_table = FALSE), "ggplot")
  # with a risk table + patchwork, the composed object is still a ggplot/patchwork
  p <- plot(km, risk_table = TRUE)
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})

test_that("the survival wrappers validate their inputs", {
  skip_if_not_installed("survival")
  expect_error(survivalFit(survival::lung, "nope", "status"), "not")
  expect_error(milestoneHazard(data.frame(time = 1, event = 1), "grp"), "group")
})


# --- regression tests for adversarial-review findings ---

test_that("coxModel reports the robust SE that drives z, p and the CI", {
  skip_if_not_installed("survival")
  d <- survival::lung
  d$inst[is.na(d$inst)] <- 1
  cx <- coxModel(d, formula = survival::Surv(time, status) ~ age + sex +
                   survival::cluster(inst))
  co <- rv(cx)$coefficients
  # the reported SE must reproduce the reported z and the HR CI
  expect_equal(co$z, co$coef / co$std_error, tolerance = 1e-6)
  expect_equal(co$hr_lower,
               exp(co$coef - stats::qnorm(0.975) * co$std_error),
               tolerance = 1e-5)
  ref <- summary(survival::coxph(survival::Surv(time, status) ~ age + sex +
                                   survival::cluster(inst), data = d))
  expect_equal(co$std_error, unname(ref$coefficients[, "robust se"]),
               tolerance = 1e-8)
})

test_that("coxModel tolerates a zero-event (no-attainer) endpoint", {
  skip_if_not_installed("survival")
  ms <- data.frame(time = c(8, 8, 8, 8), event = c(0L, 0L, 0L, 0L),
                   group = c("t", "t", "c", "c"))
  expect_warning(res <- milestoneHazard(ms, "group"), "0 events")
  expect_equal(PhysioCore::resultType(res), "milestone_hazard")
  expect_false(rv(res)$ph_violated)
  # no crash, and the hazard ratio is NA rather than a spurious number
  expect_true(all(is.na(rv(res)$coefficients$hr)))
})
