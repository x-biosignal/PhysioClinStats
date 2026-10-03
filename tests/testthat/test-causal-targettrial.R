library(testthat)
library(PhysioClinStats)

tt_rv <- function(x) PhysioCore::resultValue(x)

.tt_protocol <- function(follow_up = 2) {
  targetTrialProtocol(
    eligibility = function(data) data$eligible,
    treatment_strategies = list(never = 0, always = 1),
    assignment = "clone at eligibility",
    time_zero = 0,
    follow_up = follow_up,
    outcome = "event or recovery by horizon",
    causal_contrast = "always versus never",
    analysis_plan = "clone-censor-weight"
  )
}

.tt_data <- function(n = 120L, survival = TRUE, zero_events = FALSE) {
  set.seed(811)
  id <- seq_len(n)
  x <- stats::rnorm(n)
  a0 <- stats::rbinom(n, 1, stats::plogis(-0.2 + 0.5 * x))
  long <- do.call(rbind, lapply(id, function(i) {
    a1 <- stats::rbinom(1, 1, stats::plogis(-0.3 + 0.7 * a0[i] + 0.3 * x[i]))
    a2 <- stats::rbinom(1, 1, stats::plogis(-0.2 + 0.8 * a1 + 0.3 * x[i]))
    data.frame(
      id = i, visit = 0:2, a = c(a0[i], a1, a2),
      x = x[i], l = c(0, a0[i] + x[i], a1 + x[i]),
      eligible = TRUE, stringsAsFactors = FALSE
    )
  }))
  if (survival) {
    base_a <- a0[long$id]
    event_subject <- if (zero_events) rep(0L, n) else
      stats::rbinom(n, 1, stats::plogis(-1 + 0.9 * a0 - 0.3 * x))
    long$event_time <- 2
    long$event <- event_subject[long$id]
  } else {
    recovery <- stats::rbinom(
      n, 1, stats::plogis(-0.4 + 0.8 * a0 + 0.2 * x))
    long$recovery <- recovery[long$id]
  }
  long
}

test_that("targetTrialProtocol validates and reports static and dynamic rules", {
  protocol <- .tt_protocol()
  expect_s3_class(protocol, "target_trial_protocol")
  expect_true(protocol$portable)
  expect_match(capture.output(print(protocol))[1], "Target-trial protocol")
  expect_true(all(protocol$validation$valid))
  expect_match(protocol$protocol_id, "^ttp-")

  dynamic <- targetTrialProtocol(
    eligibility = ~ eligible,
    treatment_strategies = list(
      never = 0,
      adaptive = function(data, history) as.integer(tail(history$l, 1) > 0)),
    assignment = "clone", time_zero = 0, follow_up = 2,
    outcome = "event", causal_contrast = "adaptive versus never",
    analysis_plan = "clone-censor-weight")
  expect_false(dynamic$portable)
  expect_equal(dynamic$strategy_specs$adaptive$type, "dynamic")
  expect_true(nzchar(dynamic$strategy_specs$adaptive$hash))

  expect_error(targetTrialProtocol(
    TRUE, list(a = 0), "x", 0, 2, "y", "z", "plan"),
  "at least two")
  expect_error(targetTrialProtocol(
    TRUE, list(a = 0, b = function(x) 1), "x", 0, 2, "y", "z", "plan"),
  "data, history")
  expect_error(targetTrialProtocol(
    TRUE, list(a = 0, b = 1), "x", 0, 2, "y", "z", "plan",
    version = "one"), "semantic")
})

test_that("dynamic strategies receive no outcome or current treatment leakage", {
  data <- .tt_data(n = 80L)
  guarded_strategy <- function(data, history) {
    stopifnot(!"event_time" %in% names(data))
    stopifnot(!"event" %in% names(history))
    stopifnot(is.na(tail(history$a, 1L)))
    as.integer(data$l > 0)
  }
  protocol <- targetTrialProtocol(
    eligibility = function(data) data$eligible,
    treatment_strategies = list(never = 0, guarded = guarded_strategy),
    assignment = "clone", time_zero = 0, follow_up = 2,
    outcome = "event", causal_contrast = "guarded versus never",
    analysis_plan = "clone-censor-weight")
  expect_no_error(targetTrialEmulate(
    data, protocol, "id", "visit", "a", "event_time", event = "event",
    baseline_covariates = "x", time_varying_covariates = "l",
    estimand = "per_protocol", weight_truncation = c(0, 1)))
})

test_that("internal baseline IPTW follows ATE, ATT and stabilization formulas", {
  d <- data.frame(x = c(-2, -1, 0, 1, 2, 3))
  a <- c(0, 1, 0, 1, 0, 1)
  internal <- PhysioClinStats:::.baseline_iptw(
    a, d, estimand = "ATE", stabilized = FALSE)
  p <- internal$propensity
  expect_equal(
    internal$weights,
    ifelse(a == 1, 1 / p, 1 / (1 - p)), tolerance = 1e-12)
  stabilized <- PhysioClinStats:::.baseline_iptw(
    a, d, estimand = "ATE", stabilized = TRUE)
  expect_equal(
    stabilized$weights,
    internal$weights * ifelse(a == 1, mean(a), 1 - mean(a)),
    tolerance = 1e-12)
  att <- PhysioClinStats:::.baseline_iptw(
    a, d, estimand = "ATT", stabilized = FALSE)
  expect_equal(att$weights, ifelse(a == 1, 1, p / (1 - p)),
               tolerance = 1e-12)
  expect_true(all(c("propensity", "weight_quantiles",
                    "effective_sample_size", "balance") %in%
                  names(internal$diagnostics)))
})

test_that("baseline ATE weights match WeightIt and ipw backends", {
  set.seed(17)
  d <- data.frame(x = stats::rnorm(200), z = stats::rbinom(200, 1, 0.4))
  a <- stats::rbinom(200, 1, stats::plogis(-0.1 + 0.5 * d$x - 0.3 * d$z))
  internal <- PhysioClinStats:::.baseline_iptw(
    a, d, estimand = "ATE", stabilized = FALSE, backend = "internal")
  if (requireNamespace("WeightIt", quietly = TRUE)) {
    weightit <- PhysioClinStats:::.baseline_iptw(
      a, d, estimand = "ATE", stabilized = FALSE, backend = "WeightIt")
    expect_equal(weightit$weights, internal$weights, tolerance = 1e-6)
  }
  if (requireNamespace("ipw", quietly = TRUE)) {
    ipw <- PhysioClinStats:::.baseline_iptw(
      a, d, estimand = "ATE", stabilized = FALSE, backend = "ipw")
    expect_equal(ipw$weights, internal$weights, tolerance = 1e-6)
  }
})

test_that("pinned baseline fixture reproduces all IPTW backend weights", {
  fixture <- readRDS(test_path(
    "fixtures", "causal", "baseline-iptw.rds"))
  covariates <- fixture$data[, c("x", "z")]
  for (backend in c("internal", "WeightIt", "ipw")) {
    if (backend != "internal" &&
        !requireNamespace(backend, quietly = TRUE)) {
      skip(sprintf("%s is not installed", backend))
    }
    result <- PhysioClinStats:::.baseline_iptw(
      fixture$data$treatment, covariates,
      estimand = "ATE", stabilized = FALSE, backend = backend)
    expect_equal(result$weights, fixture$weights[[backend]],
                 tolerance = 1e-12)
    expect_equal(result$propensity, fixture$propensity,
                 tolerance = 1e-12)
  }
})

test_that("hand-worked cloning censors at first deviation without future rows", {
  data <- data.frame(
    id = c(1, 1, 1, 2, 2, 2), visit = rep(0:2, 2),
    a = c(0, 0, 1, 1, 1, 1), ytime = 2, event = 0,
    eligible = TRUE, x = 0, .source_row = 1:6)
  protocol <- .tt_protocol()
  columns <- list(
    id = "id", time = "visit", treatment = "a",
    outcome = "ytime", event = "event",
    baseline = "x", time_varying = character())
  probability <- list(
    denominator_probability = rep(0.5, 6),
    numerator_probability = rep(0.4, 6))
  clones <- PhysioClinStats:::.tt_clone(
    data, protocol, columns, probability, horizon = 2)
  expect_equal(length(unique(clones$clone_id)), 4L)
  expect_equal(nrow(clones), 8L)
  never_1 <- clones[clones$clone_id == "1::never", ]
  expect_equal(never_1$adherent, c(TRUE, TRUE, FALSE))
  expect_equal(never_1$raw_weight, c(0.8, 0.64, 0.64),
               tolerance = 1e-12)
  expect_equal(tail(never_1$censor_reason, 1), "artificial_deviation")
  always_1 <- clones[clones$clone_id == "1::always", ]
  expect_equal(nrow(always_1), 1L)
  expect_true(always_1$artificial_censor)
  expect_equal(always_1$raw_weight, 1)
  expect_false(any(clones$.source_row[clones$clone_id == "1::always"] > 1))
})

test_that("hand-auditable target-trial fixture matches every clone and Cox field", {
  skip_if_not_installed("survival")
  fixture <- readRDS(test_path(
    "fixtures", "causal", "target-trial-reference.rds"))
  data <- fixture$input
  data$.source_row <- seq_len(nrow(data))
  columns <- list(
    id = "id", time = "visit", treatment = "treatment",
    outcome = "outcome_time", event = "event",
    baseline = "x", time_varying = character())
  probability <- list(
    denominator_probability = rep(0.5, nrow(data)),
    numerator_probability = rep(0.4, nrow(data)))
  clones <- PhysioClinStats:::.tt_clone(
    data, .tt_protocol(), columns, probability, horizon = 2)
  expected <- fixture$expected_clones
  expected$censor_reason[expected$censor_reason == ""] <- NA_character_
  actual <- data.frame(
    clone_id = clones$clone_id,
    original_id = as.integer(clones$original_id),
    strategy = clones$strategy,
    visit = clones$visit,
    observed_treatment = clones$treatment,
    expected_treatment = clones$expected_treatment,
    adherent = clones$adherent,
    artificial_censor = clones$artificial_censor,
    censor_reason = clones$censor_reason,
    denominator_probability = clones$denominator_probability,
    numerator_probability = clones$numerator_probability,
    interval_weight = clones$interval_weight,
    raw_weight = clones$raw_weight,
    stringsAsFactors = FALSE)
  expect_equal(actual, expected, tolerance = 1e-12)

  clones$analysis_weight <- clones$raw_weight
  endpoint <- PhysioClinStats:::.tt_clone_endpoint(
    clones, columns, horizon = 2)
  expected_endpoint <- fixture$expected_endpoint
  endpoint <- endpoint[, names(expected_endpoint)]
  endpoint$original_id <- as.integer(endpoint$original_id)
  endpoint$strategy <- as.character(endpoint$strategy)
  expected_endpoint$strategy <- as.character(expected_endpoint$strategy)
  endpoint <- endpoint[order(endpoint$clone_id), , drop = FALSE]
  expected_endpoint <- expected_endpoint[
    order(expected_endpoint$clone_id), , drop = FALSE]
  row.names(endpoint) <- NULL
  row.names(expected_endpoint) <- NULL
  expect_equal(endpoint, expected_endpoint, tolerance = 1e-12)

  endpoint$strategy <- factor(
    endpoint$strategy, levels = c("never", "always"))
  direct <- survival::coxph(
    survival::Surv(endpoint_time, endpoint_event) ~ strategy,
    data = endpoint, weights = analysis_weight,
    robust = TRUE, cluster = original_id, model = TRUE)
  summary <- summary(direct, conf.int = 0.95)
  expect_equal(unname(stats::coef(direct)),
               fixture$survival$coefficient, tolerance = 1e-12)
  expect_equal(unname(exp(stats::coef(direct))),
               fixture$survival$hazard_ratio, tolerance = 1e-12)
  expect_equal(unname(summary$coefficients[, "robust se"]),
               fixture$survival$robust_se, tolerance = 1e-12)
})

test_that("per-protocol output retains raw weights and matches direct robust Cox", {
  skip_if_not_installed("survival")
  data <- .tt_data()
  result <- targetTrialEmulate(
    data, .tt_protocol(), id = "id", time = "visit",
    treatment = "a", outcome = "event_time", event = "event",
    baseline_covariates = "x", time_varying_covariates = "l",
    estimand = "per_protocol", stabilized = TRUE,
    weight_truncation = c(0, 1))
  value <- tt_rv(result)
  expect_false(is.null(value$clone_data))
  expect_true(all(c(
    "raw_weight", "analysis_weight", "artificial_censor",
    "censor_reason") %in% names(value$clone_data)))
  expect_equal(value$clone_data$raw_weight,
               value$clone_data$analysis_weight, tolerance = 1e-12)
  expect_gt(value$censoring_counts$count[
    value$censoring_counts$reason == "artificial_deviation"], 0)

  endpoint <- value$endpoint_data
  endpoint$strategy_factor <- factor(
    endpoint$strategy, levels = c("never", "always"))
  direct <- survival::coxph(
    survival::Surv(endpoint_time, endpoint_event) ~ strategy_factor,
    data = endpoint, weights = analysis_weight,
    robust = TRUE, cluster = original_id, model = TRUE)
  summary <- summary(direct, conf.int = 0.95)
  effect <- value$effects
  expect_equal(effect$estimate,
               as.numeric(summary$coefficients[, "exp(coef)"]),
               tolerance = 1e-10)
  expect_equal(effect$std_error,
               as.numeric(summary$coefficients[, "robust se"]),
               tolerance = 1e-10)
  expect_equal(effect$conf_low,
               as.numeric(summary$conf.int[, "lower .95"]),
               tolerance = 1e-10)

  ninety <- targetTrialEmulate(
    data, .tt_protocol(), id = "id", time = "visit",
    treatment = "a", outcome = "event_time", event = "event",
    baseline_covariates = "x", time_varying_covariates = "l",
    estimand = "per_protocol", stabilized = TRUE,
    weight_truncation = c(0, 1), conf.level = 0.9)
  ninety_endpoint <- tt_rv(ninety)$endpoint_data
  ninety_endpoint$strategy_factor <- factor(
    ninety_endpoint$strategy, levels = c("never", "always"))
  ninety_direct <- summary(survival::coxph(
    survival::Surv(endpoint_time, endpoint_event) ~ strategy_factor,
    data = ninety_endpoint, weights = analysis_weight,
    robust = TRUE, cluster = original_id), conf.int = 0.9)
  expect_equal(tt_rv(ninety)$effects$conf_low,
               as.numeric(ninety_direct$conf.int[, 3L]),
               tolerance = 1e-10)
  expect_equal(tt_rv(ninety)$effects$conf_high,
               as.numeric(ninety_direct$conf.int[, 4L]),
               tolerance = 1e-10)
})

test_that("ITT does not clone or artificially censor", {
  skip_if_not_installed("survival")
  data <- .tt_data()
  result <- targetTrialEmulate(
    data, .tt_protocol(), "id", "visit", "a",
    "event_time", event = "event",
    baseline_covariates = "x",
    estimand = "intention_to_treat",
    weight_truncation = c(0, 1))
  value <- tt_rv(result)
  expect_null(value$clone_data)
  expect_equal(nrow(value$endpoint_data), length(unique(data$id)))
  expect_equal(value$censoring_counts$count[
    value$censoring_counts$reason == "artificial_deviation"], 0L)
})

test_that("swapping strategy order reverses the hazard-ratio orientation", {
  skip_if_not_installed("survival")
  data <- .tt_data()
  forward <- targetTrialEmulate(
    data, .tt_protocol(), "id", "visit", "a",
    "event_time", event = "event",
    baseline_covariates = "x", estimand = "intention_to_treat",
    weight_truncation = c(0, 1))
  reversed_protocol <- targetTrialProtocol(
    eligibility = function(data) data$eligible,
    treatment_strategies = list(always = 1, never = 0),
    assignment = "baseline assignment", time_zero = 0, follow_up = 2,
    outcome = "event", causal_contrast = "never versus always",
    analysis_plan = "intention-to-treat Cox")
  reversed <- targetTrialEmulate(
    data, reversed_protocol, "id", "visit", "a",
    "event_time", event = "event",
    baseline_covariates = "x", estimand = "intention_to_treat",
    weight_truncation = c(0, 1))
  expect_equal(
    tt_rv(forward)$effects$estimate *
      tt_rv(reversed)$effects$estimate,
    1, tolerance = 1e-10)
  expect_equal(tt_rv(forward)$effects$std_error,
               tt_rv(reversed)$effects$std_error,
               tolerance = 1e-10)
})

test_that("binary target trials return standardized risks, RD and RR", {
  data <- .tt_data(survival = FALSE)
  result <- targetTrialEmulate(
    data, .tt_protocol(), "id", "visit", "a", "recovery",
    baseline_covariates = "x",
    estimand = "intention_to_treat",
    weight_truncation = c(0, 1))
  value <- tt_rv(result)
  expect_setequal(value$effects$measure,
                  c("risk_difference", "risk_ratio"))
  expect_equal(names(value$standardized_risks), c("never", "always"))
  expect_true(all(value$standardized_risks > 0 &
                    value$standardized_risks < 1))
  expect_false(is.null(value$robust_vcov))
})

test_that("truncation changes only analysis weights and records exact rows", {
  data <- .tt_data()
  result <- targetTrialEmulate(
    data, .tt_protocol(), "id", "visit", "a",
    "event_time", event = "event",
    baseline_covariates = "x", time_varying_covariates = "l",
    weight_truncation = c(0.2, 0.8), max_weight = 2)
  clones <- tt_rv(result)$clone_data
  changed <- clones$weight_truncated
  expect_gt(sum(changed), 0)
  expect_true(all(clones$analysis_weight[changed] !=
                    clones$raw_weight[changed]))
  expect_equal(
    tt_rv(result)$positivity$truncation$truncation_count,
    sum(changed))
  valid <- clones$adherent & is.finite(clones$raw_weight)
  expect_equal(
    unname(tt_rv(result)$positivity$truncation$thresholds),
    unname(stats::quantile(clones$raw_weight[valid], c(0.2, 0.8))),
    tolerance = 1e-12)
  expect_true(all(clones$raw_weight[changed] > 0))
})

test_that("ITT retains raw weights and applies declared truncation separately", {
  data <- .tt_data()
  result <- targetTrialEmulate(
    data, .tt_protocol(), "id", "visit", "a",
    "event_time", event = "event",
    baseline_covariates = "x",
    estimand = "intention_to_treat",
    weight_truncation = c(0.2, 0.8), max_weight = 1.8)
  value <- tt_rv(result)
  expect_equal(value$raw_weights, value$endpoint_data$raw_weight)
  expect_equal(value$analysis_weights,
               value$endpoint_data$analysis_weight)
  expect_gt(sum(value$raw_weights != value$analysis_weights), 0)
  expect_equal(
    value$positivity$truncation$truncation_count,
    sum(value$raw_weights != value$analysis_weights))
})

test_that("subject-block permutation and treatment factor order are stable", {
  data <- .tt_data()
  first <- targetTrialEmulate(
    data, .tt_protocol(), "id", "visit", "a",
    "event_time", event = "event",
    baseline_covariates = "x", estimand = "intention_to_treat")
  blocks <- split(data, data$id)
  permuted <- do.call(rbind, rev(blocks))
  row.names(permuted) <- NULL
  permuted$a <- factor(permuted$a, levels = c(1, 0))
  second <- targetTrialEmulate(
    permuted, .tt_protocol(), "id", "visit", "a",
    "event_time", event = "event",
    baseline_covariates = "x", estimand = "intention_to_treat")
  expect_equal(tt_rv(first)$effects, tt_rv(second)$effects,
               tolerance = 1e-10)
})

test_that("zero-event target trial returns structured NA effects", {
  skip_if_not_installed("survival")
  data <- .tt_data(zero_events = TRUE)
  expect_warning(result <- targetTrialEmulate(
    data, .tt_protocol(), "id", "visit", "a",
    "event_time", event = "event",
    baseline_covariates = "x", estimand = "intention_to_treat"),
  "0 events")
  expect_true(all(is.na(tt_rv(result)$effects$estimate)))
  expect_match(tt_rv(result)$warnings, "zero events")
})

test_that("target-trial input defects fail before modelling", {
  data <- .tt_data()
  protocol <- .tt_protocol()
  duplicate <- rbind(data, data[1, ])
  expect_error(targetTrialEmulate(
    duplicate, protocol, "id", "visit", "a", "event_time", "event"),
  "id,time")
  decreasing <- data
  decreasing[c(1, 2), ] <- decreasing[c(2, 1), ]
  expect_error(targetTrialEmulate(
    decreasing, protocol, "id", "visit", "a", "event_time", "event"),
  "strictly increasing")
  missing_zero <- data[data$id != 1 | data$visit != 0, ]
  expect_error(targetTrialEmulate(
    missing_zero, protocol, "id", "visit", "a", "event_time", "event"),
  "time-zero")
  post_outcome <- data
  post_outcome$event_time[post_outcome$id == 1] <- 1
  expect_error(targetTrialEmulate(
    post_outcome, protocol, "id", "visit", "a", "event_time", "event"),
  "Post-outcome")
  missing_value <- data
  missing_value$x[1] <- NA
  expect_error(targetTrialEmulate(
    missing_value, protocol, "id", "visit", "a", "event_time", "event",
    baseline_covariates = "x"), "complete")
  missing_numerator <- data
  missing_numerator$numerator <- 1
  missing_numerator$numerator[2] <- NA
  expect_error(targetTrialEmulate(
    missing_numerator, protocol, "id", "visit", "a",
    "event_time", "event", baseline_covariates = "x",
    numerator_covariates = "numerator"), "Numerator-model")
  impossible <- targetTrialProtocol(
    function(data) data$eligible, list(never = 0, impossible = 2),
    "clone", 0, 2, "event", "contrast", "plan")
  expect_error(targetTrialEmulate(
    data, impossible, "id", "visit", "a", "event_time", "event"),
  "absent from observed")
})
