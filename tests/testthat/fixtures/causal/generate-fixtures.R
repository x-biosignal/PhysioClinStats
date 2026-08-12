fixture_dir <- if (dir.exists("tests/testthat/fixtures/causal")) {
  "tests/testthat/fixtures/causal"
} else {
  "physio-ecosystem/PhysioClinStats/tests/testthat/fixtures/causal"
}

for (package in c("mediation", "WeightIt", "ipw")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Install fixture backend: ", package, call. = FALSE)
  }
}

data(jobs, package = "mediation")
formula_m <- job_seek ~ treat + econ_hard + sex + age
formula_y <- depress2 ~ treat + job_seek + econ_hard + sex + age
model_m <- stats::lm(formula_m, data = jobs)
model_y <- stats::lm(formula_y, data = jobs)
set.seed(481)
mediation_reference <- mediation::mediate(
  model_m, model_y, treat = "treat", mediator = "job_seek",
  control.value = 0, treat.value = 1, sims = 500L,
  boot = FALSE, boot.ci.type = "perc", conf.level = 0.95,
  robustSE = FALSE)
medsens_reference <- mediation::medsens(
  mediation_reference, effect.type = "both")
mediation_fixture <- list(
  provenance = list(
    source = "mediation::jobs framing example",
    R = as.character(getRversion()),
    mediation = as.character(utils::packageVersion("mediation")),
    seed = 481L, sims = 500L
  ),
  data = jobs,
  formulas = list(
    mediator = paste(deparse(formula_m), collapse = " "),
    outcome = paste(deparse(formula_y), collapse = " ")
  ),
  settings = list(
    treat = "treat", mediator = "job_seek",
    control.value = 0, treat.value = 1,
    sims = 500L, boot = FALSE, conf.level = 0.95, seed = 481L
  ),
  effects = list(
    d0 = mediation_reference$d0,
    d1 = mediation_reference$d1,
    d_avg = mediation_reference$d.avg,
    z0 = mediation_reference$z0,
    z1 = mediation_reference$z1,
    z_avg = mediation_reference$z.avg,
    total = mediation_reference$tau.coef,
    d0_ci = mediation_reference$d0.ci,
    d1_ci = mediation_reference$d1.ci,
    d_avg_ci = mediation_reference$d.avg.ci,
    z0_ci = mediation_reference$z0.ci,
    z1_ci = mediation_reference$z1.ci,
    z_avg_ci = mediation_reference$z.avg.ci,
    total_ci = mediation_reference$tau.ci
  ),
  sensitivity = list(
    rho = medsens_reference$rho,
    d0 = medsens_reference$d0,
    d1 = medsens_reference$d1,
    z0 = medsens_reference$z0,
    z1 = medsens_reference$z1,
    err_cr_d = medsens_reference$err.cr.d,
    err_cr_z = medsens_reference$err.cr.z
  )
)
saveRDS(
  mediation_fixture, file.path(fixture_dir, "mediation-jobs.rds"),
  version = 2)

set.seed(1701)
n <- 200L
iptw_data <- data.frame(
  x = stats::rnorm(n),
  z = stats::rbinom(n, 1, 0.4)
)
iptw_data$treatment <- stats::rbinom(
  n, 1, stats::plogis(-0.1 + 0.5 * iptw_data$x - 0.3 * iptw_data$z))
formula <- treatment ~ x + z
fit <- stats::glm(formula, data = iptw_data, family = stats::binomial())
propensity <- as.numeric(stats::predict(fit, type = "response"))
internal_weights <- ifelse(
  iptw_data$treatment == 1, 1 / propensity, 1 / (1 - propensity))
weightit_fit <- WeightIt::weightit(
  formula, data = iptw_data, method = "glm",
  estimand = "ATE", stabilize = FALSE)
ipw_fit <- ipw::ipwpoint(
  exposure = treatment, family = "binomial", link = "logit",
  denominator = ~ x + z, data = iptw_data)
iptw_fixture <- list(
  provenance = list(
    source = "author-generated seeded logistic example",
    R = as.character(getRversion()),
    WeightIt = as.character(utils::packageVersion("WeightIt")),
    ipw = as.character(utils::packageVersion("ipw")),
    seed = 1701L, estimand = "ATE", stabilized = FALSE,
    formula = "treatment ~ x + z"
  ),
  data = iptw_data,
  propensity = propensity,
  weights = list(
    internal = internal_weights,
    WeightIt = as.numeric(weightit_fit$weights),
    ipw = as.numeric(ipw_fit$ipw.weights)
  )
)
saveRDS(
  iptw_fixture, file.path(fixture_dir, "baseline-iptw.rds"),
  version = 2)

input <- utils::read.csv(
  file.path(fixture_dir, "target-trial-input.csv"),
  stringsAsFactors = FALSE)
endpoint <- utils::read.csv(
  file.path(fixture_dir, "target-trial-expected-endpoint.csv"),
  stringsAsFactors = FALSE)
endpoint$strategy <- factor(endpoint$strategy, levels = c("never", "always"))
cox <- survival::coxph(
  survival::Surv(endpoint_time, endpoint_event) ~ strategy,
  data = endpoint, weights = analysis_weight,
  robust = TRUE, cluster = original_id, model = TRUE)
summary <- summary(cox, conf.int = 0.95)
target_fixture <- list(
  provenance = list(
    source = "author-reconstructed hand-auditable static-strategy example",
    R = as.character(getRversion()),
    survival = as.character(utils::packageVersion("survival")),
    denominator_probability = 0.5,
    numerator_probability = 0.4,
    interval_weight = 0.8
  ),
  input = input,
  expected_clones = utils::read.csv(
    file.path(fixture_dir, "target-trial-expected-clones.csv"),
    stringsAsFactors = FALSE),
  expected_endpoint = endpoint,
  survival = list(
    coefficient = unname(stats::coef(cox)),
    hazard_ratio = unname(exp(stats::coef(cox))),
    robust_se = unname(summary$coefficients[, "robust se"]),
    conf_low = unname(summary$conf.int[, "lower .95"]),
    conf_high = unname(summary$conf.int[, "upper .95"]),
    p_value = unname(summary$coefficients[, "Pr(>|z|)"])
  )
)
saveRDS(
  target_fixture, file.path(fixture_dir, "target-trial-reference.rds"),
  version = 2)
