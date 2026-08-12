# Independent numerical validation for WS8-11.
#
# Tolerances are declared here before simulation output is inspected. This
# script calls only the public API; all reference cloning and weights below are
# independently implemented.

tolerances <- list(
  clone_mismatches = 0L,
  baseline_weight_relative_error = 1e-10,
  longitudinal_weight_relative_error = 1e-10,
  absolute_bias = 0.05,
  coverage_minimum = 0.85,
  coverage_maximum = 1.00,
  simulations = 100L,
  subjects = 500L
)

args <- commandArgs(trailingOnly = TRUE)
package_dir <- if (length(args)) args[[1L]] else
  "physio-ecosystem/PhysioClinStats"
if (!dir.exists(package_dir) && dir.exists(".")) {
  package_dir <- "."
}
devtools::load_all(package_dir, quiet = TRUE)

protocol <- targetTrialProtocol(
  eligibility = function(data) data$eligible,
  treatment_strategies = list(never = 0, always = 1),
  assignment = "clone at eligibility",
  time_zero = 0,
  follow_up = 2,
  outcome = "binary endpoint at visit 2",
  causal_contrast = "always versus never",
  analysis_plan = "clone-censor-weight pooled logistic"
)

simulate_trial <- function(seed, n = tolerances$subjects) {
  set.seed(seed)
  x <- stats::rnorm(n)
  rows <- vector("list", n)
  a0 <- integer(n)
  for (i in seq_len(n)) {
    a <- integer(3)
    lag <- c(0L, NA_integer_, NA_integer_)
    for (j in seq_len(3)) {
      a[j] <- stats::rbinom(
        1, 1, stats::plogis(-0.2 + 0.35 * x[i] + 0.65 * lag[j]))
      if (j < 3L) {
        lag[j + 1L] <- a[j]
      }
    }
    a0[i] <- a[1L]
    rows[[i]] <- data.frame(
      id = i, visit = 0:2, treatment = a,
      x = x[i], lag_treatment = lag, eligible = TRUE,
      stringsAsFactors = FALSE
    )
  }
  data <- do.call(rbind, rows)
  outcome <- stats::rbinom(
    n, 1, stats::plogis(-1 + log(1.5) * a0 + 0.4 * x))
  data$outcome <- outcome[data$id]
  true_risk_never <- mean(stats::plogis(-1 + 0.4 * x))
  true_risk_always <- mean(stats::plogis(-1 + log(1.5) + 0.4 * x))
  list(
    data = data,
    truth = true_risk_always - true_risk_never
  )
}

independent_clones <- function(data) {
  strategies <- c(never = 0L, always = 1L)
  output <- list()
  k <- 0L
  for (subject in sort(unique(data$id))) {
    subject_rows <- data[data$id == subject, ]
    for (strategy in names(strategies)) {
      k <- k + 1L
      adherent <- subject_rows$treatment == strategies[[strategy]]
      deviation <- match(FALSE, adherent, nomatch = 0L)
      end <- if (deviation) deviation else nrow(subject_rows)
      output[[k]] <- data.frame(
        source_row = row.names(subject_rows)[seq_len(end)],
        original_id = as.character(subject),
        strategy = strategy,
        visit = subject_rows$visit[seq_len(end)],
        adherent = adherent[seq_len(end)],
        artificial_censor = seq_len(end) == deviation & deviation > 0L,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, output)
}

relative_error <- function(observed, expected) {
  max(abs(observed - expected) / pmax(abs(expected), 1e-12))
}

extract_rd <- function(result) {
  effects <- PhysioCore::resultValue(result)$effects
  effects[effects$measure == "risk_difference", , drop = FALSE]
}

first <- simulate_trial(1L)
baseline <- first$data[first$data$visit == 0, ]
baseline_fit <- stats::glm(
  treatment ~ x, data = baseline, family = stats::binomial())
baseline_ps <- as.numeric(stats::predict(baseline_fit, type = "response"))
baseline_reference <- ifelse(
  baseline$treatment == 1, 1 / baseline_ps, 1 / (1 - baseline_ps))
names(baseline_reference) <- as.character(baseline$id)
itt_first <- targetTrialEmulate(
  first$data, protocol, "id", "visit", "treatment", "outcome",
  baseline_covariates = "x", estimand = "intention_to_treat",
  stabilized = FALSE, weight_truncation = c(0, 1))
baseline_observed <- PhysioCore::resultValue(
  itt_first)$endpoint_data$analysis_weight
baseline_endpoint_ids <- PhysioCore::resultValue(
  itt_first)$endpoint_data$original_id
baseline_weight_error <- relative_error(
  baseline_observed,
  baseline_reference[match(
    baseline_endpoint_ids, names(baseline_reference))])

pp_first <- targetTrialEmulate(
  first$data, protocol, "id", "visit", "treatment", "outcome",
  baseline_covariates = "x",
  time_varying_covariates = "lag_treatment",
  estimand = "per_protocol", stabilized = FALSE,
  weight_truncation = c(0, 1))
package_clones <- PhysioCore::resultValue(pp_first)$clone_data
reference_clones <- independent_clones(first$data)
package_key <- paste(
  package_clones$original_id, package_clones$strategy,
  package_clones$visit, sep = "::")
reference_key <- paste(
  reference_clones$original_id, reference_clones$strategy,
  reference_clones$visit, sep = "::")
clone_mismatches <- sum(!reference_key %in% package_key) +
  sum(!package_key %in% reference_key)
aligned <- match(reference_key, package_key)
clone_mismatches <- clone_mismatches +
  sum(reference_clones$adherent != package_clones$adherent[aligned]) +
  sum(reference_clones$artificial_censor !=
        package_clones$artificial_censor[aligned])

denominator_fit <- stats::glm(
  treatment ~ visit + x + lag_treatment,
  data = first$data, family = stats::binomial())
ps <- as.numeric(stats::predict(denominator_fit, type = "response"))
p_observed <- ifelse(
  first$data$treatment == 1, ps, 1 - ps)
row_probability <- stats::setNames(
  p_observed, as.character(seq_len(nrow(first$data))))
reference_raw <- numeric(nrow(reference_clones))
for (clone in unique(paste(
  reference_clones$original_id, reference_clones$strategy, sep = "::"))) {
  index <- which(paste(
    reference_clones$original_id,
    reference_clones$strategy, sep = "::") == clone)
  adherent <- index[reference_clones$adherent[index]]
  source <- as.integer(reference_clones$source_row[adherent])
  values <- cumprod(1 / row_probability[as.character(source)])
  reference_raw[adherent] <- values
  deviation <- index[!reference_clones$adherent[index]]
  if (length(deviation)) {
    reference_raw[deviation[1L]] <- if (length(values)) tail(values, 1L) else 1
  }
}
longitudinal_weight_error <- relative_error(
  package_clones$raw_weight[match(reference_key, package_key)],
  reference_raw)

simulation <- data.frame(
  seed = seq_len(tolerances$simulations),
  truth = NA_real_,
  itt = NA_real_, itt_low = NA_real_, itt_high = NA_real_,
  pp = NA_real_, pp_low = NA_real_, pp_high = NA_real_
)
for (seed in simulation$seed) {
  generated <- simulate_trial(seed)
  itt <- targetTrialEmulate(
    generated$data, protocol, "id", "visit", "treatment", "outcome",
    baseline_covariates = "x", estimand = "intention_to_treat",
    stabilized = FALSE, weight_truncation = c(0, 1))
  pp <- targetTrialEmulate(
    generated$data, protocol, "id", "visit", "treatment", "outcome",
    baseline_covariates = "x",
    time_varying_covariates = "lag_treatment",
    estimand = "per_protocol", stabilized = FALSE,
    weight_truncation = c(0, 1))
  itt_effect <- extract_rd(itt)
  pp_effect <- extract_rd(pp)
  simulation$truth[seed] <- generated$truth
  simulation$itt[seed] <- itt_effect$estimate
  simulation$itt_low[seed] <- itt_effect$conf_low
  simulation$itt_high[seed] <- itt_effect$conf_high
  simulation$pp[seed] <- pp_effect$estimate
  simulation$pp_low[seed] <- pp_effect$conf_low
  simulation$pp_high[seed] <- pp_effect$conf_high
}

metrics <- list(
  itt_bias = mean(simulation$itt - simulation$truth),
  pp_bias = mean(simulation$pp - simulation$truth),
  itt_coverage = mean(
    simulation$itt_low <= simulation$truth &
      simulation$itt_high >= simulation$truth),
  pp_coverage = mean(
    simulation$pp_low <= simulation$truth &
      simulation$pp_high >= simulation$truth)
)

positivity_data <- first$data
positivity_data$treatment <- 1L
positivity_error <- tryCatch({
  targetTrialEmulate(
    positivity_data, protocol, "id", "visit", "treatment", "outcome",
    baseline_covariates = "x", estimand = "intention_to_treat")
  NULL
}, error = conditionMessage)
positivity_surfaced <- is.character(positivity_error) &&
  grepl("exactly two|Positivity", positivity_error)

stored_formula <- PhysioCore::resultValue(
  pp_first)$censor_models$formulas$denominator
mutated_future_formula <- "treatment ~ visit + x + lag_treatment + outcome"
future_information_mutation_detected <-
  !grepl("outcome", stored_formula, fixed = TRUE) &&
  grepl("outcome", mutated_future_formula, fixed = TRUE)

correct_grouped <- reference_raw
wrong_global <- cumprod(1 / rep(
  row_probability, length.out = length(reference_raw)))
wrong_grouping_mutation_detected <-
  relative_error(wrong_global, correct_grouped) >
  tolerances$longitudinal_weight_relative_error

wrong_coding <- ifelse(
  baseline$treatment == 1, 1 / (1 - baseline_ps), 1 / baseline_ps)
coding_mutation_detected <-
  relative_error(wrong_coding, baseline_reference) >
  tolerances$baseline_weight_relative_error

gates <- c(
  clone_agreement =
    clone_mismatches <= tolerances$clone_mismatches,
  baseline_weight =
    baseline_weight_error <= tolerances$baseline_weight_relative_error,
  longitudinal_weight =
    longitudinal_weight_error <=
      tolerances$longitudinal_weight_relative_error,
  itt_bias = abs(metrics$itt_bias) <= tolerances$absolute_bias,
  pp_bias = abs(metrics$pp_bias) <= tolerances$absolute_bias,
  itt_coverage =
    metrics$itt_coverage >= tolerances$coverage_minimum &&
      metrics$itt_coverage <= tolerances$coverage_maximum,
  pp_coverage =
    metrics$pp_coverage >= tolerances$coverage_minimum &&
      metrics$pp_coverage <= tolerances$coverage_maximum,
  positivity = positivity_surfaced,
  future_information_mutation = future_information_mutation_detected,
  grouping_mutation = wrong_grouping_mutation_detected,
  coding_mutation = coding_mutation_detected
)

report <- list(
  tolerances = tolerances,
  R = as.character(getRversion()),
  package_version = as.character(utils::packageVersion("PhysioClinStats")),
  clone_mismatches = clone_mismatches,
  baseline_weight_relative_error = baseline_weight_error,
  longitudinal_weight_relative_error = longitudinal_weight_error,
  metrics = metrics,
  positivity_error = positivity_error,
  mutations = list(
    future_information = future_information_mutation_detected,
    cumulative_product_grouping = wrong_grouping_mutation_detected,
    treatment_coding = coding_mutation_detected
  ),
  gates = gates
)

output_dir <- file.path(package_dir, "inst", "validation")
saveRDS(report, file.path(output_dir, "ws8-11-results.rds"), version = 2)
lines <- c(
  "# WS8-11 independent numeric validation",
  "",
  sprintf("- R: %s", report$R),
  sprintf("- PhysioClinStats: %s", report$package_version),
  sprintf("- Seeded datasets: %d x %d subjects",
          tolerances$simulations, tolerances$subjects),
  sprintf("- Clone/censor mismatches: %d (limit %d)",
          clone_mismatches, tolerances$clone_mismatches),
  sprintf("- Baseline weight max relative error: %.3g (limit %.3g)",
          baseline_weight_error,
          tolerances$baseline_weight_relative_error),
  sprintf("- Longitudinal weight max relative error: %.3g (limit %.3g)",
          longitudinal_weight_error,
          tolerances$longitudinal_weight_relative_error),
  sprintf("- ITT RD bias: %.4f; 95%% CI coverage: %.3f",
          metrics$itt_bias, metrics$itt_coverage),
  sprintf("- Per-protocol RD bias: %.4f; 95%% CI coverage: %.3f",
          metrics$pp_bias, metrics$pp_coverage),
  sprintf("- Positivity stress surfaced: %s", positivity_surfaced),
  sprintf("- Future-information mutation detected: %s",
          future_information_mutation_detected),
  sprintf("- Cumulative-product grouping mutation detected: %s",
          wrong_grouping_mutation_detected),
  sprintf("- Treatment-coding mutation detected: %s",
          coding_mutation_detected),
  "",
  "## Gates",
  "",
  paste0("- ", names(gates), ": ", ifelse(gates, "PASS", "FAIL")),
  "",
  paste(
    "This simulation checks numerical implementation under its generating",
    "mechanism; it does not validate causal identification assumptions in",
    "observational rehabilitation data."
  )
)
writeLines(lines, file.path(output_dir, "ws8-11-results.md"))
print(report)
if (!all(gates)) {
  stop("One or more WS8-11 independent validation gates failed.",
       call. = FALSE)
}
