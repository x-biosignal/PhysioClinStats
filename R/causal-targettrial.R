# Longitudinal target-trial emulation.

.tt_match_arg <- function(x, choices, name) {
  if (length(x) > 1L) {
    x <- x[1L]
  }
  if (!is.character(x) || length(x) != 1L || is.na(x) || !x %in% choices) {
    stop(sprintf("'%s' must be one of: %s.", name,
                 paste(choices, collapse = ", ")), call. = FALSE)
  }
  x
}

.tt_column <- function(data, x, argument) {
  x <- .causal_scalar_string(x, argument)
  if (!x %in% names(data)) {
    stop(sprintf("Column '%s' supplied as '%s' is absent from 'data'.",
                 x, argument), call. = FALSE)
  }
  x
}

.tt_equal <- function(x, y) {
  as.character(x) == as.character(y)
}

.tt_binary_code <- function(x, name) {
  if (anyNA(x)) {
    stop(sprintf("'%s' contains missing values.", name), call. = FALSE)
  }
  levels <- sort(unique(as.character(x)))
  if (length(levels) != 2L) {
    stop(sprintf("'%s' must have exactly two observed values.", name),
         call. = FALSE)
  }
  list(value = as.integer(as.character(x) == levels[2L]), levels = levels)
}

.tt_effective_n <- function(w) {
  if (!length(w) || !all(is.finite(w)) || sum(w ^ 2) == 0) {
    return(NA_real_)
  }
  sum(w) ^ 2 / sum(w ^ 2)
}

.tt_weighted_mean <- function(x, w) {
  sum(x * w) / sum(w)
}

.tt_weighted_var <- function(x, w) {
  m <- .tt_weighted_mean(x, w)
  sum(w * (x - m) ^ 2) / sum(w)
}

.tt_balance <- function(treatment, covariates, weights = NULL) {
  if (!ncol(covariates)) {
    return(data.frame(
      covariate = character(), level = character(),
      smd = numeric(), stringsAsFactors = FALSE
    ))
  }
  a <- as.integer(treatment)
  w <- if (is.null(weights)) rep(1, length(a)) else weights
  rows <- list()
  add_row <- function(x, covariate, level = "") {
    m1 <- .tt_weighted_mean(x[a == 1L], w[a == 1L])
    m0 <- .tt_weighted_mean(x[a == 0L], w[a == 0L])
    v1 <- .tt_weighted_var(x[a == 1L], w[a == 1L])
    v0 <- .tt_weighted_var(x[a == 0L], w[a == 0L])
    denom <- sqrt((v1 + v0) / 2)
    smd <- if (is.finite(denom) && denom > 0) (m1 - m0) / denom else
      if (isTRUE(all.equal(m1, m0))) 0 else sign(m1 - m0) * Inf
    data.frame(
      covariate = covariate, level = level, smd = smd,
      stringsAsFactors = FALSE
    )
  }
  for (nm in names(covariates)) {
    x <- covariates[[nm]]
    if (is.numeric(x) || is.integer(x)) {
      rows[[length(rows) + 1L]] <- add_row(as.numeric(x), nm)
    } else {
      for (lev in sort(unique(as.character(x)))) {
        rows[[length(rows) + 1L]] <- add_row(
          as.integer(as.character(x) == lev), nm, lev)
      }
    }
  }
  do.call(rbind, rows)
}

.tt_weight_diagnostics <- function(a, ps, weights, covariates) {
  before <- .tt_balance(a, covariates)
  after <- .tt_balance(a, covariates, weights)
  if (nrow(before)) {
    balance <- merge(
      before, after, by = c("covariate", "level"),
      suffixes = c("_before", "_after"), sort = FALSE
    )
  } else {
    balance <- data.frame(
      covariate = character(), level = character(),
      smd_before = numeric(), smd_after = numeric(),
      stringsAsFactors = FALSE
    )
  }
  list(
    propensity = c(min = min(ps), max = max(ps)),
    weight_quantiles = stats::quantile(
      weights, c(0, 0.01, 0.5, 0.99, 1), names = TRUE),
    max_weight = max(weights),
    effective_sample_size = c(
      control = .tt_effective_n(weights[a == 0L]),
      treated = .tt_effective_n(weights[a == 1L])
    ),
    balance = balance
  )
}

.tt_internal_ps <- function(a, covariates) {
  d <- data.frame(.treatment = a, covariates, check.names = FALSE)
  rhs <- names(d)[-1L]
  formula <- if (length(rhs)) {
    stats::reformulate(rhs, response = ".treatment")
  } else {
    .treatment ~ 1
  }
  fit <- suppressWarnings(stats::glm(
    formula, data = d, family = stats::binomial(),
    na.action = stats::na.fail, model = TRUE, x = TRUE, y = TRUE
  ))
  list(
    propensity = as.numeric(stats::predict(fit, type = "response")),
    fit = fit,
    formula = formula
  )
}

# Private parity helper used by target-trial weighting.
.baseline_iptw <- function(
    treatment,
    covariates,
    estimand = c("ATE", "ATT"),
    stabilized = FALSE,
    backend = c("internal", "WeightIt", "ipw")) {
  estimand <- .tt_match_arg(estimand, c("ATE", "ATT"), "estimand")
  backend <- .tt_match_arg(
    backend, c("internal", "WeightIt", "ipw"), "backend")
  if (!is.data.frame(covariates)) {
    covariates <- as.data.frame(covariates, stringsAsFactors = FALSE)
  }
  if (nrow(covariates) != length(treatment)) {
    stop("'covariates' must have one row per treatment value.", call. = FALSE)
  }
  if (any(!stats::complete.cases(covariates))) {
    stop("IPTW covariates must be complete.", call. = FALSE)
  }
  code <- .tt_binary_code(treatment, "treatment")
  a <- code$value
  marginal <- mean(a)
  backend_fit <- NULL

  if (backend == "internal") {
    fitted <- .tt_internal_ps(a, covariates)
    ps <- fitted$propensity
    backend_fit <- fitted$fit
    formula <- fitted$formula
  } else if (backend == "WeightIt") {
    requireBackend("WeightIt")
    d <- data.frame(.treatment = a, covariates, check.names = FALSE)
    formula <- if (ncol(covariates)) {
      stats::reformulate(names(covariates), response = ".treatment")
    } else {
      .treatment ~ 1
    }
    backend_fit <- WeightIt::weightit(
      formula, data = d, method = "glm", estimand = estimand,
      stabilize = stabilized
    )
    weights <- as.numeric(backend_fit$weights)
    ps <- backend_fit$ps
    if (is.null(ps)) {
      ps <- ifelse(a == 1L, 1 / weights, 1 - 1 / weights)
    }
  } else {
    requireBackend("ipw")
    if (estimand != "ATE") {
      stop("The ipw backend parity path supports the ATE estimand only.",
           call. = FALSE)
    }
    d <- data.frame(.treatment = a, covariates, check.names = FALSE)
    denominator <- if (ncol(covariates)) {
      stats::reformulate(names(covariates))
    } else {
      ~ 1
    }
    ipw_args <- list(
      exposure = quote(.treatment),
      family = "binomial",
      link = "logit",
      denominator = denominator,
      data = d
    )
    if (stabilized) {
      ipw_args$numerator <- ~ 1
    }
    backend_fit <- do.call(ipw::ipwpoint, ipw_args)
    weights <- as.numeric(backend_fit$ipw.weights)
    ps <- as.numeric(stats::predict(
      backend_fit$den.mod, type = "response"))
    formula <- denominator
  }

  eps <- 1e-6
  if (any(!is.finite(ps)) || any(ps <= eps | ps >= 1 - eps)) {
    stop(sprintf(
      "Positivity failure: propensity scores must lie in (%g, %g).",
      eps, 1 - eps), call. = FALSE)
  }
  if (backend == "internal") {
    if (estimand == "ATE") {
      weights <- ifelse(a == 1L, 1 / ps, 1 / (1 - ps))
      if (stabilized) {
        weights <- weights * ifelse(a == 1L, marginal, 1 - marginal)
      }
    } else {
      weights <- ifelse(a == 1L, 1, ps / (1 - ps))
      if (stabilized) {
        weights[a == 0L] <- weights[a == 0L] *
          ((1 - marginal) / marginal)
      }
    }
  }
  if (any(!is.finite(weights)) || any(weights <= 0)) {
    stop("The IPTW backend returned non-finite or non-positive weights.",
         call. = FALSE)
  }
  list(
    weights = weights,
    propensity = as.numeric(ps),
    treatment = a,
    treatment_levels = code$levels,
    estimand = estimand,
    stabilized = stabilized,
    backend = backend,
    formula = formula,
    fit = backend_fit,
    diagnostics = .tt_weight_diagnostics(a, ps, weights, covariates)
  )
}

.tt_eval_rule <- function(rule, data, name) {
  value <- if (is.function(rule)) {
    rule(data)
  } else if (inherits(rule, "formula") && length(rule) == 2L) {
    eval(rule[[2L]], data, environment(rule))
  } else if (is.logical(rule)) {
    rule
  } else {
    stop(sprintf(
      "Protocol '%s' must be an executable function, one-sided formula, or logical rule.",
      name), call. = FALSE)
  }
  if (!is.logical(value) || length(value) != nrow(data) || anyNA(value)) {
    stop(sprintf("Protocol '%s' must return one non-missing logical per row.",
                 name), call. = FALSE)
  }
  value
}

.tt_time_zero_rows <- function(rule, data, time) {
  if (is.function(rule) ||
      (inherits(rule, "formula") && length(rule) == 2L) ||
      is.logical(rule)) {
    .tt_eval_rule(rule, data, "time_zero")
  } else if (length(rule) == 1L && is.atomic(rule) && !is.na(rule)) {
    .tt_equal(data[[time]], rule)
  } else {
    stop("Protocol 'time_zero' must be a scalar value or executable rule.",
         call. = FALSE)
  }
}

.tt_prepare_data <- function(
    data, protocol, id, time, treatment, outcome, event,
    baseline_covariates, time_varying_covariates) {
  if (!is.data.frame(data) || !nrow(data)) {
    stop("'data' must be a non-empty data.frame.", call. = FALSE)
  }
  if (!inherits(protocol, "target_trial_protocol")) {
    stop("'protocol' must be a target_trial_protocol.", call. = FALSE)
  }
  id <- .tt_column(data, id, "id")
  time <- .tt_column(data, time, "time")
  treatment <- .tt_column(data, treatment, "treatment")
  outcome <- .tt_column(data, outcome, "outcome")
  if (!is.null(event)) {
    event <- .tt_column(data, event, "event")
  }
  baseline_covariates <- unique(as.character(baseline_covariates))
  time_varying_covariates <- unique(as.character(time_varying_covariates))
  model_columns <- unique(c(
    id, time, treatment, outcome, event,
    baseline_covariates, time_varying_covariates
  ))
  missing_columns <- setdiff(model_columns, names(data))
  if (length(missing_columns)) {
    stop(sprintf("Model column(s) absent from data: %s.",
                 paste(missing_columns, collapse = ", ")), call. = FALSE)
  }
  if (any(!stats::complete.cases(data[, model_columns, drop = FALSE]))) {
    stop("Columns used by treatment, censoring, or outcome models must be complete.",
         call. = FALSE)
  }
  if (anyDuplicated(data[, c(id, time), drop = FALSE])) {
    stop("Each id,time pair must identify exactly one row.", call. = FALSE)
  }
  if (!is.numeric(data[[time]]) || any(!is.finite(data[[time]]))) {
    stop("The time column must be finite numeric.", call. = FALSE)
  }
  row_groups <- split(seq_len(nrow(data)), as.character(data[[id]]))
  increasing <- vapply(row_groups, function(i) {
    length(i) == 1L || all(diff(data[[time]][i]) > 0)
  }, logical(1))
  if (!all(increasing)) {
    stop("Rows must have strictly increasing time within each subject; input is not silently sorted.",
         call. = FALSE)
  }

  time_zero <- .tt_time_zero_rows(protocol$time_zero, data, time)
  zero_counts <- vapply(row_groups, function(i) sum(time_zero[i]), integer(1))
  if (any(zero_counts != 1L)) {
    stop("Every subject must have exactly one declared time-zero row.",
         call. = FALSE)
  }
  zero_index <- vapply(row_groups, function(i) i[time_zero[i]], integer(1))
  zero_times <- data[[time]][zero_index]
  if (length(unique(zero_times)) != 1L) {
    stop("All subjects must share a common eligible time zero.", call. = FALSE)
  }
  before_zero <- vapply(row_groups, function(i) {
    any(data[[time]][i] < data[[time]][i[time_zero[i]]])
  }, logical(1))
  if (any(before_zero)) {
    stop("Rows before the declared time zero are not allowed.", call. = FALSE)
  }
  baseline <- data[zero_index, , drop = FALSE]
  eligible <- .tt_eval_rule(protocol$eligibility, baseline, "eligibility")
  eligible_ids <- as.character(baseline[[id]][eligible])
  if (!length(eligible_ids)) {
    stop("No participant satisfies the declared eligibility rule.",
         call. = FALSE)
  }
  keep <- as.character(data[[id]]) %in% eligible_ids
  d <- data[keep, , drop = FALSE]
  d$.source_row <- which(keep)

  if (!is.numeric(d[[outcome]]) || any(!is.finite(d[[outcome]]))) {
    stop("The outcome column must be finite numeric.", call. = FALSE)
  }
  eligible_groups <- split(seq_len(nrow(d)), as.character(d[[id]]))
  constant_outcome <- vapply(eligible_groups, function(i) {
    length(unique(d[[outcome]][i])) == 1L
  }, logical(1))
  if (!all(constant_outcome)) {
    stop("Outcome time/value must be constant within subject in long data.",
         call. = FALSE)
  }
  if (!is.null(event)) {
    event_values <- unique(d[[event]])
    if (!all(event_values %in% c(0, 1, FALSE, TRUE))) {
      stop("The event column must be binary (0/1).", call. = FALSE)
    }
    constant_event <- vapply(eligible_groups, function(i) {
      length(unique(d[[event]][i])) == 1L
    }, logical(1))
    if (!all(constant_event)) {
      stop("Event status must be constant within subject in long data.",
           call. = FALSE)
    }
    if (any(d[[time]] > d[[outcome]])) {
      stop("Post-outcome rows are not allowed.", call. = FALSE)
    }
  } else {
    if (!all(unique(d[[outcome]]) %in% c(0, 1))) {
      stop("Without 'event', the outcome must be binary (0/1).",
           call. = FALSE)
    }
    if (!is.numeric(protocol$follow_up) ||
        length(protocol$follow_up) != 1L ||
        !is.finite(protocol$follow_up) || protocol$follow_up <= 0) {
      stop("Binary target-trial emulation requires a finite positive scalar follow-up horizon.",
           call. = FALSE)
    }
    if (any(d[[time]] > protocol$follow_up)) {
      stop("Rows after the protocol follow-up horizon are not allowed.",
           call. = FALSE)
    }
  }
  code <- .tt_binary_code(d[[treatment]], "treatment")
  static <- protocol$treatment_strategies[
    vapply(protocol$strategy_specs, `[[`, character(1), "type") == "static"]
  impossible <- names(static)[!vapply(static, function(x) {
    as.character(x) %in% code$levels
  }, logical(1))]
  if (length(impossible)) {
    stop(sprintf("Strategy value(s) are absent from observed treatment: %s.",
                 paste(impossible, collapse = ", ")), call. = FALSE)
  }
  list(
    data = d, baseline = baseline[eligible, , drop = FALSE],
    eligible_ids = eligible_ids,
    counts = list(
      input_subjects = length(row_groups),
      eligible_subjects = length(eligible_ids),
      excluded_subjects = length(row_groups) - length(eligible_ids)
    ),
    columns = list(
      id = id, time = time, treatment = treatment,
      outcome = outcome, event = event,
      baseline = baseline_covariates,
      time_varying = time_varying_covariates
    ),
    treatment_levels = code$levels,
    time_zero = unique(zero_times)
  )
}

.tt_expected_treatment <- function(strategy, rows, columns) {
  if (!is.function(strategy)) {
    return(rep(strategy, nrow(rows)))
  }
  available <- setdiff(
    names(rows), c(columns$outcome, columns$event))
  expected <- vector("list", nrow(rows))
  for (j in seq_len(nrow(rows))) {
    history <- rows[seq_len(j), available, drop = FALSE]
    history[[columns$treatment]][j] <- NA
    value <- strategy(
      data = history[j, , drop = FALSE],
      history = history
    )
    if (length(value) != 1L || is.na(value) || !is.atomic(value)) {
      stop("A dynamic strategy must return one non-missing treatment value per row.",
           call. = FALSE)
    }
    expected[[j]] <- value
  }
  unlist(expected, recursive = FALSE, use.names = FALSE)
}

.tt_fit_interval_probabilities <- function(
    data, columns, backend, stabilized, numerator_covariates) {
  code <- .tt_binary_code(data[[columns$treatment]], "treatment")
  denominator_names <- unique(c(
    columns$time, columns$baseline, columns$time_varying))
  denominator_covariates <- data[, denominator_names, drop = FALSE]
  denominator <- .baseline_iptw(
    code$value, denominator_covariates, estimand = "ATE",
    stabilized = FALSE, backend = backend
  )
  ps_denom <- denominator$propensity
  p_denom <- ifelse(code$value == 1L, ps_denom, 1 - ps_denom)

  numerator_covariates <- unique(as.character(numerator_covariates))
  missing_numerator <- setdiff(numerator_covariates, names(data))
  if (length(missing_numerator)) {
    stop(sprintf("Numerator covariate(s) absent from data: %s.",
                 paste(missing_numerator, collapse = ", ")), call. = FALSE)
  }
  if (stabilized) {
    numerator_names <- unique(c(columns$time, numerator_covariates))
    if (any(!stats::complete.cases(
      data[, numerator_names, drop = FALSE]))) {
      stop("Numerator-model covariates must be complete.", call. = FALSE)
    }
    numerator <- .tt_internal_ps(
      code$value, data[, numerator_names, drop = FALSE])
    ps_num <- numerator$propensity
    p_num <- ifelse(code$value == 1L, ps_num, 1 - ps_num)
  } else {
    numerator <- NULL
    p_num <- rep(1, nrow(data))
  }
  eps <- 1e-6
  if (any(!is.finite(p_denom)) || any(p_denom <= eps) ||
      any(!is.finite(p_num)) || any(p_num <= 0)) {
    stop("Positivity failure in interval adherence probabilities.",
         call. = FALSE)
  }
  list(
    denominator_probability = p_denom,
    numerator_probability = p_num,
    denominator = denominator,
    numerator = numerator,
    formulas = list(
      denominator = paste(deparse(denominator$formula), collapse = " "),
      numerator = if (is.null(numerator)) "1" else
        paste(deparse(numerator$formula), collapse = " ")
    )
  )
}

.tt_clone <- function(data, protocol, columns, probabilities, horizon) {
  groups <- split(seq_len(nrow(data)), as.character(data[[columns$id]]))
  clones <- vector("list", length(groups) * length(protocol$treatment_strategies))
  k <- 0L
  for (subject in names(groups)) {
    rows <- data[groups[[subject]], , drop = FALSE]
    for (strategy_name in names(protocol$treatment_strategies)) {
      k <- k + 1L
      strategy <- protocol$treatment_strategies[[strategy_name]]
      expected <- .tt_expected_treatment(strategy, rows, columns)
      if (any(!as.character(expected) %in%
              sort(unique(as.character(data[[columns$treatment]]))))) {
        stop(sprintf(
          "Strategy '%s' returned a treatment value absent from observed data.",
          strategy_name), call. = FALSE)
      }
      adherent <- .tt_equal(rows[[columns$treatment]], expected)
      first_deviation <- match(FALSE, adherent, nomatch = 0L)
      end <- if (first_deviation) first_deviation else nrow(rows)
      clone <- rows[seq_len(end), , drop = FALSE]
      clone$original_id <- subject
      clone$strategy <- strategy_name
      clone$clone_id <- paste(subject, strategy_name, sep = "::")
      clone$expected_treatment <- expected[seq_len(end)]
      clone$adherent <- adherent[seq_len(end)]
      clone$artificial_censor <- FALSE
      if (first_deviation) {
        clone$artificial_censor[end] <- TRUE
      }
      clone$censor_reason <- NA_character_
      event_status <- if (is.null(columns$event)) 0L else
        as.integer(rows[[columns$event]][1L])
      outcome_time <- if (is.null(columns$event)) horizon else
        rows[[columns$outcome]][1L]
      deviation_time <- if (first_deviation) {
        clone[[columns$time]][end]
      } else {
        Inf
      }
      reason <- if (first_deviation && deviation_time <= outcome_time) {
        "artificial_deviation"
      } else if (!is.null(columns$event) && event_status == 1L &&
                 outcome_time <= horizon) {
        "outcome"
      } else if (max(rows[[columns$time]]) >=
                 min(horizon, outcome_time)) {
        "administrative"
      } else {
        "loss_to_follow_up"
      }
      clone$censor_reason[end] <- reason
      clone$denominator_probability <- NA_real_
      clone$numerator_probability <- NA_real_
      clone$interval_weight <- NA_real_
      clone$raw_weight <- NA_real_
      analysis_index <- which(clone$adherent)
      if (length(analysis_index)) {
        source <- match(clone$.source_row[analysis_index], data$.source_row)
        pd <- probabilities$denominator_probability[source]
        pn <- probabilities$numerator_probability[source]
        interval <- pn / pd
        raw <- cumprod(interval)
        clone$denominator_probability[analysis_index] <- pd
        clone$numerator_probability[analysis_index] <- pn
        clone$interval_weight[analysis_index] <- interval
        clone$raw_weight[analysis_index] <- raw
      }
      if (first_deviation) {
        prior <- if (end > 1L) clone$raw_weight[end - 1L] else 1
        clone$raw_weight[end] <- prior
      }
      clones[[k]] <- clone
    }
  }
  do.call(rbind, clones)
}

.tt_apply_truncation <- function(clones, probs, max_weight) {
  if (!is.numeric(probs) || length(probs) != 2L ||
      any(!is.finite(probs)) || probs[1L] < 0 || probs[2L] > 1 ||
      probs[1L] >= probs[2L]) {
    stop("'weight_truncation' must be increasing probabilities in [0, 1].",
         call. = FALSE)
  }
  if (!is.numeric(max_weight) || length(max_weight) != 1L ||
      is.na(max_weight) || max_weight <= 0) {
    stop("'max_weight' must be one positive number or Inf.", call. = FALSE)
  }
  valid <- clones$adherent & is.finite(clones$raw_weight)
  thresholds <- stats::quantile(
    clones$raw_weight[valid], probs = probs, names = FALSE, type = 7)
  analysis <- clones$raw_weight
  analysis[valid] <- pmax(thresholds[1L],
                          pmin(thresholds[2L], analysis[valid]))
  analysis[valid] <- pmin(analysis[valid], max_weight)
  clones$analysis_weight <- analysis
  clones$weight_truncated <- valid &
    abs(clones$analysis_weight - clones$raw_weight) >
    sqrt(.Machine$double.eps)
  list(
    data = clones,
    diagnostics = list(
      probabilities = probs,
      thresholds = stats::setNames(thresholds, c("lower", "upper")),
      max_weight = max_weight,
      truncation_count = sum(clones$weight_truncated),
      truncated_rows = data.frame(
        row = which(clones$weight_truncated),
        clone_id = clones$clone_id[clones$weight_truncated],
        source_row = clones$.source_row[clones$weight_truncated],
        stringsAsFactors = FALSE
      )
    )
  )
}

.tt_clone_endpoint <- function(clones, columns, horizon) {
  groups <- split(seq_len(nrow(clones)), clones$clone_id)
  rows <- lapply(groups, function(i) {
    x <- clones[i, , drop = FALSE]
    deviation <- which(x$artificial_censor)
    deviation_time <- if (length(deviation)) {
      x[[columns$time]][deviation[1L]]
    } else {
      Inf
    }
    outcome_time <- if (is.null(columns$event)) horizon else
      x[[columns$outcome]][1L]
    event_status <- if (is.null(columns$event)) NA_integer_ else
      as.integer(x[[columns$event]][1L])
    observed_until <- max(x[[columns$time]])
    lost <- identical(
      utils::tail(x$censor_reason, 1L), "loss_to_follow_up")
    endpoint_time <- min(outcome_time, horizon, deviation_time, observed_until)
    endpoint_event <- if (is.null(columns$event)) NA_integer_ else
      as.integer(event_status == 1L && outcome_time <= horizon &&
                   outcome_time < deviation_time &&
                   outcome_time <= observed_until)
    endpoint_outcome <- if (is.null(columns$event) && !lost &&
                            deviation_time > horizon &&
                            observed_until >= horizon) {
      as.integer(x[[columns$outcome]][1L])
    } else if (is.null(columns$event)) {
      NA_integer_
    } else {
      NA_integer_
    }
    adherent_rows <- which(x$adherent & is.finite(x$analysis_weight))
    endpoint_weight <- if (length(adherent_rows)) {
      x$analysis_weight[max(adherent_rows)]
    } else {
      1
    }
    data.frame(
      original_id = x$original_id[1L],
      clone_id = x$clone_id[1L],
      strategy = x$strategy[1L],
      endpoint_time = endpoint_time,
      endpoint_event = endpoint_event,
      endpoint_outcome = endpoint_outcome,
      analysis_weight = endpoint_weight,
      censor_reason = utils::tail(stats::na.omit(x$censor_reason), 1L),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.tt_itt_endpoint <- function(prepared, protocol, horizon, weights) {
  data <- prepared$data
  columns <- prepared$columns
  groups <- split(seq_len(nrow(data)), as.character(data[[columns$id]]))
  strategy_values <- protocol$treatment_strategies
  if (any(vapply(protocol$strategy_specs, `[[`, character(1), "type") !=
          "static")) {
    stop("Intention-to-treat emulation requires static treatment strategies.",
         call. = FALSE)
  }
  rows <- lapply(seq_along(groups), function(j) {
    x <- data[groups[[j]], , drop = FALSE]
    baseline_value <- x[[columns$treatment]][1L]
    matches <- names(strategy_values)[vapply(strategy_values, function(value) {
      .tt_equal(baseline_value, value)
    }, logical(1))]
    if (length(matches) != 1L) {
      stop("Each baseline treatment must map to exactly one static strategy.",
           call. = FALSE)
    }
    outcome_time <- if (is.null(columns$event)) horizon else
      x[[columns$outcome]][1L]
    observed_until <- max(x[[columns$time]])
    data.frame(
      original_id = as.character(x[[columns$id]][1L]),
      strategy = matches,
      endpoint_time = min(outcome_time, horizon, observed_until),
      endpoint_event = if (is.null(columns$event)) NA_integer_ else
        as.integer(x[[columns$event]][1L] == 1L &&
                     outcome_time <= horizon &&
                     outcome_time <= observed_until),
      endpoint_outcome = if (is.null(columns$event) &&
                             observed_until >= horizon) {
        as.integer(x[[columns$outcome]][1L])
      } else {
        NA_integer_
      },
      analysis_weight = weights[j],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.tt_cluster_vcov_glm <- function(fit, cluster) {
  x <- stats::model.matrix(fit)
  y <- fit$y
  mu <- stats::fitted(fit)
  w <- fit$prior.weights
  information <- crossprod(x, x * as.numeric(w * mu * (1 - mu)))
  bread <- tryCatch(solve(information), error = function(e) {
    stop("Outcome model information matrix is singular.", call. = FALSE)
  })
  score <- x * as.numeric(w * (y - mu))
  cluster_score <- rowsum(score, as.character(cluster), reorder = FALSE)
  meat <- crossprod(cluster_score)
  g <- nrow(cluster_score)
  n <- nrow(x)
  p <- ncol(x)
  correction <- if (g > 1L && n > p) {
    (g / (g - 1)) * ((n - 1) / (n - p))
  } else {
    1
  }
  correction * bread %*% meat %*% bread
}

.tt_survival_effect <- function(endpoint, strategy_levels, conf.level) {
  requireBackend("survival")
  if (!any(endpoint$endpoint_event == 1L)) {
    warning("The target-trial endpoint has 0 events; hazard ratios are undefined.",
            call. = FALSE)
    effects <- data.frame(
      contrast = paste(strategy_levels[-1L], "vs", strategy_levels[1L]),
      measure = "hazard_ratio", estimate = NA_real_,
      std_error = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
      p_value = NA_real_, stringsAsFactors = FALSE
    )
    return(list(effects = effects, fit = NULL,
                warning = "zero events"))
  }
  endpoint$strategy_factor <- factor(
    endpoint$strategy, levels = strategy_levels)
  fit <- survival::coxph(
    survival::Surv(endpoint_time, endpoint_event) ~ strategy_factor,
    data = endpoint, weights = endpoint$analysis_weight,
    robust = TRUE, cluster = endpoint$original_id, model = TRUE
  )
  summary <- summary(fit, conf.int = conf.level)
  coefficients <- summary$coefficients
  intervals <- summary$conf.int
  effects <- data.frame(
    contrast = paste(strategy_levels[-1L], "vs", strategy_levels[1L]),
    measure = "hazard_ratio",
    estimate = as.numeric(coefficients[, "exp(coef)"]),
    std_error = as.numeric(coefficients[, "robust se"]),
    conf_low = as.numeric(intervals[, 3L]),
    conf_high = as.numeric(intervals[, 4L]),
    p_value = as.numeric(coefficients[, "Pr(>|z|)"]),
    stringsAsFactors = FALSE
  )
  list(effects = effects, fit = fit, warning = NULL)
}

.tt_binary_effect <- function(endpoint, strategy_levels, conf.level) {
  endpoint <- endpoint[!is.na(endpoint$endpoint_outcome), , drop = FALSE]
  if (!nrow(endpoint) ||
      length(unique(endpoint$endpoint_outcome)) < 2L) {
    warning("The fixed-horizon outcome has fewer than two observed values; contrasts are undefined.",
            call. = FALSE)
    effects <- expand.grid(
      contrast = paste(strategy_levels[-1L], "vs", strategy_levels[1L]),
      measure = c("risk_difference", "risk_ratio"),
      stringsAsFactors = FALSE
    )
    effects$estimate <- effects$std_error <- effects$conf_low <-
      effects$conf_high <- effects$p_value <- NA_real_
    return(list(effects = effects, fit = NULL, risks = NULL,
                robust_vcov = NULL,
                warning = "degenerate fixed-horizon outcome"))
  }
  endpoint$strategy_factor <- factor(
    endpoint$strategy, levels = strategy_levels)
  fit <- suppressWarnings(stats::glm(
    endpoint_outcome ~ strategy_factor,
    data = endpoint, weights = endpoint$analysis_weight,
    family = stats::binomial(), x = TRUE, y = TRUE, model = TRUE
  ))
  vcov <- .tt_cluster_vcov_glm(fit, endpoint$original_id)
  zcrit <- stats::qnorm(1 - (1 - conf.level) / 2)
  design <- stats::model.matrix(
    ~ strategy_factor,
    data = data.frame(
      strategy_factor = factor(strategy_levels, levels = strategy_levels)))
  eta <- as.numeric(design %*% stats::coef(fit))
  risks <- stats::plogis(eta)
  risk_gradient <- design * as.numeric(risks * (1 - risks))
  rows <- list()
  for (j in 2:length(strategy_levels)) {
    contrast <- paste(strategy_levels[j], "vs", strategy_levels[1L])
    gradient_rd <- risk_gradient[j, ] - risk_gradient[1L, ]
    rd <- risks[j] - risks[1L]
    se_rd <- sqrt(drop(gradient_rd %*% vcov %*% gradient_rd))
    gradient_log_rr <- risk_gradient[j, ] / risks[j] -
      risk_gradient[1L, ] / risks[1L]
    rr <- risks[j] / risks[1L]
    se_log_rr <- sqrt(drop(
      gradient_log_rr %*% vcov %*% gradient_log_rr))
    rows[[length(rows) + 1L]] <- data.frame(
      contrast = contrast, measure = "risk_difference",
      estimate = rd, std_error = se_rd,
      conf_low = rd - zcrit * se_rd,
      conf_high = rd + zcrit * se_rd,
      p_value = 2 * stats::pnorm(-abs(rd / se_rd)),
      stringsAsFactors = FALSE
    )
    rows[[length(rows) + 1L]] <- data.frame(
      contrast = contrast, measure = "risk_ratio",
      estimate = rr, std_error = se_log_rr,
      conf_low = exp(log(rr) - zcrit * se_log_rr),
      conf_high = exp(log(rr) + zcrit * se_log_rr),
      p_value = 2 * stats::pnorm(-abs(log(rr) / se_log_rr)),
      stringsAsFactors = FALSE
    )
  }
  list(
    effects = do.call(rbind, rows), fit = fit,
    risks = stats::setNames(risks, strategy_levels),
    robust_vcov = vcov, warning = NULL
  )
}

.tt_censor_counts <- function(clones) {
  tab <- table(
    factor(stats::na.omit(clones$censor_reason),
           levels = c("administrative", "outcome", "loss_to_follow_up",
                      "artificial_deviation")))
  data.frame(
    reason = names(tab), count = as.integer(tab),
    stringsAsFactors = FALSE
  )
}

.tt_strategy_balance <- function(
    endpoint, baseline, columns, strategy_levels) {
  if (!length(columns$baseline)) {
    return(data.frame(
      contrast = character(), covariate = character(),
      level = character(), smd_before = numeric(),
      smd_after = numeric(), stringsAsFactors = FALSE
    ))
  }
  baseline_key <- as.character(baseline[[columns$id]])
  endpoint_covariates <- baseline[
    match(endpoint$original_id, baseline_key),
    columns$baseline, drop = FALSE]
  rows <- list()
  for (j in 2:length(strategy_levels)) {
    keep <- endpoint$strategy %in% c(strategy_levels[1L],
                                    strategy_levels[j])
    a <- as.integer(endpoint$strategy[keep] == strategy_levels[j])
    covariates <- endpoint_covariates[keep, , drop = FALSE]
    before <- .tt_balance(a, covariates)
    after <- .tt_balance(
      a, covariates, endpoint$analysis_weight[keep])
    table <- merge(
      before, after, by = c("covariate", "level"),
      suffixes = c("_before", "_after"), sort = FALSE)
    table$contrast <- paste(strategy_levels[j], "vs", strategy_levels[1L])
    rows[[j - 1L]] <- table[
      , c("contrast", "covariate", "level",
          "smd_before", "smd_after"), drop = FALSE]
  }
  do.call(rbind, rows)
}

#' Emulate a declared target trial
#'
#' Applies eligibility at the declared time zero and estimates either an
#' intention-to-treat contrast or a per-protocol clone-censor-weight contrast.
#' The result remains conditional on consistency, exchangeability, positivity,
#' correct treatment/censoring/outcome models, and a correctly aligned time
#' zero. A Cox hazard ratio is non-collapsible and is not a marginal risk ratio.
#'
#' @param data Long person-period data.
#' @param protocol A [targetTrialProtocol()] object.
#' @param id,time,treatment,outcome Column names. For survival outcomes,
#'   `outcome` is the subject-level outcome/censoring time repeated on each long
#'   row and `event` is its repeated status. Without `event`, `outcome` is the
#'   repeated fixed-horizon binary outcome.
#' @param event Optional binary event-status column.
#' @param baseline_covariates,time_varying_covariates Declared predictors for
#'   treatment/adherence models.
#' @param estimand `"per_protocol"` or `"intention_to_treat"`.
#' @param weight_backend Propensity backend. Optional backends are validation
#'   paths and never silent fallbacks.
#' @param stabilized Use time-based numerator probabilities.
#' @param numerator_covariates Additional numerator-model covariates.
#' @param weight_truncation Lower and upper quantiles applied after raw
#'   cumulative weights are retained.
#' @param max_weight Additional positive upper cap applied after quantile
#'   truncation.
#' @param conf.level Confidence level.
#' @return An `AnalysisResult` containing effects, models, clone rows, both raw
#'   and analysis weights, diagnostics, censoring counts, and provenance.
#' @export
targetTrialEmulate <- function(
    data,
    protocol,
    id,
    time,
    treatment,
    outcome,
    event = NULL,
    baseline_covariates = NULL,
    time_varying_covariates = NULL,
    estimand = c("per_protocol", "intention_to_treat"),
    weight_backend = c("internal", "WeightIt", "ipw"),
    stabilized = TRUE,
    numerator_covariates = NULL,
    weight_truncation = c(0.01, 0.99),
    max_weight = Inf,
    conf.level = 0.95) {
  estimand <- .tt_match_arg(
    estimand, c("per_protocol", "intention_to_treat"), "estimand")
  weight_backend <- .tt_match_arg(
    weight_backend, c("internal", "WeightIt", "ipw"), "weight_backend")
  if (!is.logical(stabilized) || length(stabilized) != 1L ||
      is.na(stabilized)) {
    stop("'stabilized' must be one logical value.", call. = FALSE)
  }
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("'conf.level' must be finite and strictly between 0 and 1.",
         call. = FALSE)
  }
  prepared <- .tt_prepare_data(
    data, protocol, id, time, treatment, outcome, event,
    baseline_covariates, time_varying_covariates)
  columns <- prepared$columns
  horizon <- if (is.numeric(protocol$follow_up) &&
                 length(protocol$follow_up) == 1L &&
                 is.finite(protocol$follow_up)) {
    protocol$follow_up
  } else {
    max(prepared$data[[columns$outcome]])
  }
  if (horizon <= prepared$time_zero) {
    stop("Protocol follow-up must occur after time zero.", call. = FALSE)
  }
  warnings <- character()

  if (estimand == "per_protocol") {
    probabilities <- .tt_fit_interval_probabilities(
      prepared$data, columns, weight_backend, stabilized,
      numerator_covariates)
    clone_data <- .tt_clone(
      prepared$data, protocol, columns, probabilities, horizon)
    truncated <- .tt_apply_truncation(
      clone_data, weight_truncation, max_weight)
    clone_data <- truncated$data
    endpoint <- .tt_clone_endpoint(clone_data, columns, horizon)
    censor_models <- list(
      denominator = probabilities$denominator$fit,
      numerator = if (is.null(probabilities$numerator)) NULL else
        probabilities$numerator$fit,
      formulas = probabilities$formulas
    )
    censoring_counts <- .tt_censor_counts(clone_data)
    positivity <- list(
      interval_probability = c(
        min = min(clone_data$denominator_probability, na.rm = TRUE),
        max = max(clone_data$denominator_probability, na.rm = TRUE)
      ),
      raw_weight = stats::quantile(
        clone_data$raw_weight[is.finite(clone_data$raw_weight)],
        c(0, 0.01, 0.5, 0.99, 1), names = TRUE),
      analysis_weight = stats::quantile(
        clone_data$analysis_weight[is.finite(clone_data$analysis_weight)],
        c(0, 0.01, 0.5, 0.99, 1), names = TRUE),
      truncation = truncated$diagnostics
    )
    balance <- list(
      baseline = .tt_strategy_balance(
        endpoint, prepared$baseline, columns,
        names(protocol$treatment_strategies)),
      effective_sample_size_by_strategy = vapply(
        split(endpoint$analysis_weight, endpoint$strategy),
        .tt_effective_n, numeric(1))
    )
  } else {
    clone_data <- NULL
    censor_models <- NULL
    censoring_counts <- data.frame(
      reason = c("administrative", "outcome", "loss_to_follow_up",
                 "artificial_deviation"),
      count = c(0L, 0L, 0L, 0L), stringsAsFactors = FALSE
    )
    baseline_cov <- prepared$baseline[, columns$baseline, drop = FALSE]
    if (ncol(baseline_cov)) {
      baseline_iptw <- .baseline_iptw(
        prepared$baseline[[columns$treatment]], baseline_cov,
        estimand = "ATE", stabilized = stabilized,
        backend = weight_backend)
      raw_weights <- baseline_iptw$weights
      censor_models <- list(baseline = baseline_iptw$fit,
                            formula = baseline_iptw$formula)
    } else {
      raw_weights <- rep(1, nrow(prepared$baseline))
    }
    weight_audit <- data.frame(
      clone_id = as.character(prepared$baseline[[columns$id]]),
      .source_row = seq_len(nrow(prepared$baseline)),
      adherent = TRUE,
      raw_weight = raw_weights,
      stringsAsFactors = FALSE)
    truncated <- .tt_apply_truncation(
      weight_audit, weight_truncation, max_weight)
    weights <- truncated$data$analysis_weight
    endpoint <- .tt_itt_endpoint(
      prepared, protocol, horizon, weights)
    endpoint$raw_weight <- raw_weights
    baseline_code <- .tt_binary_code(
      prepared$baseline[[columns$treatment]], "baseline treatment")
    if (ncol(baseline_cov)) {
      positivity <- .tt_weight_diagnostics(
        baseline_code$value, baseline_iptw$propensity,
        weights, baseline_cov)
      positivity$truncation <- truncated$diagnostics
      balance <- list(baseline = positivity$balance)
    } else {
      positivity <- list(
        message = "unweighted: no baseline covariates declared",
        effective_sample_size = length(weights),
        truncation = truncated$diagnostics)
      balance <- list(baseline = data.frame())
    }
  }

  strategy_levels <- names(protocol$treatment_strategies)
  usable <- if (is.null(columns$event)) {
    !is.na(endpoint$endpoint_outcome)
  } else {
    endpoint$endpoint_time > 0
  }
  missing_strategy <- setdiff(
    strategy_levels, unique(endpoint$strategy[usable]))
  if (length(missing_strategy)) {
    stop(sprintf(
      "Positivity failure: no analysable endpoint remains for strategy: %s.",
      paste(missing_strategy, collapse = ", ")), call. = FALSE)
  }
  fitted <- if (is.null(columns$event)) {
    .tt_binary_effect(endpoint, strategy_levels, conf.level)
  } else {
    .tt_survival_effect(endpoint, strategy_levels, conf.level)
  }
  if (!is.null(fitted$warning)) {
    warnings <- c(warnings, fitted$warning)
  }
  propensity_range <- if (estimand == "per_protocol") {
    positivity$interval_probability
  } else if (!is.null(positivity$propensity)) {
    positivity$propensity
  } else {
    NULL
  }
  if (length(propensity_range) &&
      (min(propensity_range) < 0.01 || max(propensity_range) > 0.99)) {
    warnings <- c(
      warnings,
      "Near-positivity violation: fitted probabilities extend outside [0.01, 0.99]."
    )
  }
  effects <- fitted$effects
  estimates <- stats::setNames(
    effects$estimate,
    paste(effects$contrast, effects$measure, sep = ":"))
  assumptions <- c(
    "consistency",
    "conditional exchangeability",
    "positivity",
    "correct treatment, censoring, and outcome model specification",
    "correctly aligned time zero"
  )
  backend_versions <- c(
    survival = if (is.null(columns$event)) NA_character_ else
      as.character(utils::packageVersion("survival")),
    weighting = if (weight_backend == "internal") {
      as.character(getRversion())
    } else {
      as.character(utils::packageVersion(weight_backend))
    }
  )
  analysis_hash <- .causal_hash(list(
    protocol_hash = protocol$protocol_hash,
    columns = columns,
    estimand = estimand,
    weight_formulas = if (is.null(censor_models)) {
      NULL
    } else if (!is.null(censor_models$formulas)) {
      censor_models$formulas
    } else {
      censor_models$formula
    },
    strategy_representations = lapply(
      protocol$strategy_specs,
      function(x) x[c("type", "representation", "hash")])
  ))
  interpretation <- if (is.null(columns$event)) {
    "Risk differences and risk ratios are standardized fixed-horizon contrasts."
  } else {
    paste(
      "Cox hazard ratios are non-collapsible, conditional rate contrasts;",
      "they are not marginal risk ratios."
    )
  }

  PhysioCore::AnalysisResult(
    type = "target_trial",
    estimate = estimates,
    uncertainty = list(
      type = "analytic", level = conf.level,
      lower = stats::setNames(
        effects$conf_low, paste(effects$contrast, effects$measure, sep = ":")),
      upper = stats::setNames(
        effects$conf_high, paste(effects$contrast, effects$measure, sep = ":"))
    ),
    method = if (estimand == "per_protocol") {
      "target-trial clone-censor-weight emulation"
    } else {
      "target-trial intention-to-treat emulation"
    },
    estimand = list(
      treatment = paste(strategy_levels, collapse = " versus "),
      population = "participants eligible at the declared time zero",
      endpoint = protocol$outcome,
      causal_contrast = protocol$causal_contrast,
      analysis = estimand,
      assumptions = assumptions,
      interpretation = interpretation
    ),
    result = list(
      protocol = protocol,
      effects = effects,
      fitted_outcome_model = fitted$fit,
      standardized_risks = fitted$risks,
      robust_vcov = fitted$robust_vcov,
      endpoint_data = endpoint,
      clone_data = clone_data,
      raw_weights = if (estimand == "per_protocol") {
        clone_data$raw_weight
      } else {
        endpoint$raw_weight
      },
      analysis_weights = if (estimand == "per_protocol") {
        clone_data$analysis_weight
      } else {
        endpoint$analysis_weight
      },
      censor_models = censor_models,
      positivity = positivity,
      balance = balance,
      censoring_counts = censoring_counts,
      eligibility = prepared$counts,
      warnings = warnings,
      interpretation = interpretation,
      analysis_hash = analysis_hash,
      assumptions = assumptions
    ),
    parameters = list(
      estimand = estimand,
      weight_backend = weight_backend,
      stabilized = stabilized,
      weight_truncation = weight_truncation,
      max_weight = max_weight,
      conf.level = conf.level,
      treatment_coding = prepared$treatment_levels,
      time_zero = prepared$time_zero,
      horizon = horizon
    ),
    provenance = data.frame(
      step = "targetTrialEmulate",
      protocol_id = protocol$protocol_id,
      protocol_hash = protocol$protocol_hash,
      analysis_hash = analysis_hash,
      weight_backend = weight_backend,
      weight_backend_version = backend_versions[["weighting"]],
      outcome_backend = if (is.null(columns$event)) "stats::glm" else
        "survival::coxph",
      outcome_backend_version = backend_versions[["survival"]],
      timestamp = NA_character_,
      stringsAsFactors = FALSE
    )
  )
}
