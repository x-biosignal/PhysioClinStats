# Kaplan-Meier estimation, Cox proportional-hazards modelling (with a
# proportional-hazards check) and milestone hazard ratios, wrapping the survival
# package.

.need_survival <- function(fn) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop(sprintf("%s needs the survival package.", fn), call. = FALSE)
  }
}

#' Kaplan-Meier survival fit with attainment-time summary
#'
#' Fits a Kaplan-Meier estimator with `survival::survfit`, optionally stratified
#' by a grouping variable, and returns the median and quartile attainment times
#' with their confidence intervals.
#'
#' @param data A data.frame with the survival columns.
#' @param time,event Column names of the follow-up time and the event indicator
#'   (1 = event, 0 = censored).
#' @param group Optional grouping-variable column name for stratified curves.
#' @param conf_level Confidence level (default 0.95).
#' @return An object of class `physio_km` (a list) with `fit` (the `survfit`),
#'   `quantiles` (median / Q1 / Q3 attainment times with CIs, per group) and
#'   `summary` table, plus `print` and `plot` methods.
#' @references Kaplan EL, Meier P (1958). Nonparametric estimation from
#'   incomplete observations. JASA, 53(282), 457-481.
#' @seealso [coxModel()], [timeToMilestone()]
#' @export
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   survivalFit(survival::lung, "time", "status", group = "sex")
#' }
survivalFit <- function(data, time = "time", event = "event", group = NULL,
                        conf_level = 0.95) {
  .need_survival("survivalFit")
  stopifnot(is.data.frame(data), all(c(time, event) %in% names(data)))
  rhs <- if (is.null(group)) "1" else group
  form <- stats::as.formula(sprintf("survival::Surv(%s, %s) ~ %s",
                                    time, event, rhs))
  fit <- survival::survfit(form, data = data, conf.int = conf_level)

  qtab <- stats::quantile(fit, probs = c(0.25, 0.5, 0.75), conf.int = TRUE)
  # normalise the median/quartile output into a tidy per-group data.frame
  q <- qtab$quantile; lo <- qtab$lower; hi <- qtab$upper
  if (is.null(dim(q))) {
    quant <- data.frame(group = if (is.null(group)) "all" else NA,
                        Q1 = q[["25"]], median = q[["50"]], Q3 = q[["75"]],
                        median_lower = lo[["50"]], median_upper = hi[["50"]],
                        stringsAsFactors = FALSE)
  } else {
    quant <- data.frame(group = rownames(q),
                        Q1 = q[, "25"], median = q[, "50"], Q3 = q[, "75"],
                        median_lower = lo[, "50"], median_upper = hi[, "50"],
                        row.names = NULL, stringsAsFactors = FALSE)
  }

  # numbers-at-risk table on an evenly-spaced time grid, for the risk table
  grid <- pretty(c(0, max(fit$time, na.rm = TRUE)), n = 6)
  grid <- grid[grid >= 0 & grid <= max(fit$time, na.rm = TRUE)]
  sm <- summary(fit, times = grid, extend = TRUE)
  strata_lab <- if (!is.null(sm$strata)) as.character(sm$strata) else
    if (is.null(group)) "all" else NA
  risk_table <- data.frame(time = sm$time, group = strata_lab,
                           n_risk = sm$n.risk, stringsAsFactors = FALSE)

  structure(list(fit = fit, quantiles = quant, risk_table = risk_table,
                 group = group, conf_level = conf_level, formula = form),
            class = "physio_km")
}

#' @export
print.physio_km <- function(x, ...) {
  cat("Kaplan-Meier fit\n")
  print(x$fit)
  cat("\nAttainment times (median / quartiles):\n")
  q <- x$quantiles
  num <- vapply(q, is.numeric, logical(1)); q[num] <- lapply(q[num], round, 2)
  print(q, row.names = FALSE)
  invisible(x)
}

#' @export
plot.physio_km <- function(x, risk_table = TRUE, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot.physio_km needs the ggplot2 package.", call. = FALSE)
  }
  fit <- x$fit
  strata <- if (!is.null(fit$strata)) rep(names(fit$strata), fit$strata) else
    rep("all", length(fit$time))
  df <- data.frame(time = fit$time, surv = fit$surv,
                   lower = fit$lower, upper = fit$upper,
                   strata = strata, stringsAsFactors = FALSE)
  # begin each curve at (0, 1)
  starts <- do.call(rbind, lapply(unique(df$strata), function(s)
    data.frame(time = 0, surv = 1, lower = 1, upper = 1, strata = s)))
  df <- rbind(starts, df)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$time, y = .data$surv,
                                        colour = .data$strata)) +
    ggplot2::geom_step() +
    ggplot2::labs(x = "Time", y = "Survival probability", colour = "Group",
                  title = "Kaplan-Meier estimate") +
    ggplot2::ylim(0, 1) + ggplot2::theme_minimal()

  # Compose a numbers-at-risk table beneath the curves when patchwork is
  # available; otherwise return the curves alone (the table is in x$risk_table).
  if (isTRUE(risk_table) && !is.null(x$risk_table) &&
      requireNamespace("patchwork", quietly = TRUE)) {
    rt <- x$risk_table
    tbl <- ggplot2::ggplot(rt, ggplot2::aes(x = .data$time, y = .data$group,
                                            colour = .data$group)) +
      ggplot2::geom_text(ggplot2::aes(label = .data$n_risk), size = 3) +
      ggplot2::labs(x = "Time", y = NULL, title = "Number at risk") +
      ggplot2::guides(colour = "none") + ggplot2::theme_minimal()
    return(patchwork::wrap_plots(p, tbl, ncol = 1,
                                 heights = c(4, 1)))
  }
  p
}

#' Cox proportional-hazards model with a PH check
#'
#' Fits a Cox proportional-hazards model with `survival::coxph`, tests the
#' proportional-hazards assumption with `survival::cox.zph`, and warns when the
#' global test is significant. Returns hazard ratios with confidence intervals.
#'
#' @param data A data.frame.
#' @param formula A model formula whose response is `survival::Surv(time,
#'   event)`, or supply `time`/`event`/`covariates` instead.
#' @param time,event,covariates Alternative to `formula`: the time and event
#'   column names and a character vector of covariate columns.
#' @param conf_level Confidence level (default 0.95).
#' @param ph_alpha Significance level for the global PH warning (default 0.05).
#' @return An `AnalysisResult` (`type = "cox_model"`) whose `estimate` is the
#'   hazard-ratio vector, with `result$coefficients` (coef, HR, SE, z, p, CI),
#'   `result$ph_test` (the `cox.zph` table), `result$ph_violated` and
#'   `result$fit`.
#' @references Cox DR (1972). Regression models and life-tables. JRSS B, 34(2).
#'   Grambsch PM, Therneau TM (1994). Proportional hazards tests and diagnostics
#'   based on weighted residuals. Biometrika, 81(3).
#' @seealso [survivalFit()], [milestoneHazard()]
#' @export
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   coxModel(survival::lung, time = "time", event = "status",
#'            covariates = c("age", "sex", "ph.ecog"))
#' }
coxModel <- function(data, formula = NULL, time = "time", event = "event",
                     covariates = NULL, conf_level = 0.95, ph_alpha = 0.05) {
  .need_survival("coxModel")
  stopifnot(is.data.frame(data))
  if (is.null(formula)) {
    if (is.null(covariates)) {
      stop("Supply either 'formula' or 'covariates'.", call. = FALSE)
    }
    formula <- stats::as.formula(sprintf(
      "survival::Surv(%s, %s) ~ %s", time, event,
      paste(covariates, collapse = " + ")))
  }
  # model = TRUE stores the model frame in the fit so cox.zph does not have to
  # re-evaluate the data from the (possibly foreign) formula environment.
  fit <- survival::coxph(formula, data = data, model = TRUE)
  sm <- summary(fit, conf.int = conf_level)
  cc <- sm$coefficients; ci <- sm$conf.int
  # a robust / cluster() fit adds a "robust se" column that drives z, p and the
  # CI; report that SE so the row is self-consistent.
  se_col <- if ("robust se" %in% colnames(cc)) "robust se" else "se(coef)"
  coefs <- data.frame(
    term = rownames(cc), coef = cc[, "coef"], hr = cc[, "exp(coef)"],
    std_error = cc[, se_col], z = cc[, "z"],
    p_value = cc[, ncol(cc)],
    hr_lower = ci[, 3], hr_upper = ci[, 4],
    row.names = NULL, stringsAsFactors = FALSE)

  # cox.zph fails on a zero-event fit; a milestone with no attainers is a
  # documented input, so tolerate it (like survivalFit) with an NA-HR result.
  if (isTRUE(fit$nevent == 0L)) {
    warning("The endpoint has 0 events (no subject attained the milestone); ",
            "the hazard ratio is undefined.", call. = FALSE)
    ph_tab <- data.frame()
    ph_violated <- FALSE
  } else {
    zph <- survival::cox.zph(fit)
    ph_tab <- as.data.frame(zph$table)
    ph_violated <- isTRUE(ph_tab["GLOBAL", "p"] < ph_alpha)
    if (ph_violated) {
      warning(sprintf(paste0("Proportional-hazards assumption is violated ",
                            "(global cox.zph p = %.4g < %.2g); the hazard ",
                            "ratios are time-averaged."),
                      ph_tab["GLOBAL", "p"], ph_alpha), call. = FALSE)
    }
  }

  est <- stats::setNames(coefs$hr, coefs$term)
  PhysioCore::AnalysisResult(
    type = "cox_model", estimate = est, method = "Cox proportional hazards",
    result = list(coefficients = coefs, ph_test = ph_tab,
                  ph_violated = ph_violated, fit = fit),
    parameters = list(conf_level = conf_level, ph_alpha = ph_alpha),
    provenance = data.frame(step = "coxModel", timestamp = NA_character_,
                            stringsAsFactors = FALSE))
}

#' Hazard ratio for attaining a milestone
#'
#' Fits a Cox model for the time-to-milestone endpoint against a grouping
#' variable and returns the group hazard ratio - the relative rate of attaining
#' the milestone - with its confidence interval.
#'
#' @param milestone A time-to-milestone data.frame from [timeToMilestone()]
#'   (columns `time`, `event`), joined with the grouping variable.
#' @param group Column name of the grouping variable.
#' @param conf_level Confidence level (default 0.95).
#' @return An `AnalysisResult` (`type = "milestone_hazard"`) as [coxModel()].
#' @seealso [timeToMilestone()], [coxModel()]
#' @export
milestoneHazard <- function(milestone, group, conf_level = 0.95) {
  .need_survival("milestoneHazard")
  stopifnot(is.data.frame(milestone),
            all(c("time", "event", group) %in% names(milestone)))
  res <- coxModel(milestone, time = "time", event = "event",
                  covariates = group, conf_level = conf_level)
  res@type <- "milestone_hazard"
  res
}
