# Therapy dose-response modelling. Relates a continuous / ordinal rehabilitation
# dose (hours, repetitions, sessions) to an outcome (or its change from
# baseline), with confounder adjustment and optional non-linearity. Every
# coefficient (slope, ED50, Emax) is ESTIMATED FROM THE SUPPLIED DATA - nothing
# is hard-coded; the estimate is associational unless the dose was randomised
# (see the package's estimand tooling for causal framing).

# Tidy an lm coefficient table.
.dr_lm_coef <- function(fit) {
  co <- summary(fit)$coefficients
  data.frame(term = rownames(co), estimate = co[, 1], se = co[, 2],
             statistic = co[, 3], p = co[, 4], row.names = NULL,
             stringsAsFactors = FALSE)
}

# Covariate value for the prediction grid: mean (numeric) or reference/modal
# level (factor/character).
.dr_ref_value <- function(col) {
  if (is.numeric(col)) return(mean(col, na.rm = TRUE))
  if (is.factor(col)) return(levels(col)[1])
  names(sort(table(col), decreasing = TRUE))[1]
}

#' Therapy dose-response model
#'
#' Fits an outcome (or change from baseline) as a function of a continuous /
#' ordinal therapy dose, with optional confounder adjustment, subject random
#' intercepts and a non-linear dose form. Returns the fitted dose effect, a
#' predicted dose-response curve and the coefficient table.
#'
#' @param data A `data.frame` with the columns named below.
#' @param outcome Name of the outcome column.
#' @param dose Name of the therapy-dose column (hours / repetitions / sessions).
#' @param covariates Optional character vector of confounder columns to adjust
#'   for (added linearly).
#' @param subject Optional subject-id column; when given, a random intercept per
#'   subject is fit via [fitMixedModel()] (needs `lme4`).
#' @param baseline Optional baseline column; when given it is added as an
#'   adjustment covariate (ANCOVA-style change analysis).
#' @param form Dose form: `"linear"`, `"log"` (needs positive dose),
#'   `"spline"` (natural cubic spline, `df` terms) or `"emax"`
#'   (`E0 + Emax*dose/(ED50+dose)`, fit by [stats::nls()]; no subject/covariates).
#' @param df Spline degrees of freedom for `form = "spline"` (default 3).
#' @param n_grid Number of dose points in the predicted curve (default 50).
#' @return A [PhysioCore::AnalysisResult] of `type = "dose_response"`: `estimate`
#'   is the dose effect (slope, spline coefficients, or `E0`/`Emax`/`ED50`);
#'   `result` holds `coefficients`, the `curve` (`dose`, `predicted`), the `fit`
#'   and the `form`.
#' @references Ruberg SJ (1995) Dose response studies. *J Biopharm Stat* 5:1-14.
#' @seealso [fitMixedModel()], [proportionalRecoveryRule()]
#' @importFrom stats as.formula lm nls predict median sd
#' @importFrom splines ns
#' @export
#' @examples
#' set.seed(1)
#' d <- data.frame(dose = rep(c(0, 5, 10, 20), each = 15))
#' d$change <- 0.4 * d$dose + rnorm(nrow(d), 0, 2)
#' dr <- doseResponse(d, outcome = "change", dose = "dose", form = "linear")
#' PhysioCore::resultValue(dr)$coefficients
doseResponse <- function(data, outcome, dose, covariates = NULL, subject = NULL,
                         baseline = NULL, form = c("linear", "log", "spline",
                                                   "emax"), df = 3, n_grid = 50) {
  form <- match.arg(form)
  stopifnot(is.data.frame(data), outcome %in% names(data), dose %in% names(data))
  d <- data
  d[[outcome]] <- as.numeric(d[[outcome]])
  d[[dose]] <- as.numeric(d[[dose]])
  if (!is.null(baseline)) {
    stopifnot(baseline %in% names(data))
    covariates <- unique(c(covariates, baseline))
  }
  if (!is.null(covariates)) stopifnot(all(covariates %in% names(data)))
  keep <- stats::complete.cases(d[, unique(c(outcome, dose, covariates, subject)),
                                   drop = FALSE])
  d <- d[keep, , drop = FALSE]
  if (nrow(d) < 4L) stop("too few complete rows to fit a dose-response.", call. = FALSE)

  grid <- seq(min(d[[dose]]), max(d[[dose]]), length.out = n_grid)
  newdata <- data.frame(x = grid); names(newdata) <- dose
  for (cv in covariates) newdata[[cv]] <- .dr_ref_value(d[[cv]])

  if (form == "emax") {
    if (!is.null(subject) || !is.null(covariates)) {
      warning("the 'emax' form ignores 'subject'/'covariates' (fixed nls only).",
              call. = FALSE)
    }
    y <- d[[outcome]]; x <- d[[dose]]
    start <- list(E0 = min(y), Emax = diff(range(y)),
                  ED50 = max(stats::median(x[x > 0]), .Machine$double.eps))
    fit <- tryCatch(
      stats::nls(stats::as.formula(sprintf("%s ~ E0 + Emax * %s/(ED50 + %s)",
                                           outcome, dose, dose)),
                 data = d, start = start),
      error = function(e) stop("emax nls did not converge: ",
                               conditionMessage(e), call. = FALSE))
    co <- summary(fit)$coefficients
    coef_df <- data.frame(term = rownames(co), estimate = co[, 1], se = co[, 2],
                          statistic = co[, 3], p = co[, 4], row.names = NULL,
                          stringsAsFactors = FALSE)
    est <- stats::setNames(co[, 1], rownames(co))
    curve <- data.frame(dose = grid,
                        predicted = as.numeric(stats::predict(fit, newdata)))
  } else {
    dose_term <- switch(form,
      linear = dose,
      log = { if (any(d[[dose]] <= 0)) stop("form='log' needs positive dose.",
                                            call. = FALSE); sprintf("log(%s)", dose) },
      spline = sprintf("ns(%s, df = %d)", dose, df))
    rhs <- paste(c(dose_term, covariates), collapse = " + ")
    fml <- stats::as.formula(paste(outcome, "~", rhs))
    if (!is.null(subject)) {
      mm <- fitMixedModel(d, fml,
                          random = stats::as.formula(sprintf("~ (1 | %s)", subject)))
      coef_df <- PhysioCore::resultValue(mm)$fixed
      fit <- PhysioCore::resultValue(mm)$fit
      curve <- data.frame(dose = grid, predicted = as.numeric(
        stats::predict(fit, newdata, re.form = NA)))
    } else {
      fit <- stats::lm(fml, data = d)
      coef_df <- .dr_lm_coef(fit)
      curve <- data.frame(dose = grid,
                          predicted = as.numeric(stats::predict(fit, newdata)))
    }
    dterm <- grep(if (form == "linear") sprintf("^%s$", dose) else "dose|ns\\(",
                  coef_df$term)
    est <- stats::setNames(coef_df$estimate[dterm], coef_df$term[dterm])
  }

  PhysioCore::AnalysisResult(
    type = "dose_response", estimate = est, method = form,
    result = list(coefficients = coef_df, curve = curve, fit = fit, form = form),
    parameters = list(outcome = outcome, dose = dose, covariates = covariates,
                      subject = subject, baseline = baseline, form = form, df = df),
    provenance = data.frame(step = "doseResponse", timestamp = NA_character_,
                            stringsAsFactors = FALSE))
}
