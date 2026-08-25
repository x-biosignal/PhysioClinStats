# Distribution-free predictive uncertainty via conformal prediction
# (Vovk et al.; Lei et al. 2018 split conformal; Barber et al. 2021 jackknife+).

#' Build the provenance log for an uncertainty result
#' @keywords internal
#' @noRd
.uncertaintyProvenance <- function(step, seed, extra = list()) {
  base <- list(step = step, timestamp = NA_character_,
               seed = if (is.null(seed)) NA_integer_ else as.integer(seed))
  do.call(data.frame, c(base, extra, stringsAsFactors = FALSE))
}

#' Conformal prediction interval for a fitted model
#'
#' Distribution-free predictive intervals with finite-sample marginal coverage
#' at least `1 - alpha`, requiring only exchangeability. `"split"` conformal
#' (Lei et al. 2018) scores an independent calibration set with the fitted
#' model and takes the conformal quantile of the absolute residuals;
#' `"jackknife_plus"` (Barber et al. 2021) refits the model leave-one-out over
#' the calibration data and combines the leave-one-out predictions and
#' residuals.
#'
#' The coverage guarantee requires the calibration data to be exchangeable with,
#' and (for `"split"`) independent of, the data the model was trained on. When
#' the model uses a transformed response (e.g. `log(y) ~ x`) the interval is on
#' that transformed scale.
#'
#' @param model A fitted model supporting `predict(model, newdata=)` (and, for
#'   `"jackknife_plus"`, `update(model, data=)`), e.g. an `lm`.
#' @param calib A data frame with the response and predictors: an independent
#'   calibration set for `"split"`, or the full data to refit over for
#'   `"jackknife_plus"`.
#' @param newx A data frame of predictor values to predict for.
#' @param alpha Miscoverage level (default 0.1 for 90% intervals).
#' @param type `"split"` (default) or `"jackknife_plus"`.
#' @param seed Optional integer random seed recorded in the provenance.
#'
#' @return An `AnalysisResult` (from PhysioCore) with `type = "conformal"` whose
#'   `result` holds `point`, `lower`, `upper` (one per `newx` row), the
#'   `alpha`/`method`/`quantile`, and a provenance log capturing the seed.
#' @references Lei, J. et al. (2018). JASA 113(523):1094-1111. Barber, R.F. et
#'   al. (2021). Ann Statist 49(1):486-507.
#' @seealso [conformalPredict()], [bayesEstimate()]
#' @export
#' @examples
#' set.seed(1)
#' train <- data.frame(x = rnorm(60)); train$y <- 2 * train$x + rnorm(60)
#' calib <- data.frame(x = rnorm(60)); calib$y <- 2 * calib$x + rnorm(60)
#' fit <- lm(y ~ x, data = train)
#' res <- conformalInterval(fit, calib, data.frame(x = c(-1, 0, 1)))
#' PhysioCore::resultValue(res)$lower
conformalInterval <- function(model, calib, newx, alpha = 0.1,
                              type = c("split", "jackknife_plus"), seed = NULL) {
  type <- match.arg(type)
  stopifnot(is.numeric(alpha), length(alpha) == 1L, alpha > 0, alpha < 1)
  if (!is.data.frame(calib) || !is.data.frame(newx)) {
    stop("`calib` and `newx` must be data frames.", call. = FALSE)
  }
  resp <- all.vars(stats::formula(model))[1]
  if (!resp %in% names(calib)) {
    stop(sprintf("the response '%s' is not a column of `calib`.", resp),
         call. = FALSE)
  }
  # complete calibration rows only: a missing value would be dropped from the
  # response (na.omit) but padded in prediction (na.pass), mis-aligning the
  # scores and silently voiding the coverage guarantee
  model_vars <- intersect(all.vars(stats::formula(model)), names(calib))
  if (anyNA(calib[, model_vars, drop = FALSE])) {
    stop("`calib` has missing values in the model variables; conformal ",
         "scoring requires complete calibration rows.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(as.integer(seed))
  # response on the MODEL'S scale (so a transformed response such as log(y) is
  # scored against log-scale predictions, not the raw response)
  y <- as.numeric(stats::model.response(
    stats::model.frame(stats::formula(model), data = calib)))
  n <- length(y)
  m <- nrow(newx)
  if (n < 2L) stop("`calib` needs at least 2 rows.", call. = FALSE)

  if (type == "split") {
    scores <- abs(y - as.numeric(stats::predict(model, newdata = calib)))
    k <- ceiling((n + 1) * (1 - alpha))
    qhat <- if (k > n) Inf else sort(scores)[k]
    point <- as.numeric(stats::predict(model, newdata = newx))
    lower <- point - qhat
    upper <- point + qhat
  } else {
    loo_pred <- matrix(NA_real_, n, m)     # leave-one-out prediction at newx
    loo_res <- numeric(n)                  # leave-one-out absolute residual
    for (i in seq_len(n)) {
      fit_i <- tryCatch(
        stats::update(model, data = calib[-i, , drop = FALSE]),
        error = function(e) stop("jackknife_plus could not refit the model ",
          "leave-one-out via update(); the model call likely references ",
          "variables (weights/subset/offset) not in `calib`. Original error: ",
          conditionMessage(e), call. = FALSE))
      loo_res[i] <- abs(y[i] -
        as.numeric(stats::predict(fit_i, newdata = calib[i, , drop = FALSE])))
      loo_pred[i, ] <- as.numeric(stats::predict(fit_i, newdata = newx))
    }
    lo_idx <- floor(alpha * (n + 1))
    hi_idx <- ceiling((1 - alpha) * (n + 1))
    lower <- upper <- numeric(m)
    for (j in seq_len(m)) {
      lo_vals <- sort(loo_pred[, j] - loo_res)
      hi_vals <- sort(loo_pred[, j] + loo_res)
      lower[j] <- if (lo_idx < 1L) -Inf else lo_vals[lo_idx]
      upper[j] <- if (hi_idx > n) Inf else hi_vals[hi_idx]
    }
    point <- colMeans(loo_pred)
    qhat <- NA_real_
  }

  PhysioCore::AnalysisResult(
    type = "conformal",
    result = list(point = point, lower = lower, upper = upper,
                  alpha = alpha, method = type, quantile = qhat),
    parameters = list(alpha = alpha, type = type, response = resp,
                      n_calib = n),
    provenance = .uncertaintyProvenance("conformalInterval", seed,
                                        list(type = type, alpha = alpha,
                                             n_calib = n)))
}

#' Conformal prediction band for a single new patient
#'
#' Convenience wrapper around split-conformal [conformalInterval()] for one new
#' observation, returning the guaranteed-coverage predicted band.
#'
#' @param model A fitted model (see [conformalInterval()]).
#' @param calib The calibration data frame (response + predictors).
#' @param newx A one-row data frame of predictors for the new patient.
#' @param alpha Miscoverage level (default 0.1).
#' @param seed Optional integer random seed.
#' @return An `AnalysisResult` with `type = "conformal"` for the single patient.
#' @seealso [conformalInterval()]
#' @export
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rnorm(80)); d$y <- 1.5 * d$x + rnorm(80)
#' fit <- lm(y ~ x, data = d[1:40, ])
#' conformalPredict(fit, d[41:80, ], data.frame(x = 0.5))
conformalPredict <- function(model, calib, newx, alpha = 0.1, seed = NULL) {
  if (!is.data.frame(newx) || nrow(newx) != 1L) {
    stop("`newx` must be a one-row data frame (a single new patient).",
         call. = FALSE)
  }
  conformalInterval(model, calib, newx, alpha = alpha, type = "split",
                    seed = seed)
}
