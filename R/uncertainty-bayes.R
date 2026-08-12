# Bayesian per-estimate uncertainty: a thin Stan wrapper with an analytic
# Normal-Normal conjugate fallback, and credible intervals (HDI / quantile).

#' Bayesian estimate of a mean (conjugate Normal-Normal or Stan)
#'
#' For a numeric sample, returns the analytic Normal-Normal conjugate posterior
#' for the mean: a normal prior `N(prior_mean, prior_sd^2)` combined with a
#' normal likelihood of known scale `sigma` gives a normal posterior in closed
#' form (a flat prior, `prior_sd = Inf`, reproduces the sampling posterior
#' `N(mean(y), sigma^2/n)`). For a model `formula` the estimation is delegated
#' to rstanarm/brms when installed, and errors informatively otherwise -- the
#' analytic numeric-vector path always runs.
#'
#' @param y A numeric vector (analytic path) or a model `formula` (Stan path).
#' @param data Data frame for the `formula` path.
#' @param prior_mean,prior_sd Normal prior mean and SD (default 0 and `Inf`, a
#'   flat prior).
#' @param sigma Known likelihood SD; `NULL` uses the sample SD of `y`.
#' @param level Credible level (default 0.95).
#' @param seed Optional integer random seed recorded in the provenance.
#'
#' @return An `AnalysisResult` (from PhysioCore) with `type = "bayes"` whose
#'   `result` holds `posterior_mean`, `posterior_sd`, `ci_lower`, `ci_upper`,
#'   the `level`/`method`, plus a provenance log with the seed.
#' @references Gelman, A. et al. (2013). Bayesian Data Analysis, 3rd ed.
#' @seealso [credibleInterval()], [conformalInterval()]
#' @export
#' @examples
#' set.seed(1)
#' bayesEstimate(rnorm(30, mean = 5), prior_mean = 0, prior_sd = 10)
bayesEstimate <- function(y, data = NULL, prior_mean = 0, prior_sd = Inf,
                          sigma = NULL, level = 0.95, seed = NULL) {
  stopifnot(is.numeric(level), length(level) == 1L, level > 0, level < 1)
  if (!is.null(seed)) set.seed(as.integer(seed))

  if (inherits(y, "formula")) {
    have_rstanarm <- requireNamespace("rstanarm", quietly = TRUE)
    have_brms <- requireNamespace("brms", quietly = TRUE)
    if (!have_rstanarm && !have_brms) {
      stop("bayesEstimate() with a model formula needs the rstanarm or brms ",
           "package. For a mean, pass a numeric vector to use the analytic ",
           "Normal-Normal fallback.", call. = FALSE)
    }
    # dispatch to whichever backend is actually installed
    if (have_rstanarm) {
      fit <- rstanarm::stan_glm(y, data = data, seed = seed %||% 1L,
                                refresh = 0)
      draws <- as.matrix(fit)[, "(Intercept)"]
      backend <- "rstanarm"
    } else {
      fit <- brms::brm(y, data = data, seed = seed %||% 1L, refresh = 0,
                       silent = 2)
      draws <- as.data.frame(fit)[["b_Intercept"]]
      backend <- "brms"
    }
    ci <- credibleInterval(draws, level = level, method = "quantile")
    return(PhysioCore::AnalysisResult(
      type = "bayes",
      result = list(posterior_mean = mean(draws), posterior_sd = stats::sd(draws),
                    ci_lower = ci[1], ci_upper = ci[2], level = level,
                    method = "stan"),
      parameters = list(level = level, backend = backend),
      provenance = .uncertaintyProvenance("bayesEstimate", seed,
                                          list(method = "stan"))))
  }

  y <- as.numeric(y); y <- y[is.finite(y)]
  n <- length(y)
  if (n < 1L) stop("`y` must contain at least one finite value.", call. = FALSE)
  stopifnot(is.numeric(prior_mean), length(prior_mean) == 1L,
            is.finite(prior_mean))
  stopifnot(is.numeric(prior_sd), length(prior_sd) == 1L, prior_sd > 0)
  sig <- if (is.null(sigma)) stats::sd(y) else sigma
  stopifnot(is.numeric(sig), length(sig) == 1L, is.finite(sig), sig > 0)

  prior_prec <- if (is.infinite(prior_sd)) 0 else 1 / prior_sd^2
  data_prec <- n / sig^2
  post_prec <- prior_prec + data_prec
  post_var <- 1 / post_prec
  post_mean <- post_var * (prior_prec * prior_mean + data_prec * mean(y))
  z <- stats::qnorm((1 + level) / 2)
  post_sd <- sqrt(post_var)

  PhysioCore::AnalysisResult(
    type = "bayes",
    result = list(posterior_mean = post_mean, posterior_sd = post_sd,
                  ci_lower = post_mean - z * post_sd,
                  ci_upper = post_mean + z * post_sd, level = level,
                  method = "normal-normal-analytic"),
    parameters = list(level = level, prior_mean = prior_mean,
                      prior_sd = prior_sd, sigma = sig, n = n),
    provenance = .uncertaintyProvenance("bayesEstimate", seed,
                                        list(method = "analytic", n = n)))
}

#' Credible interval from posterior samples
#'
#' The highest-density interval (HDI, the narrowest interval carrying `level`
#' probability) or the equal-tailed quantile interval, from a vector of
#' posterior draws.
#'
#' @param x Numeric vector of posterior samples.
#' @param level Credible level (default 0.95).
#' @param method `"hdi"` (default, narrowest) or `"quantile"` (equal-tailed).
#' @return Numeric `c(lower, upper)`.
#' @seealso [bayesEstimate()]
#' @export
#' @examples
#' set.seed(1)
#' credibleInterval(rnorm(2000), method = "hdi")
credibleInterval <- function(x, level = 0.95, method = c("hdi", "quantile")) {
  method <- match.arg(method)
  stopifnot(is.numeric(level), length(level) == 1L, level > 0, level < 1)
  x <- as.numeric(x); x <- x[is.finite(x)]
  n <- length(x)
  if (n < 2L) stop("`x` must contain at least 2 finite samples.", call. = FALSE)
  if (method == "quantile") {
    return(unname(stats::quantile(x, c((1 - level) / 2, (1 + level) / 2))))
  }
  xs <- sort(x)
  w <- max(1L, as.integer(ceiling(level * n)))
  if (w >= n) return(c(xs[1], xs[n]))
  starts <- seq_len(n - w + 1L)
  widths <- xs[starts + w - 1L] - xs[starts]
  i <- which.min(widths)
  c(xs[i], xs[i + w - 1L])
}

`%||%` <- function(a, b) if (is.null(a)) b else a
