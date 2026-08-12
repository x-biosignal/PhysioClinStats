# Population recovery-trajectory models: a nonlinear mixed-effects (NLME) fit
# with partial pooling across subjects (nlme::nlme). Each subject shares the
# population mean curve but has its own random asymptote / rate, which stabilises
# the per-subject asymptote and time-to-90% estimates versus independent NLS.

#' Fit a population recovery-trajectory NLME
#'
#' Fits a nonlinear mixed-effects recovery curve to longitudinal panel data with
#' subject-level random effects (partial pooling), and returns the population
#' fixed effects together with per-subject predicted asymptote, rate, and
#' time-to-90%-recovery.
#'
#' Models: \code{"exponential"}/\code{"asymptotic"} use the self-starting
#' asymptotic curve \eqn{y = A + (R_0 - A)\,e^{-e^{lrc} t}} (\code{SSasymp});
#' \code{"logistic"} uses \code{SSlogis}.
#'
#' @param data A long-format data frame.
#' @param subject,time,outcome Column names (character) for the grouping factor,
#'   the time variable, and the response.
#' @param model \code{"exponential"} (default), \code{"asymptotic"} (an alias),
#'   or \code{"logistic"}.
#' @return An \code{AnalysisResult} (type \code{"recovery_trajectory"}) whose
#'   \code{result} holds the \code{fixed_effects}, \code{random_effects}, a
#'   per-subject data frame (\code{asymptote}, \code{rate}, \code{time_to_90}),
#'   and the fitted \code{nlme} object.
#' @references Pinheiro & Bates 2000 (nlme); Lindstrom & Bates 1990.
#' @seealso [proportionalRecoveryRule()], [latentClassGrowth()]
#' @export
#' @examples
#' set.seed(1)
#' df <- do.call(rbind, lapply(1:12, function(s) {
#'   A <- 60 + rnorm(1, 0, 5); rate <- 0.3 * exp(rnorm(1, 0, 0.2))
#'   t <- 0:8; data.frame(subject = s, time = t,
#'     y = A * (1 - exp(-rate * t)) + rnorm(9, 0, 2))
#' }))
#' fit <- recoveryTrajectoryLME(df, "subject", "time", "y")
recoveryTrajectoryLME <- function(data, subject, time, outcome,
                                  model = c("exponential", "asymptotic",
                                            "logistic")) {
  model <- match.arg(model)
  if (!requireNamespace("nlme", quietly = TRUE)) {
    stop("recoveryTrajectoryLME() requires the 'nlme' package.", call. = FALSE)
  }
  for (col in c(subject, time, outcome)) {
    if (!col %in% names(data)) stop("column '", col, "' not in data.",
                                    call. = FALSE)
  }
  d <- data.frame(.subj = factor(data[[subject]]),
                  .t = as.numeric(data[[time]]),
                  .y = as.numeric(data[[outcome]]))
  d <- d[stats::complete.cases(d), , drop = FALSE]
  if (nlevels(d$.subj) < 2L) {
    stop("recoveryTrajectoryLME() needs at least 2 subjects.", call. = FALSE)
  }
  gd <- nlme::groupedData(.y ~ .t | .subj, data = d)

  if (model == "logistic") {
    pars <- c("Asym", "xmid", "scal")
    form <- .y ~ SSlogis(.t, Asym, xmid, scal)
  } else {
    pars <- c("Asym", "R0", "lrc")
    form <- .y ~ SSasymp(.t, Asym, R0, lrc)
  }
  fixed_f <- stats::as.formula(paste(paste(pars, collapse = " + "), "~ 1"))
  rand_f <- stats::as.formula(paste(paste(pars, collapse = " + "), "~ 1 | .subj"))

  # full random structure, then progressively simpler on convergence failure
  fit <- tryCatch(
    nlme::nlme(form, data = gd, fixed = fixed_f, random = rand_f,
               control = nlme::nlmeControl(maxIter = 200, returnObject = TRUE)),
    error = function(e) NULL)
  if (is.null(fit)) {
    rand1 <- stats::as.formula(paste(pars[1], "~ 1 | .subj"))
    fit <- tryCatch(
      nlme::nlme(form, data = gd, fixed = fixed_f, random = rand1,
                 control = nlme::nlmeControl(maxIter = 200, returnObject = TRUE)),
      error = function(e) NULL)
  }
  if (is.null(fit)) {
    stop("the NLME fit did not converge for this data.", call. = FALSE)
  }

  fe <- nlme::fixef(fit)
  cf <- stats::coef(fit)                                  # per-subject params
  per_subj <- .recovery_per_subject(cf, model)

  res <- list(
    fixed_effects = fe,
    random_effects = as.data.frame(nlme::ranef(fit)),
    per_subject = per_subj,
    model = model, converged = TRUE, fit = fit,
    logLik = as.numeric(stats::logLik(fit)), AIC = stats::AIC(fit))

  PhysioCore::AnalysisResult(
    type = "recovery_trajectory", result = res,
    parameters = list(model = model, n_subjects = nlevels(d$.subj)),
    estimate = fe, method = "nlme")
}

# per-subject asymptote / rate / time-to-90% from the fitted coefficients
.recovery_per_subject <- function(cf, model) {
  subj <- rownames(cf)
  if (model == "logistic") {
    Asym <- cf[["Asym"]]; scal <- cf[["scal"]]; xmid <- cf[["xmid"]]
    rate <- 1 / scal                                     # logistic growth rate
    # y = 0.9*Asym  =>  t = xmid - scal*log(1/0.9 - 1)
    t90 <- xmid - scal * log(1 / 0.9 - 1)
  } else {
    Asym <- cf[["Asym"]]; lrc <- cf[["lrc"]]
    rate <- exp(lrc)                                     # rate constant
    t90 <- log(10) / rate                                # 90% of A - R0 recovery
  }
  data.frame(subject = subj, asymptote = Asym, rate = rate,
             time_to_90 = t90, row.names = NULL, stringsAsFactors = FALSE)
}
