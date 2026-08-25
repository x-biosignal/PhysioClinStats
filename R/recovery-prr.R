# Proportional recovery rule (PRR): after stroke, many patients recover about
# 70% of their maximum potential (Prabhakaran 2008; Winters 2015). Because the
# change (followup - initial) and the potential (max - initial) both contain
# -initial, they are mathematically coupled, so a ~0.7 slope can arise as a
# regression-to-the-mean artefact even without true proportional recovery
# (Hawe/Hope 2019). This function estimates the rule AND guards that artefact.

#' Evaluate the proportional recovery rule with an artefact guard
#'
#' Tests whether recovery follows the proportional recovery rule
#' \eqn{\Delta \approx p \,(\mathrm{max} - \mathrm{initial})} (default
#' \eqn{p = 0.7}), classifies fitters vs non-fitters, handles ceiling effects,
#' and runs two diagnostics against the known statistical artefact: a
#' mixture-vs-single-line comparison and a regression-to-the-mean (label-shuffle)
#' control.
#'
#' A genuine proportional rule implies the follow-up regresses on the initial
#' score with slope \eqn{1 - p}; a pure regression-to-the-mean / coupling
#' artefact instead reproduces the \eqn{\Delta}-vs-potential slope when the
#' follow-up scores are shuffled across subjects. When the observed slope is not
#' distinguishable from that shuffled null, \code{artefact_suspected} is set.
#'
#' @param initial Numeric baseline scores.
#' @param followup Numeric follow-up scores (same length).
#' @param max_score Scalar maximum attainable score (the recovery ceiling).
#' @param expected_prop Expected recovery proportion (default 0.7).
#' @param fitter_tol Half-width around \code{expected_prop} for a "fitter"
#'   (default 0.2).
#' @param ceiling_frac Subjects with \code{initial >= ceiling_frac * max_score}
#'   are treated as ceiling cases (default 0.9) and excluded from the slope fit.
#' @param n_shuffle Label-shuffles for the regression-to-the-mean control
#'   (default 999).
#' @param seed Optional RNG seed.
#' @return An \code{AnalysisResult} (type \code{"proportional_recovery"}) whose
#'   \code{result} holds the fitted \code{slope} (+ CI and test vs
#'   \code{expected_prop}), per-subject \code{classification}, the
#'   \code{mixture_preferred} flag, and \code{artefact_suspected}.
#' @references Prabhakaran 2008; Winters 2015; Hawe, Scott & Dukelow 2019
#'   (proportional-recovery artefact).
#' @seealso [recoveryTrajectoryLME()], [latentClassGrowth()]
#' @export
#' @examples
#' set.seed(1)
#' init <- runif(60, 5, 55)
#' fu <- init + 0.7 * (66 - init) + rnorm(60, 0, 3)   # genuine 0.7 rule
#' proportionalRecoveryRule(init, fu, max_score = 66)
proportionalRecoveryRule <- function(initial, followup, max_score,
                                     expected_prop = 0.7, fitter_tol = 0.2,
                                     ceiling_frac = 0.9, n_shuffle = 999,
                                     seed = NULL) {
  initial <- as.numeric(initial); followup <- as.numeric(followup)
  if (length(initial) != length(followup)) {
    stop("'initial' and 'followup' must have the same length.", call. = FALSE)
  }
  if (length(max_score) != 1L || !is.finite(max_score)) {
    stop("'max_score' must be a single finite value.", call. = FALSE)
  }
  ok <- is.finite(initial) & is.finite(followup)
  initial <- initial[ok]; followup <- followup[ok]
  n <- length(initial)
  if (n < 5L) stop("need at least 5 complete subjects.", call. = FALSE)
  n_shuffle <- as.integer(n_shuffle)
  if (is.na(n_shuffle) || n_shuffle < 1L) {
    stop("'n_shuffle' must be a positive integer.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  potential <- max_score - initial
  delta <- followup - initial
  ceiling <- initial >= ceiling_frac * max_score          # little room to recover
  obs_prop <- ifelse(potential > 0, delta / potential, NA_real_)

  fit_idx <- which(!ceiling & potential > 0)
  if (length(fit_idx) < 3L) {
    stop("too few non-ceiling subjects to fit the rule.", call. = FALSE)
  }
  # slope of the change on the potential, through the origin (the PRR slope)
  pot <- potential[fit_idx]; dlt <- delta[fit_idx]
  lm_fit <- stats::lm(dlt ~ 0 + pot)
  slope <- unname(stats::coef(lm_fit)[1])
  ci <- stats::confint(lm_fit)[1, ]
  # test the slope against expected_prop
  se <- summary(lm_fit)$coefficients[1, 2]
  t_stat <- (slope - expected_prop) / se
  p_vs_expected <- 2 * stats::pt(-abs(t_stat), df = length(fit_idx) - 1)

  # per-subject classification
  classification <- rep("non_fitter", n)
  classification[ceiling] <- "ceiling"
  is_fitter <- !ceiling & is.finite(obs_prop) &
    abs(obs_prop - expected_prop) <= fitter_tol
  classification[is_fitter] <- "fitter"

  # --- artefact guard 1: regression-to-the-mean / coupling control ---
  # shuffling followup across subjects breaks any true initial->followup link but
  # preserves the marginal distributions and the -initial coupling; the shuffled
  # slope is the artefact baseline.
  shuffled <- vapply(seq_len(n_shuffle), function(i) {
    fu_s <- sample(followup)
    d_s <- (fu_s - initial)[fit_idx]
    unname(stats::coef(stats::lm(d_s ~ 0 + pot))[1])
  }, numeric(1))
  shuffle_ci <- stats::quantile(shuffled, c(0.025, 0.975), names = FALSE)
  artefact_suspected <- slope >= shuffle_ci[1] && slope <= shuffle_ci[2]

  # --- artefact guard 2: mixture-vs-single-line (are there really 2 subgroups?)
  mix <- .prr_mixture_diagnostic(pot, dlt)

  # follow-up regressed on initial: a genuine rule gives slope ~ (1 - prop)
  fu_init_slope <- unname(stats::coef(stats::lm(followup ~ initial))[2])

  res <- list(
    slope = slope, slope_ci = unname(ci), expected_prop = expected_prop,
    p_vs_expected = p_vs_expected,
    followup_on_initial_slope = fu_init_slope,
    n = n, n_ceiling = sum(ceiling), n_fitter = sum(is_fitter),
    classification = classification, observed_prop = obs_prop,
    mixture_preferred = mix$mixture_preferred, mixture_bic = mix$bic,
    shuffle_slope_ci = shuffle_ci, artefact_suspected = artefact_suspected)

  PhysioCore::AnalysisResult(
    type = "proportional_recovery", result = res,
    parameters = list(max_score = max_score, expected_prop = expected_prop,
                      ceiling_frac = ceiling_frac),
    estimate = slope, method = "proportional_recovery_rule")
}

# Two-component vs single-line diagnostic: fit a 1- and a 2-component regression
# mixture of delta ~ 0 + potential (flexmix) and compare BIC. A properly
# penalised mixture BIC avoids the false positive of splitting a single noisy
# line, so this flags only a genuine fitter/non-fitter subgroup structure.
.prr_mixture_diagnostic <- function(pot, dlt) {
  if (!requireNamespace("flexmix", quietly = TRUE)) {
    return(list(bic = c(single = NA_real_, mixture = NA_real_),
                mixture_preferred = NA))
  }
  dat <- data.frame(dlt = dlt, pot = pot)
  m1 <- tryCatch(flexmix::flexmix(dlt ~ 0 + pot, data = dat, k = 1),
                 error = function(e) NULL)
  m2 <- tryCatch(flexmix::flexmix(dlt ~ 0 + pot, data = dat, k = 2,
                                  control = list(minprior = 0.1)),
                 error = function(e) NULL)
  if (is.null(m1)) {
    return(list(bic = c(single = NA_real_, mixture = NA_real_),
                mixture_preferred = NA))
  }
  bic1 <- stats::BIC(m1)
  # a converged 2-component model that did not collapse to a single cluster
  bic2 <- if (!is.null(m2) && length(unique(flexmix::clusters(m2))) == 2L) {
    stats::BIC(m2)
  } else Inf
  list(bic = c(single = bic1, mixture = bic2),
       mixture_preferred = is.finite(bic2) && bic2 < bic1 - 2)
}
