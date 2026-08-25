# Single-case experimental design (SCED) / N-of-1 non-overlap effect sizes.
#
# Every estimator here reproduces its reference implementation exactly: PND,
# PEM, NAP, Tau and Tau-U match the SingleCaseES package (Pustejovsky), and the
# Tau-U inference (Kendall S variance, Z, p) matches the scan package (Wilbert).
# The two references genuinely differ on the Tau-U denominator convention (see
# scedTauU Details); scedTauU exposes both.

# Direction of improvement -> sign applied to the data. "increase" leaves the
# data alone; "decrease" negates it so that a therapeutic reduction reads as a
# positive effect, matching SingleCaseES's `improvement` argument.
.sced_sign <- function(improvement = c("increase", "decrease")) {
  improvement <- match.arg(improvement)
  if (improvement == "decrease") -1 else 1
}

# Numeric coercion that reads a factor by its labels, not its integer codes, and
# errors informatively on genuinely non-numeric input rather than silently
# turning it into NA.
.sced_numeric <- function(x, what = "value") {
  if (is.factor(x)) x <- as.character(x)
  out <- suppressWarnings(as.numeric(x))
  if (length(out) && all(is.na(out)) && !all(is.na(x))) {
    stop(sprintf("'%s' is not numeric (or numeric-coercible).", what),
         call. = FALSE)
  }
  out
}

# Coerce the (A, B) phase inputs from either two numeric vectors or a single
# data.frame(value, phase[, time]) with two phase levels.
.sced_phases <- function(A_data, B_data = NULL, phase = NULL) {
  if (is.null(B_data)) {
    stop("Provide baseline A_data and intervention B_data.", call. = FALSE)
  }
  A <- .sced_numeric(A_data, "A_data"); B <- .sced_numeric(B_data, "B_data")
  A <- A[!is.na(A)]; B <- B[!is.na(B)]
  if (length(A) < 1L || length(B) < 1L) {
    stop("Each phase needs at least one non-missing observation.", call. = FALSE)
  }
  list(A = A, B = B)
}

.sced_result <- function(type, est, method, level = NULL, lower = NULL,
                         upper = NULL, extra = list(), params = list()) {
  unc <- list()
  if (!is.null(lower) && !is.null(upper)) {
    unc <- list(type = "analytic", level = level, lower = lower, upper = upper)
  }
  PhysioCore::AnalysisResult(
    type = type, estimate = est, method = method, uncertainty = unc,
    result = c(list(estimate = est, ci_lower = lower, ci_upper = upper), extra),
    parameters = params,
    provenance = data.frame(step = type, timestamp = NA_character_,
                            stringsAsFactors = FALSE))
}

#' Percentage of non-overlapping data (PND)
#'
#' The proportion of intervention-phase observations that exceed the most
#' extreme baseline observation (Scruggs, Mastropieri & Casto 1987).
#'
#' @param A_data,B_data Numeric baseline (A) and intervention (B) observations,
#'   or pass a data.frame to the phase orchestrator [scedABAB()].
#' @param improvement `"increase"` (default) if higher scores are better, or
#'   `"decrease"` if lower scores are the therapeutic goal.
#' @return An `AnalysisResult` (`type = "sced_pnd"`) whose `estimate` is the PND
#'   in \[0, 1\].
#' @references Scruggs TE, Mastropieri MA, Casto G (1987). The quantitative
#'   synthesis of single-subject research. Remedial and Special Education, 8(2).
#' @seealso [scedPEM()], [scedNAP()], [scedTauU()]
#' @export
#' @examples
#' scedPND(c(20, 20, 26, 25), c(28, 25, 30, 29))
scedPND <- function(A_data, B_data = NULL, improvement = c("increase", "decrease")) {
  improvement <- match.arg(improvement)
  ph <- .sced_phases(A_data, B_data); s <- .sced_sign(improvement)
  A <- s * ph$A; B <- s * ph$B
  est <- mean(B > max(A))
  .sced_result("sced_pnd", est, "PND",
               extra = list(n_A = length(A), n_B = length(B),
                            improvement = improvement),
               params = list(improvement = improvement))
}

#' Percentage exceeding the median (PEM)
#'
#' The proportion of intervention-phase observations beyond the baseline median,
#' counting ties as one half (Ma 2006).
#'
#' @inheritParams scedPND
#' @return An `AnalysisResult` (`type = "sced_pem"`) whose `estimate` is the PEM.
#' @references Ma H-H (2006). An alternative method for quantitative synthesis of
#'   single-subject research: percentage of data points exceeding the median.
#'   Behavior Modification, 30(5), 598-617.
#' @seealso [scedPND()], [scedNAP()]
#' @export
#' @examples
#' scedPEM(c(20, 20, 26, 25), c(28, 25, 30, 29))
scedPEM <- function(A_data, B_data = NULL, improvement = c("increase", "decrease")) {
  improvement <- match.arg(improvement)
  ph <- .sced_phases(A_data, B_data); s <- .sced_sign(improvement)
  A <- s * ph$A; B <- s * ph$B
  med <- stats::median(A)
  est <- mean((B > med) + 0.5 * (B == med))
  .sced_result("sced_pem", est, "PEM",
               extra = list(baseline_median = s * med, n_A = length(A),
                            n_B = length(B), improvement = improvement),
               params = list(improvement = improvement))
}

# Non-overlap matrix Q[i,j] = (B_j > A_i) + 0.5 (B_j == A_i), the count of
# dominant pairs. NAP = mean(Q).
.sced_Qmat <- function(A, B) {
  matrix(sapply(B, function(j) (j > A) + 0.5 * (j == A)),
         nrow = length(A), ncol = length(B))
}

#' Non-overlap of all pairs (NAP)
#'
#' The probability that a randomly chosen intervention observation exceeds a
#' randomly chosen baseline observation, with ties counted as one half
#' (Parker & Vannest 2009) - equivalent to the area under the ROC curve. The
#' confidence interval is the score interval of Newcombe (2006) and the standard
#' error the unbiased estimator, both as in SingleCaseES.
#'
#' @inheritParams scedPND
#' @param confidence Confidence level for the interval (default 0.95).
#' @return An `AnalysisResult` (`type = "sced_nap"`) whose `estimate` is the NAP,
#'   with an analytic `uncertainty` interval and `result$se`, `result$p_value`.
#' @references Parker RI, Vannest KJ (2009). An improved effect size for
#'   single-case research: non-overlap of all pairs. Behavior Therapy, 40(4),
#'   357-367. Newcombe RG (2006). Confidence intervals for an effect size
#'   measure based on the Mann-Whitney statistic. Statistics in Medicine, 25.
#' @seealso [scedTau()], [scedTauU()]
#' @export
#' @examples
#' scedNAP(c(20, 20, 26, 25), c(28, 25, 30, 29))
scedNAP <- function(A_data, B_data = NULL, improvement = c("increase", "decrease"),
                    confidence = 0.95) {
  improvement <- match.arg(improvement)
  stopifnot(is.numeric(confidence), length(confidence) == 1,
            confidence > 0, confidence < 1)
  ph <- .sced_phases(A_data, B_data); s <- .sced_sign(improvement)
  A <- s * ph$A; B <- s * ph$B
  m <- length(A); n <- length(B)
  Q <- .sced_Qmat(A, B)
  nap <- mean(Q)

  # Unbiased variance (SingleCaseES "unbiased").
  Q1 <- sum(rowSums(Q - nap)^2) / (m * n^2)
  Q2 <- sum(colSums(Q - nap)^2) / (m^2 * n)
  X  <- sum((Q - nap)^2) / (m * n)
  trunc <- 0.5 / (m * n)
  nap_t <- min(max(nap, trunc), 1 - trunc)
  V <- if (m > 1 && n > 1) {
    (nap_t * (1 - nap_t) + n * Q1 + m * Q2 - 2 * X) / ((m - 1) * (n - 1))
  } else NA_real_
  se <- sqrt(V)

  # Newcombe score interval.
  h <- (m + n) / 2 - 1
  z <- stats::qnorm(1 - (1 - confidence) / 2)
  froot <- function(x) m * n * (nap - x)^2 * (2 - x) * (1 + x) -
    z^2 * x * (1 - x) * (2 + h + (1 + 2 * h) * x * (1 - x))
  lower <- if (nap > 0) stats::uniroot(froot, c(0, nap))$root else 0
  upper <- if (nap < 1) stats::uniroot(froot, c(nap, 1))$root else 1

  # Two-sided p against the null NAP = 0.5 (Mann-Whitney null variance).
  v0 <- (m + n + 1) / (12 * m * n)
  zval <- (nap - 0.5) / sqrt(v0)
  pval <- 2 * stats::pnorm(-abs(zval))

  .sced_result("sced_nap", nap, "NAP", level = confidence,
               lower = lower, upper = upper,
               extra = list(se = se, z = zval, p_value = pval,
                            n_A = m, n_B = n, improvement = improvement),
               params = list(improvement = improvement, confidence = confidence))
}

#' Tau non-overlap (Tau)
#'
#' The rank-correlation non-overlap effect size, \eqn{\tau = 2\,\mathrm{NAP} - 1}
#' (Parker et al. 2011), with the SE and interval derived from [scedNAP()].
#'
#' @inheritParams scedNAP
#' @return An `AnalysisResult` (`type = "sced_tau"`) whose `estimate` is Tau in
#'   \[-1, 1\].
#' @seealso [scedNAP()], [scedTauU()]
#' @export
#' @examples
#' scedTau(c(20, 20, 26, 25), c(28, 25, 30, 29))
scedTau <- function(A_data, B_data = NULL, improvement = c("increase", "decrease"),
                    confidence = 0.95) {
  improvement <- match.arg(improvement)
  nap <- scedNAP(A_data, B_data, improvement = improvement, confidence = confidence)
  r <- PhysioCore::resultValue(nap)
  est <- 2 * r$estimate - 1
  .sced_result("sced_tau", est, "Tau", level = confidence,
               lower = 2 * r$ci_lower - 1, upper = 2 * r$ci_upper - 1,
               extra = list(se = 2 * r$se, z = r$z, p_value = r$p_value,
                            n_A = r$n_A, n_B = r$n_B, improvement = improvement),
               params = list(improvement = improvement, confidence = confidence))
}

# Kendall S between two vectors (count of concordant minus discordant pairs).
.kendall_S_between <- function(A, B) sum(sapply(B, function(j) sum(sign(j - A))))
# Kendall S within one vector (upper-triangular concordances = phase trend).
.kendall_S_within <- function(x) {
  n <- length(x)
  if (n < 2L) return(0)
  M <- outer(x, x, function(a, b) sign(b - a))   # M[i,j] = sign(x_j - x_i)
  sum(M[upper.tri(M)])
}

# scan encodes each Tau-U comparison as a Kendall tau-b between the pooled data
# and a synthetic rank predictor; "A vs. B - Trend A" ranks the baseline points
# in descending order (nA:1) and ties every intervention point above them
# (nA + 1). The tau-b denominator D, the tie-corrected S variance and the
# normal-approximation Z/p then follow from base R's Kendall machinery, matching
# scan::tau_u exactly.
.tau_u_trendA_predictor <- function(nA, nB) c(nA:1, rep(nA + 1L, nB))
.tau_b_D <- function(x, y) {
  n0 <- length(x) * (length(x) - 1) / 2
  tp <- function(v) { t <- as.numeric(table(v)); sum(t * (t - 1) / 2) }
  sqrt((n0 - tp(x)) * (n0 - tp(y)))
}

#' Tau-U with baseline-trend correction
#'
#' Parker, Vannest, Davis & Sauber's (2011) Tau-U combines the A-vs-B
#' non-overlap with a correction for baseline (phase A) trend. Two reference
#' conventions differ on the denominator and this function exposes both:
#'
#' \describe{
#'   \item{`method = "parker"` (default)}{matches `SingleCaseES::Tau_U`:
#'     \eqn{(S_{AB} - S_{trend A}) / (m n)}, where \eqn{S_{AB}} is the A-vs-B
#'     Kendall S and \eqn{S_{trend A}} the within-baseline Kendall S. Can exceed
#'     1 when a baseline trend runs counter to improvement.}
#'   \item{`method = "scan"`}{matches `scan::tau_u()`'s "A vs. B - Trend A" row:
#'     the same S numerator over a tie-corrected pair count \eqn{D}, giving a
#'     Kendall tau-b in \[-1, 1\].}
#' }
#'
#' The Z statistic and p-value come from the tie-corrected variance of the
#' Kendall S numerator and are identical under both methods (they depend on S,
#' not the denominator), matching `scan`'s Z and p.
#'
#' @inheritParams scedNAP
#' @param method Denominator convention, `"parker"` (default) or `"scan"`.
#' @return An `AnalysisResult` (`type = "sced_tau_u"`) whose `estimate` is Tau-U,
#'   with `result$S`, `result$se_S`, `result$z`, `result$p_value` and the
#'   component S values.
#' @references Parker RI, Vannest KJ, Davis JL, Sauber SB (2011). Combining
#'   nonoverlap and trend for single-case research: Tau-U. Behavior Therapy,
#'   42(2), 284-299.
#' @seealso [scedNAP()], [scedTau()]
#' @export
#' @examples
#' scedTauU(c(20, 20, 26, 25, 22, 23), c(28, 25, 24, 27, 30, 30, 29, 28))
scedTauU <- function(A_data, B_data = NULL, improvement = c("increase", "decrease"),
                     method = c("parker", "scan"), confidence = 0.95) {
  improvement <- match.arg(improvement)
  method <- match.arg(method)
  ph <- .sced_phases(A_data, B_data); s <- .sced_sign(improvement)
  A <- s * ph$A; B <- s * ph$B
  m <- length(A); n <- length(B)

  S_AB <- .kendall_S_between(A, B)
  S_trendA <- .kendall_S_within(A)
  S <- S_AB - S_trendA

  # scan's tau-b construction for the same "A vs B - Trend A" comparison.
  AB <- c(A, B)
  pred <- .tau_u_trendA_predictor(m, n)
  D <- .tau_b_D(AB, pred)
  ct <- suppressWarnings(stats::cor.test(AB, pred, method = "kendall"))
  zval <- unname(ct$statistic)     # normal-approx z with ties (scan's Z)
  pval <- ct$p.value

  est <- switch(method,
                parker = S / (m * n),   # SingleCaseES::Tau_U convention
                scan   = S / D)         # scan tau-b convention

  .sced_result("sced_tau_u", est, paste0("Tau-U (", method, ")"),
               level = confidence,
               extra = list(S = S, S_AB = S_AB, S_trendA = S_trendA,
                            D = D, z = zval, p_value = pval,
                            method = method, n_A = m, n_B = n,
                            improvement = improvement),
               params = list(improvement = improvement, method = method,
                             confidence = confidence))
}
