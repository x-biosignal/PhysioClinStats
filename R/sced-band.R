# SCED visual-analysis aids: the two-standard-deviation band (Nelson) and the
# celeration line (White & Haring split-middle, plus an OLS trend).

#' Two-standard-deviation band decision rule
#'
#' Draws the Nelson (1984) 2-SD band from the baseline (A) phase - the baseline
#' mean plus and minus two baseline standard deviations - and flags the
#' intervention (B) phase for a systematic shift: the rule fires when a run of
#' `consecutive` or more successive B points falls on the same side beyond the
#' band.
#'
#' @param A_data,B_data Numeric baseline and intervention observations.
#' @param improvement `"increase"` (default) or `"decrease"`; sets which side of
#'   the band counts as improvement.
#' @param k Band half-width in baseline SDs (default 2).
#' @param consecutive Number of successive out-of-band points that trigger the
#'   flag (default 2, Nelson's rule).
#' @return An `AnalysisResult` (`type = "sced_2sd"`) whose `estimate` is the
#'   logical decision, with `result$mean`, `result$sd`, `result$upper`,
#'   `result$lower`, `result$outside` (per-B-point side: -1/0/+1) and
#'   `result$first_run_at` (index of the first triggering point, or `NA`).
#' @references Nelson LS (1984). The Shewhart control chart - tests for special
#'   causes. Journal of Quality Technology, 16(4), 237-239. Applied to SCED in
#'   Gast & Ledford (2014).
#' @seealso [scedCelerationLine()], [scedABAB()]
#' @export
#' @examples
#' scedTwoSDBand(c(10, 12, 11, 9, 10), c(15, 16, 17, 16, 18))
scedTwoSDBand <- function(A_data, B_data = NULL,
                          improvement = c("increase", "decrease"),
                          k = 2, consecutive = 2L) {
  improvement <- match.arg(improvement)
  stopifnot(is.numeric(k), length(k) == 1, k > 0,
            consecutive >= 1L)
  ph <- .sced_phases(A_data, B_data)
  A <- ph$A; B <- ph$B
  mA <- mean(A); sA <- stats::sd(A)
  upper <- mA + k * sA; lower <- mA - k * sA

  side <- ifelse(B > upper, 1L, ifelse(B < lower, -1L, 0L))
  good <- if (improvement == "increase") 1L else -1L

  # first run of `consecutive` successive points on the improvement side
  run_at <- NA_integer_
  if (length(side) >= consecutive) {
    hit <- side == good
    r <- 0L
    for (i in seq_along(hit)) {
      r <- if (isTRUE(hit[i])) r + 1L else 0L
      if (r >= consecutive) { run_at <- i - consecutive + 1L; break }
    }
  }
  decision <- !is.na(run_at)

  .sced_result("sced_2sd", decision, "2-SD band",
               extra = list(mean = mA, sd = sA, upper = upper, lower = lower,
                            outside = side, first_run_at = run_at,
                            k = k, consecutive = as.integer(consecutive),
                            improvement = improvement),
               params = list(k = k, consecutive = as.integer(consecutive),
                             improvement = improvement))
}

# Split-middle line for one phase (White & Haring 1980):
#   1. order by x; split into first/second half (drop the middle point if odd)
#   2. in each half take the intersection of the median x and median y
#   3. the line through those two points gives the slope (celeration)
#   4. slide it vertically so the median residual is zero (equal points each side)
.split_middle <- function(x, y) {
  o <- order(x); x <- x[o]; y <- y[o]
  n <- length(x)
  if (n < 2L) stop("Need at least 2 points for a celeration line.", call. = FALSE)
  half <- n %/% 2L
  first <- seq_len(half)
  second <- (n - half + 1L):n
  mx1 <- stats::median(x[first]);  my1 <- stats::median(y[first])
  mx2 <- stats::median(x[second]); my2 <- stats::median(y[second])
  slope <- (my2 - my1) / (mx2 - mx1)
  intercept <- my1 - slope * mx1
  # split-middle adjustment: zero the median residual
  intercept <- intercept + stats::median(y - (slope * x + intercept))
  list(slope = slope, intercept = intercept)
}

#' Celeration (trend) line for a single-case phase
#'
#' Fits a trend line to one phase, either by the White & Haring split-middle
#' method (the standard SCED hand method) or by ordinary least squares.
#'
#' @param value Numeric phase observations, or pass `time`/`value` explicitly.
#' @param time Optional numeric session index (default `seq_along(value)`).
#' @param method `"split_middle"` (default) or `"ols"`.
#' @return An `AnalysisResult` (`type = "sced_celeration"`) whose `estimate` is
#'   the slope (celeration per session), with `result$intercept`,
#'   `result$fitted`, `result$bounce` (the ratio of the largest positive to the
#'   largest negative residual magnitude) and `result$celeration_ratio` (the
#'   multiplicative change over `ratio_period` sessions).
#' @param ratio_period Sessions spanned by the reported celeration ratio
#'   (default 7, a weekly ratio for daily data).
#' @references White OR, Haring NG (1980). Exceptional Teaching, 2nd ed.
#'   Columbus, OH: Merrill.
#' @seealso [scedTwoSDBand()], [scedABAB()]
#' @export
#' @examples
#' scedCelerationLine(c(4, 6, 5, 8, 7, 10, 9))
scedCelerationLine <- function(value, time = NULL,
                               method = c("split_middle", "ols"),
                               ratio_period = 7) {
  method <- match.arg(method)
  y <- as.numeric(value)
  x <- if (is.null(time)) seq_along(y) else as.numeric(time)
  stopifnot(length(x) == length(y))
  keep <- !is.na(x) & !is.na(y); x <- x[keep]; y <- y[keep]
  if (length(y) < 2L) stop("Need at least 2 observations.", call. = FALSE)

  fit <- if (method == "split_middle") {
    .split_middle(x, y)
  } else {
    cf <- stats::coef(stats::lm(y ~ x))
    list(slope = unname(cf[2]), intercept = unname(cf[1]))
  }
  fitted <- fit$slope * x + fit$intercept
  resid <- y - fitted
  pos <- max(resid[resid > 0], 0); neg <- max(-resid[resid < 0], 0)
  bounce <- if (neg > 0) pos / neg else NA_real_
  # celeration ratio: multiplicative change over ratio_period sessions, using
  # the fitted values at the ends of one period (guards against zero/negative).
  y0 <- fit$intercept + fit$slope * x[1]
  yP <- y0 + fit$slope * ratio_period
  cel_ratio <- if (is.finite(y0) && y0 != 0) yP / y0 else NA_real_

  .sced_result("sced_celeration", fit$slope, paste0("celeration (", method, ")"),
               extra = list(intercept = fit$intercept, fitted = fitted,
                            residuals = resid, bounce = bounce,
                            celeration_ratio = cel_ratio, method = method,
                            ratio_period = ratio_period),
               params = list(method = method, ratio_period = ratio_period))
}
