# Rubin's rules for pooling estimates across multiply-imputed data sets, with
# the Barnard-Rubin (1999) small-sample degrees-of-freedom adjustment - matching
# mice::pool term by term.

# per-fit (term, estimate, variance, dfcom) from any model with coef()/vcov()
.pool_extract <- function(fit) {
  est <- stats::coef(fit)
  v <- tryCatch(diag(as.matrix(stats::vcov(fit))), error = function(e) NULL)
  if (is.null(v) || length(v) != length(est)) {
    stop("cannot extract a variance for every coefficient of a fit.",
         call. = FALSE)
  }
  dfcom <- tryCatch(stats::df.residual(fit), error = function(e) NULL)
  if (is.null(dfcom) || !is.finite(dfcom)) dfcom <- Inf
  data.frame(term = names(est), estimate = as.numeric(est),
             variance = as.numeric(v), dfcom = dfcom,
             row.names = NULL, stringsAsFactors = FALSE)
}

# Barnard & Rubin (1999) pooled degrees of freedom, matching mice 3.x's
# single-expression form (no lambda floor; finite at lambda = 0).
.barnard_rubin <- function(m, b, t, dfcom) {
  lambda <- (1 + 1 / m) * b / t
  dfold <- (m - 1) / lambda^2
  if (is.infinite(dfcom)) return(dfold)
  tmp <- (1 - lambda) * (1 + dfcom) * dfcom
  (m - 1) * tmp / ((dfcom + 3) * (m - 1) + lambda^2 * tmp)
}

#' Pool estimates across imputations (Rubin's rules)
#'
#' Combines the coefficient estimates from a list of models fitted to multiply
#' imputed data sets using Rubin's rules with the Barnard-Rubin degrees-of-
#' freedom adjustment. The pooled estimate, total variance, df, confidence
#' interval, and fraction of missing information reproduce \code{mice::pool}.
#'
#' @param fits A list of fitted models (each with \code{coef()} and
#'   \code{vcov()}), or a list of tidy data frames with columns \code{term},
#'   \code{estimate}, and \code{std.error} (and optionally \code{df.residual}).
#' @param conf_level Confidence level for the pooled interval (default 0.95).
#' @return An \code{AnalysisResult} (type \code{"pooled_estimates"}) whose
#'   \code{result$estimates} is a data frame of pooled \code{estimate},
#'   \code{std.error}, \code{df}, \code{statistic}, \code{p.value},
#'   \code{conf.low}, \code{conf.high}, and \code{fmi} per term.
#' @references Rubin 1987; Barnard & Rubin 1999. \code{mice::pool}.
#' @seealso [multipleImputation()], [analyseEstimand()]
#' @export
#' @examples
#' set.seed(1)
#' fits <- lapply(1:5, function(i) lm(mpg ~ hp + wt,
#'   data = mtcars[sample(nrow(mtcars), replace = TRUE), ]))
#' poolEstimates(fits)
poolEstimates <- function(fits, conf_level = 0.95) {
  if (!is.list(fits) || !length(fits)) {
    stop("'fits' must be a non-empty list of fitted models.", call. = FALSE)
  }
  m <- length(fits)
  tidy <- lapply(fits, function(f) {
    if (is.data.frame(f)) {
      miss <- setdiff(c("term", "estimate", "std.error"), names(f))
      if (length(miss)) {
        stop("tidy data frame input is missing required column(s): ",
             paste(miss, collapse = ", "), ".", call. = FALSE)
      }
      f$variance <- f$std.error^2
      f$dfcom <- if (!is.null(f$df.residual)) f$df.residual else Inf
      f[, c("term", "estimate", "variance", "dfcom")]
    } else .pool_extract(f)
  })
  terms <- tidy[[1]]$term
  if (!all(vapply(tidy, function(t) identical(t$term, terms), logical(1)))) {
    stop("all fits must share the same coefficient terms.", call. = FALSE)
  }

  out <- lapply(seq_along(terms), function(j) {
    qhat <- vapply(tidy, function(t) t$estimate[j], numeric(1))
    u <- vapply(tidy, function(t) t$variance[j], numeric(1))
    dfcom <- tidy[[1]]$dfcom[j]
    qbar <- mean(qhat); ubar <- mean(u)
    b <- if (m > 1L) stats::var(qhat) else 0
    tvar <- ubar + (1 + 1 / m) * b
    # an aliased/rank-deficient term (NA coefficient in some fit) gives b = NA;
    # keep it out of the logical guard so it yields an NA row (as mice::pool does)
    df <- if (m == 1L) dfcom
      else if (is.finite(b)) .barnard_rubin(m, b, tvar, dfcom)
      else NA_real_
    riv <- (1 + 1 / m) * b / ubar
    fmi <- (riv + 2 / (df + 3)) / (riv + 1)
    se <- sqrt(tvar)
    tcrit <- stats::qt(1 - (1 - conf_level) / 2, df)
    data.frame(
      term = terms[j], estimate = qbar, std.error = se, df = df,
      statistic = qbar / se, p.value = 2 * stats::pt(-abs(qbar / se), df),
      conf.low = qbar - tcrit * se, conf.high = qbar + tcrit * se,
      fmi = fmi, row.names = NULL, stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, out)

  # wrap the table in a list: the AnalysisResult slot coerces `result` with
  # as.list(), which would otherwise flatten a data frame into its columns.
  PhysioCore::AnalysisResult(
    type = "pooled_estimates", result = list(estimates = res),
    parameters = list(m = m, conf_level = conf_level),
    estimate = stats::setNames(res$estimate, res$term), method = "rubin_pool")
}

#' Analyse an estimand under multiple imputation for dropout
#'
#' End-to-end analysis of a declared \code{\link{defineEstimand}} estimand:
#' multiply-impute the dropout-missing longitudinal response (in wide form, so
#' the within-subject correlation is respected), fit an MMRM
#' (\code{\link{fitMMRM}}) to each completed data set, and pool the fixed
#' effects with Rubin's rules ([poolEstimates()]). Treatment-policy and
#' hypothetical strategies use the standard MAR imputation; other strategies
#' warn that they need a bespoke imputation model.
#'
#' @param data Long-format data (subject x time rows) with the response,
#'   treatment, time, subject, and any (subject-level) covariates.
#' @param estimand An \code{"estimand"} from [defineEstimand()].
#' @param response,treatment,time,subject Column names.
#' @param covariates Optional subject-level covariate column names.
#' @param m Number of imputations (default 20).
#' @param method mice method (default \code{"norm"}, appropriate for a
#'   continuous longitudinal response).
#' @param seed Optional integer seed.
#' @param covariance MMRM covariance structure (default \code{"unstructured"}).
#' @return An \code{AnalysisResult} (type \code{"estimand_analysis"}) carrying
#'   the estimand attributes and the pooled fixed-effect table.
#' @references ICH E9(R1); Rubin 1987; Mallinckrodt 2008 (MMRM).
#' @seealso [defineEstimand()], [multipleImputation()], [poolEstimates()]
#' @export
analyseEstimand <- function(data, estimand, response, treatment, time, subject,
                            covariates = NULL, m = 20, method = "norm",
                            seed = NULL, covariance = "unstructured") {
  if (!inherits(estimand, "estimand")) {
    stop("'estimand' must come from defineEstimand().", call. = FALSE)
  }
  if (!requireNamespace("mice", quietly = TRUE)) {
    stop("analyseEstimand() requires the 'mice' package.", call. = FALSE)
  }
  if (!estimand$strategy %in% c("treatment-policy", "hypothetical")) {
    warning("strategy '", estimand$strategy, "' needs a bespoke imputation ",
            "model; using the standard MAR imputation as an approximation.",
            call. = FALSE)
  }
  d <- data.frame(.subj = data[[subject]], .time = data[[time]],
                  .y = as.numeric(data[[response]]),
                  .trt = data[[treatment]])
  for (cv in covariates) d[[cv]] <- data[[cv]]

  # wide: one row per subject, one response column per time (dropout -> NA)
  wide <- stats::reshape(d, idvar = ".subj", timevar = ".time",
                         v.names = ".y", direction = "wide")
  ycols <- grep("^\\.y\\.", names(wide), value = TRUE)

  # Impute separately within each treatment arm (Carpenter & Kenward): a pooled
  # additive imputation model would shrink the treatment-by-time difference
  # (uncongenial with the MMRM analysis), biasing the effect toward the null.
  arms <- split(wide, wide$.trt)
  imp_cols <- setdiff(names(wide), c(".subj", ".trt"))
  imps <- lapply(arms, function(w)
    multipleImputation(w[, imp_cols, drop = FALSE], method = method, m = m,
                       seed = seed))

  fits <- lapply(seq_len(m), function(i) {
    comp <- do.call(rbind, lapply(names(arms), function(a) {
      ci <- mice::complete(imps[[a]], i)
      ci$.subj <- arms[[a]]$.subj; ci$.trt <- a; ci
    }))
    long <- stats::reshape(comp, idvar = ".subj", varying = ycols,
                           v.names = ".y", timevar = ".time",
                           times = sub("^\\.y\\.", "", ycols), direction = "long")
    long$.time <- factor(long$.time)
    fitMMRM(long, ".y", ".trt", ".time", ".subj",
            covariates = covariates, covariance = covariance)@result$fit
  })
  pooled <- poolEstimates(fits)

  PhysioCore::AnalysisResult(
    type = "estimand_analysis",
    result = list(estimand = estimand, pooled = pooled@result$estimates,
                  strategy = estimand$strategy, m = m),
    parameters = list(m = m, method = method, covariance = covariance),
    estimand = estimand[c("treatment", "population", "endpoint", "strategy",
                          "summary_measure")],
    method = "mi_mmrm_rubin")
}
