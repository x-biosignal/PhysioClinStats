# Mixed models for repeated measures (MMRM): the longitudinal-RCT standard with
# an unstructured within-subject covariance and Kenward-Roger / Satterthwaite
# denominator df, via the mmrm package, with an nlme::gls fallback.

.mmrm_cov_fun <- c(unstructured = "us", `ar1` = "ar1", `compound-symmetry` = "cs",
                   toeplitz = "toep")

#' Fit a mixed model for repeated measures (MMRM)
#'
#' Fits the Mallinckrodt MMRM for a longitudinal randomised trial: fixed effects
#' for treatment, categorical time and their interaction (plus any covariates),
#' with a within-subject covariance over the repeated time points. Uses the
#' `mmrm` package when available - with Kenward-Roger or Satterthwaite adjusted
#' degrees of freedom - and falls back to `nlme::gls` otherwise. The fallback
#' reproduces the `mmrm` fixed effects for the unstructured covariance (its
#' primary use) and approximates the homogeneous `cs`/`ar1`/`toeplitz`
#' structures; it reports between-within (containment) degrees of freedom and
#' emits a message, since Kenward-Roger / Satterthwaite df need `mmrm`.
#'
#' @param data A long-format data.frame with one row per subject-time.
#' @param response,treatment,time,subject Column names of the outcome, the
#'   treatment factor, the categorical time factor and the subject id.
#' @param covariates Optional character vector of extra fixed-effect covariate
#'   columns (e.g. baseline value, stratifiers).
#' @param covariance Within-subject covariance: `"unstructured"` (default),
#'   `"ar1"`, `"compound-symmetry"` or `"toeplitz"`.
#' @param df Denominator-df method: `"kenward-roger"` (default) or
#'   `"satterthwaite"` (used by the `mmrm` path only).
#' @return An `AnalysisResult` (`type = "mmrm"`) whose `estimate` is the
#'   fixed-effect coefficient vector, with `result$coefficients` (estimate, SE,
#'   df, t, p), `result$backend` (`"mmrm"` or `"gls"`), `result$fit` and
#'   `result$formula`.
#' @references Mallinckrodt CH et al. (2008). Recommendations for the primary
#'   analysis of continuous endpoints in longitudinal clinical trials. Drug
#'   Information Journal, 42. Sabanes Bove D et al. (2023). mmrm: Mixed Models
#'   for Repeated Measures. R package.
#' @seealso [fitMixedModel()], [estimatedMarginalMeans()]
#' @export
#' @examples
#' if (requireNamespace("mmrm", quietly = TRUE)) {
#'   data(fev_data, package = "mmrm")
#'   fitMMRM(fev_data, "FEV1", "ARMCD", "AVISIT", "USUBJID",
#'           covariates = c("RACE", "SEX"))
#' }
fitMMRM <- function(data, response, treatment, time, subject,
                    covariates = NULL,
                    covariance = c("unstructured", "ar1",
                                   "compound-symmetry", "toeplitz"),
                    df = c("kenward-roger", "satterthwaite")) {
  covariance <- match.arg(covariance)
  df <- match.arg(df)
  stopifnot(is.data.frame(data))
  cols <- c(response, treatment, time, subject, covariates)
  miss <- setdiff(cols, names(data))
  if (length(miss)) {
    stop(sprintf("Columns not in 'data': %s", paste(miss, collapse = ", ")),
         call. = FALSE)
  }
  # time and subject must be factors for a repeated-measures structure
  data[[time]] <- as.factor(data[[time]])
  data[[subject]] <- as.factor(data[[subject]])

  fixed_terms <- c(covariates, sprintf("%s*%s", treatment, time))

  if (requireNamespace("mmrm", quietly = TRUE)) {
    cov_fun <- .mmrm_cov_fun[[covariance]]
    cov_term <- sprintf("%s(%s | %s)", cov_fun, time, subject)
    form <- stats::reformulate(c(fixed_terms, cov_term), response = response)
    meth <- if (df == "kenward-roger") "Kenward-Roger" else "Satterthwaite"
    fit <- mmrm::mmrm(form, data = data, method = meth)
    co <- summary(fit)$coefficients
    coefs <- data.frame(term = rownames(co), estimate = co[, "Estimate"],
                        std_error = co[, "Std. Error"], df = co[, "df"],
                        statistic = co[, "t value"], p_value = co[, "Pr(>|t|)"],
                        row.names = NULL, stringsAsFactors = FALSE)
    backend <- "mmrm"
  } else {
    if (!requireNamespace("nlme", quietly = TRUE)) {
      stop("fitMMRM needs the mmrm package (or nlme for the fallback).",
           call. = FALSE)
    }
    form <- stats::reformulate(fixed_terms, response = response)
    cor_struct <- switch(covariance,
      unstructured = nlme::corSymm(form = stats::as.formula(
        sprintf("~ as.integer(%s) | %s", time, subject))),
      ar1 = nlme::corAR1(form = stats::as.formula(
        sprintf("~ as.integer(%s) | %s", time, subject))),
      `compound-symmetry` = nlme::corCompSymm(form = stats::as.formula(
        sprintf("~ as.integer(%s) | %s", time, subject))),
      toeplitz = nlme::corARMA(form = stats::as.formula(
        sprintf("~ as.integer(%s) | %s", time, subject)),
        p = length(levels(data[[time]])) - 1L, q = 0L))
    # Only the unstructured covariance has per-time variances; mmrm's cs / ar1 /
    # toep are single-variance (homogeneous) structures, so applying varIdent to
    # them would silently fit the heterogeneous variant (csh / ar1h / toeph).
    weights <- if (covariance == "unstructured") {
      nlme::varIdent(form = stats::as.formula(sprintf("~ 1 | %s", time)))
    } else NULL
    fit <- nlme::gls(form, data = data, correlation = cor_struct,
                     weights = weights, method = "REML")
    tt <- summary(fit)$tTable
    # gls has no adjusted df, so use between-within (containment) denominator df:
    # a term constant within every subject (e.g. the treatment main effect) is a
    # between-subject effect tested on ~ n_subject df, the rest on the residual
    # within-subject df. This keeps the between-subject inference from being
    # anti-conservative (mmrm's Kenward-Roger df are unavailable here).
    mm <- stats::model.matrix(form, data)
    subj <- as.factor(data[[subject]])
    between <- apply(mm, 2, function(col)
      all(tapply(col, subj, function(v) diff(range(v)) == 0)))
    between <- between[rownames(tt)]
    n_subj <- nlevels(subj)
    term_df <- ifelse(between, n_subj - sum(between),
                      nrow(data) - n_subj - sum(!between))
    term_df <- pmax(term_df, 1)
    pvals <- 2 * stats::pt(-abs(tt[, "t-value"]), df = term_df)
    coefs <- data.frame(term = rownames(tt), estimate = tt[, "Value"],
                        std_error = tt[, "Std.Error"], df = unname(term_df),
                        statistic = tt[, "t-value"], p_value = unname(pvals),
                        row.names = NULL, stringsAsFactors = FALSE)
    backend <- "gls"
    message("fitMMRM: mmrm not installed; using an nlme::gls approximation with ",
            "between-within degrees of freedom. Install 'mmrm' for the ",
            "unstructured Kenward-Roger/Satterthwaite analysis.")
  }

  est <- stats::setNames(coefs$estimate, coefs$term)
  PhysioCore::AnalysisResult(
    type = "mmrm", estimate = est,
    method = sprintf("MMRM (%s, %s df, %s)", covariance, df, backend),
    result = list(coefficients = coefs, backend = backend, fit = fit,
                  formula = form, covariance = covariance),
    parameters = list(response = response, treatment = treatment, time = time,
                      subject = subject, covariance = covariance, df = df),
    provenance = data.frame(step = "fitMMRM", timestamp = NA_character_,
                            stringsAsFactors = FALSE))
}
