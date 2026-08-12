# Linear mixed-effects models for longitudinal / clustered rehabilitation data,
# wrapping lme4 with tidy fixed- and random-effect extraction into the
# ecosystem's AnalysisResult carrier.

# Tidy the fixed-effect table. When the model is an lmerModLmerTest (fitted via
# lmerTest::lmer) the summary carries Satterthwaite / Kenward-Roger denominator
# df and t tests; otherwise fall back to the normal approximation.
.tidy_fixed_lmer <- function(model, ddf) {
  co <- if (inherits(model, "lmerModLmerTest")) {
    stats::coef(summary(model, ddf = ddf))
  } else {
    stats::coef(summary(model))
  }
  has_df <- "df" %in% colnames(co)
  data.frame(
    term = rownames(co), estimate = co[, "Estimate"],
    std_error = co[, "Std. Error"],
    df = if (has_df) co[, "df"] else NA_real_,
    statistic = co[, "t value"],
    p_value = if (has_df) co[, "Pr(>|t|)"] else
      2 * stats::pnorm(-abs(co[, "t value"])),
    row.names = NULL, stringsAsFactors = FALSE)
}

#' Fit a linear mixed-effects model
#'
#' Fits a linear mixed model with `lme4::lmer` and returns tidy fixed-effect and
#' random-effect tables in an `AnalysisResult`. Denominator degrees of freedom
#' (and the associated t tests) use `lmerTest` when it is installed - either the
#' Satterthwaite or the Kenward-Roger approximation.
#'
#' @param data A data.frame in long format.
#' @param formula A model formula. It may carry the random-effects terms
#'   directly (e.g. `y ~ x + (x | subject)`), or supply them via `random`.
#' @param random Optional one-sided random-effects formula (e.g. `~ (1 |
#'   subject)`) appended to `formula` when the latter has none.
#' @param method `"REML"` (default) or `"ML"`.
#' @param df Denominator-df method for the fixed-effect tests: `"satterthwaite"`
#'   (default) or `"kenward-roger"` (both need `lmerTest`).
#' @return An `AnalysisResult` (`type = "mixed_model"`) whose `estimate` is the
#'   fixed-effect coefficient vector, with `result$fixed` (the tidy fixed-effect
#'   table), `result$random` (variance components), `result$sigma`, `result$fit`
#'   (the fitted model, for [estimatedMarginalMeans()]) and `result$formula`.
#' @references Bates D, Maechler M, Bolker B, Walker S (2015). Fitting linear
#'   mixed-effects models using lme4. Journal of Statistical Software, 67(1).
#'   Kuznetsova A, Brockhoff PB, Christensen RHB (2017). lmerTest package.
#'   Journal of Statistical Software, 82(13).
#' @seealso [fitMMRM()], [estimatedMarginalMeans()]
#' @export
#' @examples
#' if (requireNamespace("lme4", quietly = TRUE)) {
#'   data(sleepstudy, package = "lme4")
#'   fitMixedModel(sleepstudy, Reaction ~ Days + (Days | Subject))
#' }
fitMixedModel <- function(data, formula, random = NULL,
                          method = c("REML", "ML"),
                          df = c("satterthwaite", "kenward-roger")) {
  method <- match.arg(method)
  df <- match.arg(df)
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("fitMixedModel needs the lme4 package.", call. = FALSE)
  }
  stopifnot(is.data.frame(data), inherits(formula, "formula"))

  # A random-effects term uses the `( ... | group )` bar syntax; detect it
  # without lme4::findbars (which is being deprecated upstream).
  has_bar <- grepl("\\|", paste(deparse(formula), collapse = " "))
  if (!has_bar) {
    if (is.null(random)) {
      stop("No random-effects term in 'formula'; supply 'random' (e.g. ",
           "~ (1 | subject)).", call. = FALSE)
    }
    rhs <- paste(deparse(random[[length(random)]]), collapse = " ")
    formula <- stats::update(formula, paste("~ . +", rhs))
  }

  ddf <- if (df == "kenward-roger") "Kenward-Roger" else "Satterthwaite"
  # lmerTest::lmer returns an lmerModLmerTest carrying what the Satterthwaite /
  # Kenward-Roger df need; plain lme4::lmer otherwise (normal-approx tests).
  model <- if (requireNamespace("lmerTest", quietly = TRUE)) {
    lmerTest::lmer(formula, data = data, REML = (method == "REML"))
  } else {
    lme4::lmer(formula, data = data, REML = (method == "REML"))
  }
  fixed <- .tidy_fixed_lmer(model, ddf)
  vc <- as.data.frame(lme4::VarCorr(model))
  random <- data.frame(group = vc$grp, term1 = vc$var1, term2 = vc$var2,
                       variance = vc$vcov, sd = vc$sdcor,
                       stringsAsFactors = FALSE)

  est <- stats::setNames(fixed$estimate, fixed$term)
  PhysioCore::AnalysisResult(
    type = "mixed_model", estimate = est, method = method,
    result = list(fixed = fixed, random = random, sigma = stats::sigma(model),
                  fit = model, formula = formula, df_method = ddf),
    parameters = list(method = method, df = df),
    provenance = data.frame(step = "fitMixedModel", timestamp = NA_character_,
                            stringsAsFactors = FALSE))
}

#' Reshape a longitudinal PhysioExperiment to long format
#'
#' Flattens an epoched or longitudinal `PhysioExperiment` into a long
#' `data.frame(subject, time, channel, value)` suitable for [fitMixedModel()],
#' [fitMMRM()] and the marginal-means helpers.
#'
#' @param pe A `PhysioExperiment`. A 2D assay is treated as time x channel for a
#'   single subject; a 3D assay as time x channel x subject/trial.
#' @param assay Assay name (default the first).
#' @param value_name Name for the value column (default `"value"`).
#' @param subjects Optional character/factor labels for the third dimension
#'   (subjects/trials); defaults to the `colnames`/index.
#' @return A long `data.frame` with columns `subject`, `time`, `channel` and the
#'   value column.
#' @seealso [fitMixedModel()], [fitMMRM()]
#' @export
peToLong <- function(pe, assay = NULL, value_name = "value", subjects = NULL) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("peToLong needs the SummarizedExperiment package.", call. = FALSE)
  }
  stopifnot(inherits(pe, "PhysioExperiment"))
  a <- if (is.null(assay)) SummarizedExperiment::assayNames(pe)[1] else assay
  arr <- SummarizedExperiment::assay(pe, a)
  cd <- SummarizedExperiment::colData(pe)
  chan <- if ("label" %in% names(cd)) as.character(cd$label) else
    as.character(seq_len(ncol(arr)))

  d <- dim(arr)
  if (length(d) == 2L) {
    long <- expand.grid(time = seq_len(d[1]), channel = chan,
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    long$subject <- "1"
    long[[value_name]] <- as.vector(arr)
  } else if (length(d) == 3L) {
    subj <- if (!is.null(subjects)) as.character(subjects) else
      as.character(seq_len(d[3]))
    long <- expand.grid(time = seq_len(d[1]), channel = chan,
                        subject = subj, KEEP.OUT.ATTRS = FALSE,
                        stringsAsFactors = FALSE)
    long[[value_name]] <- as.vector(arr)
  } else {
    stop("Assay must be a 2D or 3D array.", call. = FALSE)
  }
  long[, c("subject", "time", "channel", value_name)]
}
