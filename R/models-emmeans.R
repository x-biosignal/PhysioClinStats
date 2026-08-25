# Estimated marginal means and contrasts, wrapping emmeans over the fitted
# mixed / MMRM models.

# Pull the underlying model object out of an AnalysisResult (or accept a raw fit).
.emm_model <- function(model) {
  if (inherits(model, "AnalysisResult")) {
    fit <- PhysioCore::resultValue(model)$fit
    if (is.null(fit)) {
      stop("This AnalysisResult carries no fitted model; pass a fit from ",
           "fitMixedModel()/fitMMRM().", call. = FALSE)
    }
    fit
  } else {
    model
  }
}

#' Estimated marginal means
#'
#' Computes estimated marginal (least-squares) means from a fitted mixed or MMRM
#' model with `emmeans`, optionally with a follow-up contrast, returning the
#' marginal-mean table (and contrasts) in an `AnalysisResult`.
#'
#' @param model An `AnalysisResult` from [fitMixedModel()] / [fitMMRM()], or a
#'   raw model object `emmeans` understands.
#' @param specs A formula or character spec for the marginal means, e.g.
#'   `~ treatment | time`.
#' @param contrasts Optional contrast method passed to `emmeans::contrast`
#'   (e.g. `"pairwise"`, `"trt.vs.ctrl"`); `NULL` (default) returns means only.
#' @param level Confidence level for the intervals (default 0.95).
#' @param ... Further arguments forwarded to `emmeans::emmeans`.
#' @return An `AnalysisResult` (`type = "emmeans"`) with `result$emmeans` (the
#'   marginal-mean table) and, when requested, `result$contrasts`.
#' @references Lenth RV (2024). emmeans: Estimated Marginal Means. R package.
#' @seealso [pairwiseContrasts()], [fitMMRM()]
#' @export
#' @examples
#' if (requireNamespace("mmrm", quietly = TRUE) &&
#'     requireNamespace("emmeans", quietly = TRUE)) {
#'   data(fev_data, package = "mmrm")
#'   fit <- fitMMRM(fev_data, "FEV1", "ARMCD", "AVISIT", "USUBJID",
#'                  covariates = c("RACE", "SEX"))
#'   estimatedMarginalMeans(fit, ~ ARMCD | AVISIT)
#' }
estimatedMarginalMeans <- function(model, specs, contrasts = NULL,
                                   level = 0.95, ...) {
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    stop("estimatedMarginalMeans needs the emmeans package.", call. = FALSE)
  }
  fit <- .emm_model(model)
  emm <- emmeans::emmeans(fit, specs = specs, level = level, ...)
  emm_df <- as.data.frame(summary(emm, level = level))

  con_df <- NULL
  if (!is.null(contrasts)) {
    con <- emmeans::contrast(emm, method = contrasts)
    con_df <- as.data.frame(summary(con, level = level))
  }

  PhysioCore::AnalysisResult(
    type = "emmeans", estimate = emm_df, method = "estimated marginal means",
    result = list(emmeans = emm_df, contrasts = con_df, emm_object = emm),
    parameters = list(specs = paste(deparse(specs), collapse = " "),
                      contrasts = contrasts, level = level),
    provenance = data.frame(step = "estimatedMarginalMeans",
                            timestamp = NA_character_, stringsAsFactors = FALSE))
}

#' Pairwise contrasts of estimated marginal means
#'
#' Convenience wrapper computing all pairwise differences of the estimated
#' marginal means for `specs`.
#'
#' @inheritParams estimatedMarginalMeans
#' @param adjust Multiplicity adjustment passed to `emmeans` (default
#'   `"tukey"`).
#' @return An `AnalysisResult` (`type = "emmeans_contrasts"`) with
#'   `result$contrasts`.
#' @seealso [estimatedMarginalMeans()]
#' @export
pairwiseContrasts <- function(model, specs, adjust = "tukey", level = 0.95, ...) {
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    stop("pairwiseContrasts needs the emmeans package.", call. = FALSE)
  }
  fit <- .emm_model(model)
  emm <- emmeans::emmeans(fit, specs = specs, level = level, ...)
  con <- emmeans::contrast(emm, method = "pairwise", adjust = adjust)
  con_df <- as.data.frame(summary(con, level = level))

  PhysioCore::AnalysisResult(
    type = "emmeans_contrasts", estimate = con_df,
    method = sprintf("pairwise contrasts (%s)", adjust),
    result = list(contrasts = con_df, emm_object = emm),
    parameters = list(specs = paste(deparse(specs), collapse = " "),
                      adjust = adjust, level = level),
    provenance = data.frame(step = "pairwiseContrasts",
                            timestamp = NA_character_, stringsAsFactors = FALSE))
}
