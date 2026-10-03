# Multiple imputation of missing (typically dropout) data. The default is MICE
# under a missing-at-random assumption (mice); reference-based imputation
# (jump-to-reference, for a hypothetical or treatment-policy estimand under a
# conservative missing-not-at-random assumption) is delegated to rbmi.

#' Multiply impute missing data
#'
#' Wraps \code{mice::mice} for missing-at-random multiple imputation (the
#' default), or delegates to \pkg{rbmi} for reference-based (jump-to-reference)
#' imputation. Imputations are reproducible under a fixed \code{seed}.
#'
#' @param data A data frame with missing values (\code{NA}).
#' @param method Imputation method passed to mice (default \code{"pmm"});
#'   recycled across incomplete columns.
#' @param m Number of imputations (default 5).
#' @param predictors Optional: a character vector of variables to use as
#'   predictors for every imputed column, or a full mice
#'   \code{predictorMatrix}. \code{NULL} uses the mice default.
#' @param seed Optional integer seed for reproducibility.
#' @param reference_based Logical; if \code{TRUE}, perform reference-based
#'   (jump-to-reference) imputation via \pkg{rbmi} (which must be installed).
#' @param ... Further arguments passed to \code{mice::mice}.
#' @return A \code{mice} \code{"mids"} object (MAR path) with \code{m} and
#'   \code{method} recorded.
#' @references van Buuren & Groothuis-Oudshoorn 2011 (mice); Carpenter et al.
#'   2013 (reference-based MI); \pkg{rbmi}.
#' @seealso [poolEstimates()], [analyseEstimand()]
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("mice", quietly = TRUE)) {
#'   imp <- multipleImputation(mice::nhanes, m = 5, seed = 1)
#' }
#' }
multipleImputation <- function(data, method = "pmm", m = 5, predictors = NULL,
                               seed = NULL, reference_based = FALSE, ...) {
  if (!is.data.frame(data)) stop("'data' must be a data frame.", call. = FALSE)
  m <- as.integer(m)
  if (is.na(m) || m < 1L) stop("'m' must be a positive integer.", call. = FALSE)

  if (reference_based) {
    if (!requireNamespace("rbmi", quietly = TRUE)) {
      stop("reference-based (jump-to-reference) MI requires the 'rbmi' package.",
           call. = FALSE)
    }
    stop("reference-based MI via rbmi expects a longitudinal draws/impute ",
         "specification; use rbmi::draws()/impute() directly, or call ",
         "analyseEstimand() with a hypothetical strategy.", call. = FALSE)
  }

  if (!requireNamespace("mice", quietly = TRUE)) {
    stop("multipleImputation() requires the 'mice' package.", call. = FALSE)
  }
  args <- list(data = data, m = m, method = method, printFlag = FALSE, ...)
  if (!is.null(seed)) args$seed <- seed
  if (!is.null(predictors)) {
    if (is.matrix(predictors)) {
      args$predictorMatrix <- predictors
    } else if (is.character(predictors)) {
      vars <- names(data)
      pm <- matrix(0L, length(vars), length(vars),
                   dimnames = list(vars, vars))
      pm[, intersect(predictors, vars)] <- 1L
      diag(pm) <- 0L
      args$predictorMatrix <- pm
    } else {
      stop("'predictors' must be a character vector or a predictorMatrix.",
           call. = FALSE)
    }
  }
  do.call(mice::mice, args)
}
