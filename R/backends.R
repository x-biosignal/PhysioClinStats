#' Which optional modelling backends are available
#'
#' The inference functions load heavy modelling packages only when present
#' (guarded via \code{requireNamespace}). This reports availability so callers
#' can fail early with a helpful message.
#'
#' @return A named logical vector for \code{mmrm}, \code{emmeans}, \code{lme4},
#'   and \code{SingleCaseES}.
#' @examples
#' clinStatsBackends()
#' @export
clinStatsBackends <- function() {
  pkgs <- c("mmrm", "emmeans", "lme4", "SingleCaseES")
  vapply(pkgs, function(p) requireNamespace(p, quietly = TRUE), logical(1))
}

#' Require an optional backend or stop with guidance
#'
#' @param pkg Backend package name.
#' @return Invisibly \code{TRUE} if available; otherwise an error.
#' @export
requireBackend <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("This function needs the optional '%s' package. Install it with install.packages('%s').",
                 pkg, pkg), call. = FALSE)
  }
  invisible(TRUE)
}
