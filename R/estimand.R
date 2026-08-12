# ICH E9(R1) estimand framework. An estimand is defined by five attributes -
# treatment, population, endpoint (variable), the handling of each intercurrent
# event (one of five strategies), and the population-level summary - which
# together pin down "what is being estimated" before any analysis is chosen.

#' Valid ICH E9(R1) intercurrent-event strategies
#' @export
ESTIMAND_STRATEGIES <- c("treatment-policy", "hypothetical", "composite",
                         "while-on-treatment", "principal-stratum")

# so an "estimand" object satisfies the (S4) "list"-typed AnalysisResult slot
methods::setOldClass(c("estimand", "list"))

# each strategy's analysis recipe (documented mapping used by analyseEstimand)
.ESTIMAND_RECIPES <- list(
  "treatment-policy" =
    "Use all observed data regardless of the intercurrent event; standard MAR multiple imputation.",
  "hypothetical" =
    "Estimate the value that would have occurred absent the event; impute post-event data under the no-event (MAR) model.",
  "composite" =
    "Incorporate the intercurrent event into the endpoint (e.g. a non-responder / failure category).",
  "while-on-treatment" =
    "Restrict the endpoint to the period before the intercurrent event.",
  "principal-stratum" =
    "Estimate the effect within the latent subpopulation defined by potential intercurrent-event status.")

#' Define an ICH E9(R1) estimand
#'
#' Constructs a validated estimand from its five ICH E9(R1) attributes. The
#' returned object is a fully-named list (class \code{"estimand"}) that
#' round-trips through the \code{estimand} slot of an
#' \code{\link[PhysioCore]{AnalysisResult}}.
#'
#' @param treatment The treatment condition(s) being compared.
#' @param population The target population.
#' @param endpoint The endpoint / variable of interest.
#' @param intercurrent A named list with \code{event} (the intercurrent event)
#'   and \code{strategy} (one of [ESTIMAND_STRATEGIES]).
#' @param summary_measure The population-level summary (e.g.
#'   \code{"difference in means"}).
#' @return An \code{"estimand"} object (a named list of the five attributes plus
#'   the strategy's analysis \code{recipe}).
#' @references ICH E9(R1) addendum on estimands and sensitivity analysis (2019).
#' @seealso [analyseEstimand()], [multipleImputation()]
#' @importFrom methods setOldClass
#' @export
#' @examples
#' defineEstimand(
#'   treatment = "rehab protocol A vs B", population = "post-stroke",
#'   endpoint = "6-month gait speed",
#'   intercurrent = list(event = "treatment discontinuation",
#'                       strategy = "treatment-policy"),
#'   summary_measure = "difference in means")
defineEstimand <- function(treatment, population, endpoint,
                           intercurrent = list(event = NA_character_,
                                               strategy = "treatment-policy"),
                           summary_measure = "difference in means") {
  if (!is.list(intercurrent) || is.null(intercurrent$strategy)) {
    stop("'intercurrent' must be a list with an 'event' and a 'strategy'.",
         call. = FALSE)
  }
  strategy <- intercurrent$strategy
  if (length(strategy) != 1L || !strategy %in% ESTIMAND_STRATEGIES) {
    stop("invalid intercurrent-event 'strategy': ", strategy,
         ". Must be one of: ", paste(ESTIMAND_STRATEGIES, collapse = ", "), ".",
         call. = FALSE)
  }
  structure(
    list(treatment = treatment, population = population, endpoint = endpoint,
         intercurrent_event = intercurrent$event %||% NA_character_,
         strategy = strategy, summary_measure = summary_measure,
         recipe = .ESTIMAND_RECIPES[[strategy]]),
    class = "estimand")
}

#' @export
print.estimand <- function(x, ...) {
  cat("<estimand> ICH E9(R1)\n")
  cat("  treatment:  ", x$treatment, "\n")
  cat("  population: ", x$population, "\n")
  cat("  endpoint:   ", x$endpoint, "\n")
  cat("  intercurrent event:", x$intercurrent_event, "\n")
  cat("  strategy:   ", x$strategy, "\n")
  cat("  summary:    ", x$summary_measure, "\n")
  invisible(x)
}
