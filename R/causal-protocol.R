# Target-trial protocol declaration and validation.

.causal_scalar_string <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(trimws(x))) {
    stop(sprintf("'%s' must be one non-empty character value.", name),
         call. = FALSE)
  }
  x
}

.causal_serializable <- function(x) {
  tryCatch({
    serialize(x, NULL, version = 2)
    TRUE
  }, error = function(e) FALSE)
}

.causal_repr <- function(x) {
  if (is.function(x)) {
    paste(deparse(x, width.cutoff = 500L, control = "all"), collapse = "\n")
  } else if (inherits(x, "formula")) {
    paste(deparse(x, width.cutoff = 500L), collapse = " ")
  } else {
    paste(utils::capture.output(dput(x)), collapse = "\n")
  }
}

.causal_hash <- function(x) {
  digest::digest(x, algo = "sha256", serialize = TRUE)
}

.nonempty_protocol_component <- function(x) {
  if (is.null(x) || !length(x) || !.causal_serializable(x)) {
    return(FALSE)
  }
  if (is.character(x)) {
    return(all(!is.na(x) & nzchar(trimws(x))))
  }
  TRUE
}

.strategy_spec <- function(strategy, name) {
  if (is.function(strategy)) {
    f <- names(formals(strategy))
    if (length(f) < 2L || !identical(f[1:2], c("data", "history"))) {
      stop(sprintf(
        "Dynamic strategy '%s' must have first arguments (data, history).",
        name), call. = FALSE)
    }
    representation <- .causal_repr(strategy)
    return(list(
      name = name, type = "dynamic", portable = FALSE,
      representation = representation,
      hash = .causal_hash(representation)
    ))
  }
  if (length(strategy) != 1L || !is.atomic(strategy) || is.na(strategy) ||
      is.list(strategy)) {
    stop(sprintf(
      "Static strategy '%s' must be one non-missing atomic treatment value.",
      name), call. = FALSE)
  }
  representation <- .causal_repr(strategy)
  list(
    name = name, type = "static", portable = TRUE,
    representation = representation,
    hash = .causal_hash(representation)
  )
}

#' Declare a target-trial protocol
#'
#' Defines the seven protocol components that must be fixed before emulating a
#' target trial. `time_zero` is a separate alignment rule and is never inferred
#' from the observed data. Static strategies are scalar sustained treatment
#' values. Dynamic strategy functions receive the current row as `data` and
#' the participant's rows through the current time as `history`; they must not
#' access later rows.
#'
#' @param eligibility Eligibility rule. Use a function of baseline data or a
#'   one-sided formula for an executable rule.
#' @param treatment_strategies Named list of at least two static scalar or
#'   dynamic function strategies.
#' @param assignment,outcome,causal_contrast,analysis_plan Non-empty,
#'   serializable protocol declarations.
#' @param time_zero Explicit scalar time value or alignment rule.
#' @param follow_up Non-empty follow-up declaration. A finite positive scalar
#'   is used as the fixed horizon by [targetTrialEmulate()].
#' @param protocol_id Optional stable identifier. When omitted, one is derived
#'   from the protocol declarations.
#' @param version Semantic protocol version.
#' @return An object of class `target_trial_protocol`.
#' @export
#' @examples
#' protocol <- targetTrialProtocol(
#'   eligibility = function(data) data$eligible,
#'   treatment_strategies = list(never = 0, always = 1),
#'   assignment = "cloning at eligibility",
#'   time_zero = 0,
#'   follow_up = 12,
#'   outcome = "binary recovery by week 12",
#'   causal_contrast = "always versus never",
#'   analysis_plan = "clone-censor-weight"
#' )
#' protocol
targetTrialProtocol <- function(
    eligibility,
    treatment_strategies,
    assignment,
    time_zero,
    follow_up,
    outcome,
    causal_contrast,
    analysis_plan,
    protocol_id = NULL,
    version = "1.0.0") {
  components <- list(
    eligibility = eligibility,
    treatment_strategies = treatment_strategies,
    assignment = assignment,
    time_zero = time_zero,
    follow_up = follow_up,
    outcome = outcome,
    causal_contrast = causal_contrast,
    analysis_plan = analysis_plan
  )
  required <- setdiff(names(components), "treatment_strategies")
  invalid <- required[!vapply(components[required],
                              .nonempty_protocol_component, logical(1))]
  if (length(invalid)) {
    stop(sprintf("Protocol component(s) must be non-empty and serializable: %s.",
                 paste(invalid, collapse = ", ")), call. = FALSE)
  }
  if (!is.list(treatment_strategies) || length(treatment_strategies) < 2L ||
      is.null(names(treatment_strategies)) ||
      any(!nzchar(names(treatment_strategies))) ||
      anyDuplicated(names(treatment_strategies))) {
    stop("'treatment_strategies' must be a uniquely named list with at least two strategies.",
         call. = FALSE)
  }
  specs <- Map(.strategy_spec, treatment_strategies,
               names(treatment_strategies))
  names(specs) <- names(treatment_strategies)

  version <- .causal_scalar_string(version, "version")
  if (!grepl("^[0-9]+\\.[0-9]+\\.[0-9]+([+-][0-9A-Za-z.-]+)?$", version)) {
    stop("'version' must use semantic version form, for example '1.0.0'.",
         call. = FALSE)
  }
  if (!is.null(protocol_id)) {
    protocol_id <- .causal_scalar_string(protocol_id, "protocol_id")
  }

  provenance_components <- lapply(components, .causal_repr)
  provenance_components$treatment_strategies <- lapply(
    specs, function(x) x[c("type", "representation", "hash")])
  protocol_hash <- .causal_hash(list(
    version = version, components = provenance_components
  ))
  if (is.null(protocol_id)) {
    protocol_id <- paste0("ttp-", substr(protocol_hash, 1L, 12L))
  }
  report <- data.frame(
    component = c(names(components), "version", "protocol_id"),
    valid = TRUE,
    message = c(rep("declared", length(components)),
                "semantic version", "stable identifier"),
    stringsAsFactors = FALSE
  )

  structure(
    c(components, list(
      protocol_id = protocol_id,
      version = version,
      strategy_specs = specs,
      portable = all(vapply(specs, `[[`, logical(1), "portable")),
      protocol_hash = protocol_hash,
      validation = report
    )),
    class = "target_trial_protocol"
  )
}

#' @export
print.target_trial_protocol <- function(x, ...) {
  cat("Target-trial protocol:", x$protocol_id, "(version", x$version, ")\n")
  cat("Strategies:",
      paste(sprintf("%s [%s]", names(x$strategy_specs),
                    vapply(x$strategy_specs, `[[`, character(1), "type")),
            collapse = ", "), "\n")
  cat("Time zero:", .causal_repr(x$time_zero), "\n")
  cat("Eligibility:", .causal_repr(x$eligibility), "\n")
  cat("Assignment:", .causal_repr(x$assignment), "\n")
  cat("Follow-up:", .causal_repr(x$follow_up), "\n")
  cat("Outcome:", .causal_repr(x$outcome), "\n")
  cat("Causal contrast:", .causal_repr(x$causal_contrast), "\n")
  cat("Analysis plan:", .causal_repr(x$analysis_plan), "\n")
  cat("Portable:", if (isTRUE(x$portable)) "yes" else
    "no (contains a dynamic function strategy)", "\n")
  cat("Validation:", sum(x$validation$valid), "of",
      nrow(x$validation), "components valid\n")
  invisible(x)
}
