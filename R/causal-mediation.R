# Guarded causal-mediation workflow.

.causal_model <- function(model, argument) {
  if (inherits(model, "AnalysisResult")) {
    model <- PhysioCore::resultValue(model)$fit
    if (is.null(model)) {
      stop(sprintf(
        "'%s' is an AnalysisResult without result$fit; supply a fitted model.",
        argument), call. = FALSE)
    }
  }
  frame <- tryCatch(stats::model.frame(model), error = function(e) NULL)
  formula <- tryCatch(stats::formula(model), error = function(e) NULL)
  if (is.null(frame) || is.null(formula)) {
    stop(sprintf(
      "'%s' must be a supported fitted model that retains its model frame.",
      argument), call. = FALSE)
  }
  list(fit = model, frame = frame, formula = formula)
}

.validate_covariate_values <- function(covariates, frames) {
  if (is.null(covariates)) {
    return(invisible(TRUE))
  }
  if (!is.list(covariates) || !length(covariates) ||
      is.null(names(covariates)) || any(!nzchar(names(covariates))) ||
      anyDuplicated(names(covariates))) {
    stop("'covariates' must be NULL or a uniquely named, non-empty list.",
         call. = FALSE)
  }
  for (nm in names(covariates)) {
    observed <- NULL
    for (frame in frames) {
      if (nm %in% names(frame)) {
        observed <- frame[[nm]]
        break
      }
    }
    if (is.null(observed)) {
      stop(sprintf("Covariate '%s' is absent from both model frames.", nm),
           call. = FALSE)
    }
    value <- covariates[[nm]]
    if (length(value) != 1L || is.na(value)) {
      stop(sprintf("Covariate '%s' must have one non-missing value.", nm),
           call. = FALSE)
    }
    if (is.numeric(value) && !is.finite(value)) {
      stop(sprintf("Covariate '%s' must be finite.", nm), call. = FALSE)
    }
    if (is.factor(observed) && !as.character(value) %in% levels(observed)) {
      stop(sprintf("Covariate '%s' value is not a fitted factor level.", nm),
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

.with_preserved_rng <- function(seed, expr) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  force(expr)
}

.mediate_effect_table <- function(object) {
  mapping <- list(
    acme_control = c("d0", "d0.ci", "d0.p"),
    acme_treated = c("d1", "d1.ci", "d1.p"),
    acme_average = c("d.avg", "d.avg.ci", "d.avg.p"),
    ade_control = c("z0", "z0.ci", "z0.p"),
    ade_treated = c("z1", "z1.ci", "z1.p"),
    ade_average = c("z.avg", "z.avg.ci", "z.avg.p"),
    total_effect = c("tau.coef", "tau.ci", "tau.p"),
    proportion_mediated_control = c("n0", "n0.ci", "n0.p"),
    proportion_mediated_treated = c("n1", "n1.ci", "n1.p"),
    proportion_mediated_average = c("n.avg", "n.avg.ci", "n.avg.p")
  )
  rows <- lapply(names(mapping), function(effect) {
    fields <- mapping[[effect]]
    missing <- fields[!fields %in% names(object)]
    if (length(missing)) {
      stop(sprintf("mediation backend omitted expected field(s): %s.",
                   paste(missing, collapse = ", ")), call. = FALSE)
    }
    ci <- as.numeric(object[[fields[2L]]])
    data.frame(
      effect = effect,
      estimate = as.numeric(object[[fields[1L]]])[1L],
      conf_low = ci[1L],
      conf_high = ci[2L],
      p_value = as.numeric(object[[fields[3L]]])[1L],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.medsens_result <- function(object, supported, reason = NULL) {
  if (!supported) {
    return(list(status = "not_applicable", reason = reason))
  }
  sensitivity <- object
  curves <- data.frame(
    rho = as.numeric(sensitivity$rho),
    acme_control = as.numeric(sensitivity$d0),
    acme_treated = as.numeric(sensitivity$d1),
    ade_control = as.numeric(sensitivity$z0),
    ade_treated = as.numeric(sensitivity$z1),
    stringsAsFactors = FALSE
  )
  list(
    status = "computed",
    curves = curves,
    zero_crossing = list(
      acme = sensitivity$err.cr.d,
      ade = sensitivity$err.cr.z
    ),
    settings = list(
      rho.by = sensitivity$rho.by,
      sims = sensitivity$sims,
      effect.type = sensitivity$effect.type
    ),
    backend = sensitivity
  )
}

#' Model-based causal mediation analysis
#'
#' Calls [mediation::mediate()] after checking that the mediator and outcome
#' models describe the same ordered observations. Estimates rely on sequential
#' ignorability, consistency, positivity, absence of a treatment-induced
#' mediator-outcome confounder, and correct mediator and outcome models. An
#' indirect association, even when statistically significant, is not evidence
#' that a biological mechanism has been established.
#'
#' @param model_m,model_y Supported fitted mediator and outcome models, or
#'   `AnalysisResult` objects whose `result$fit` contains those models.
#' @param treat,mediator Treatment and mediator variable names.
#' @param control.value,treat.value Treatment contrast supplied to the backend.
#' @param covariates Optional named list fixing covariates for a conditional
#'   mediation estimand.
#' @param sims Positive integer number of simulations or bootstrap replicates.
#' @param boot Use the nonparametric bootstrap.
#' @param boot.ci.type Bootstrap confidence-interval type.
#' @param conf.level Confidence level in `(0, 1)`.
#' @param robustSE Request heteroskedasticity-consistent uncertainty for
#'   supported `lm`/`glm` models.
#' @param sensitivity Run [mediation::medsens()] for an all-linear model pair.
#' @param seed Optional integer seed. The caller's global RNG state is restored.
#' @param ... Further arguments passed to [mediation::mediate()].
#' @return An `AnalysisResult` with a complete ACME/ADE/total/proportion table,
#'   the backend object, settings, formulas, assumptions, and sensitivity
#'   diagnostics.
#' @export
causalMediation <- function(
    model_m,
    model_y,
    treat,
    mediator,
    control.value = 0,
    treat.value = 1,
    covariates = NULL,
    sims = 1000L,
    boot = FALSE,
    boot.ci.type = "perc",
    conf.level = 0.95,
    robustSE = FALSE,
    sensitivity = TRUE,
    seed = NULL,
    ...) {
  requireBackend("mediation")
  treat <- .causal_scalar_string(treat, "treat")
  mediator <- .causal_scalar_string(mediator, "mediator")
  mm <- .causal_model(model_m, "model_m")
  my <- .causal_model(model_y, "model_y")

  if (!treat %in% names(mm$frame) || !treat %in% names(my$frame)) {
    stop(sprintf("Treatment '%s' must occur in both model frames.", treat),
         call. = FALSE)
  }
  mediator_response <- all.vars(mm$formula)[1L]
  if (!identical(mediator_response, mediator)) {
    stop(sprintf("Mediator '%s' must be the mediator-model response.", mediator),
         call. = FALSE)
  }
  if (!mediator %in% names(my$frame)) {
    stop(sprintf("Mediator '%s' must occur in the outcome model.", mediator),
         call. = FALSE)
  }
  if (!identical(row.names(mm$frame), row.names(my$frame)) ||
      !identical(as.character(mm$frame[[treat]]),
                 as.character(my$frame[[treat]])) ||
      !isTRUE(all.equal(as.vector(mm$frame[[mediator]]),
                        as.vector(my$frame[[mediator]]),
                        check.attributes = FALSE))) {
    stop("Mediator and outcome model rows must identify the same observations in the same order.",
         call. = FALSE)
  }
  observed <- unique(mm$frame[[treat]])
  if (length(control.value) != 1L || is.na(control.value) ||
      length(treat.value) != 1L || is.na(treat.value)) {
    stop("'control.value' and 'treat.value' must be non-missing scalars.",
         call. = FALSE)
  }
  if (length(observed) < 2L ||
      !as.character(control.value) %in% as.character(observed) ||
      !as.character(treat.value) %in% as.character(observed) ||
      identical(as.character(control.value), as.character(treat.value))) {
    stop("'control.value' and 'treat.value' must be distinct observed treatment levels.",
         call. = FALSE)
  }
  if (length(observed) > 2L &&
      (missing(control.value) || missing(treat.value))) {
    stop("Non-binary treatment requires explicit 'control.value' and 'treat.value'.",
         call. = FALSE)
  }
  .validate_covariate_values(covariates, list(mm$frame, my$frame))

  if (length(sims) != 1L || is.na(sims) || !is.finite(sims) ||
      sims < 1 || sims > .Machine$integer.max || sims != floor(sims)) {
    stop("'sims' must be one positive integer.", call. = FALSE)
  }
  sims <- as.integer(sims)
  if (length(conf.level) != 1L || !is.finite(conf.level) ||
      conf.level <= 0 || conf.level >= 1) {
    stop("'conf.level' must be finite and strictly between 0 and 1.",
         call. = FALSE)
  }
  if (!is.null(seed) && (length(seed) != 1L || is.na(seed) ||
                         !is.finite(seed) ||
                         abs(seed) > .Machine$integer.max ||
                         seed != trunc(seed))) {
    stop("'seed' must be NULL or one integer.", call. = FALSE)
  }
  if (!is.null(seed)) {
    seed <- as.integer(seed)
  }
  if (!is.logical(boot) || length(boot) != 1L || is.na(boot) ||
      !is.logical(robustSE) || length(robustSE) != 1L || is.na(robustSE) ||
      !is.logical(sensitivity) || length(sensitivity) != 1L ||
      is.na(sensitivity)) {
    stop("'boot', 'robustSE', and 'sensitivity' must be single logical values.",
         call. = FALSE)
  }
  if (!is.character(boot.ci.type) || length(boot.ci.type) != 1L ||
      is.na(boot.ci.type) || !boot.ci.type %in% c("perc", "bca")) {
    stop("'boot.ci.type' must be 'perc' or 'bca'.", call. = FALSE)
  }
  if (boot && robustSE) {
    stop("'robustSE = TRUE' is not available with bootstrap inference.",
         call. = FALSE)
  }
  robust_supported <- inherits(mm$fit, c("lm", "glm")) &&
    inherits(my$fit, c("lm", "glm"))
  if (robustSE && !robust_supported) {
    stop("'robustSE = TRUE' is supported only for lm/glm model pairs.",
         call. = FALSE)
  }

  mediate_args <- c(list(
    model.m = mm$fit, model.y = my$fit, treat = treat, mediator = mediator,
    covariates = covariates, sims = sims, boot = boot,
    boot.ci.type = boot.ci.type, conf.level = conf.level,
    control.value = control.value, treat.value = treat.value,
    robustSE = robustSE
  ), list(...))
  sensitivity_supported <- inherits(mm$fit, "lm") &&
    !inherits(mm$fit, "glm") && inherits(my$fit, "lm") &&
    !inherits(my$fit, "glm")
  computed <- .with_preserved_rng(seed, {
    backend <- do.call(mediation::mediate, mediate_args)
    sensitivity_backend <- if (sensitivity && sensitivity_supported) {
      mediation::medsens(backend, effect.type = "both")
    } else {
      NULL
    }
    list(backend = backend, sensitivity = sensitivity_backend)
  })
  backend <- computed$backend
  effects <- .mediate_effect_table(backend)

  sensitivity_result <- if (!sensitivity) {
    list(status = "not_requested", reason = "sensitivity = FALSE")
  } else {
    .medsens_result(
      computed$sensitivity, sensitivity_supported,
      "medsens requires linear mediator and outcome models"
    )
  }
  primary_names <- c("acme_average", "ade_average", "total_effect")
  primary <- stats::setNames(
    effects$estimate[match(primary_names, effects$effect)], primary_names)
  primary_rows <- effects[match(primary_names, effects$effect), ]
  assumptions <- c(
    "sequential ignorability",
    "consistency",
    "positivity",
    "no treatment-induced mediator-outcome confounder",
    "correct mediator and outcome model specification"
  )

  PhysioCore::AnalysisResult(
    type = "causal_mediation",
    estimate = primary,
    uncertainty = list(
      type = if (boot) "bootstrap" else "analytic",
      level = conf.level,
      lower = stats::setNames(primary_rows$conf_low, primary_names),
      upper = stats::setNames(primary_rows$conf_high, primary_names)
    ),
    method = "model-based causal mediation (mediation::mediate)",
    estimand = list(
      treatment = sprintf("%s versus %s", treat.value, control.value),
      mediator = mediator,
      effects = "average ACME, average ADE, and total effect",
      population = "observations retained in both fitted models",
      assumptions = assumptions,
      interpretation = paste(
        "ACME/ADE estimates depend on the stated identification assumptions;",
        "a statistically significant indirect effect does not prove a",
        "biological mechanism."
      )
    ),
    result = list(
      effects = effects,
      sensitivity = sensitivity_result,
      backend = backend,
      settings = list(
        sims = sims, boot = boot, boot.ci.type = boot.ci.type,
        conf.level = conf.level, robustSE = robustSE, seed = seed
      ),
      treatment = list(
        variable = treat, control = control.value, treated = treat.value
      ),
      mediator = mediator,
      sample_size = nrow(mm$frame),
      formulas = list(
        mediator = paste(deparse(mm$formula), collapse = " "),
        outcome = paste(deparse(my$formula), collapse = " ")
      ),
      assumptions = assumptions
    ),
    parameters = list(
      treat = treat, mediator = mediator, control.value = control.value,
      treat.value = treat.value, covariates = covariates
    ),
    provenance = data.frame(
      step = "causalMediation",
      backend = "mediation",
      backend_version = as.character(utils::packageVersion("mediation")),
      seed = if (is.null(seed)) NA_integer_ else seed,
      timestamp = NA_character_,
      stringsAsFactors = FALSE
    )
  )
}
