# Latent-class growth analysis (LCGA): a finite mixture of trajectory shapes.
# Each subject is assigned to one latent class with a class-specific growth curve
# (flexmix), the class count is chosen by BIC, and class separation is summarised
# by the normalised entropy.

#' Latent-class growth analysis of longitudinal trajectories
#'
#' Fits a finite mixture of polynomial growth trajectories to panel data (each
#' subject belongs to one latent class), selecting the number of classes by BIC.
#' Returns per-subject class assignments, posterior class probabilities, and the
#' normalised entropy of the classification.
#'
#' @param data A long-format data frame.
#' @param subject,time,outcome Column names for the grouping factor, the time
#'   variable, and the response.
#' @param n_classes Integer vector of class counts to compare (default
#'   \code{1:3}); the BIC-minimising count is selected.
#' @param degree Polynomial degree of the per-class growth curve (default 1,
#'   linear).
#' @param seed Optional RNG seed (flexmix uses random starts).
#' @return An \code{AnalysisResult} (type \code{"latent_class_growth"}) whose
#'   \code{result} holds the selected \code{n_classes}, the \code{bic} table, the
#'   per-subject \code{assignment}, the \code{posterior} probabilities, the
#'   \code{entropy}, and the fitted flexmix \code{model}.
#' @references Nagin 2005 (group-based trajectory modelling); Gruen & Leisch 2008
#'   (flexmix).
#' @seealso [recoveryTrajectoryLME()], [proportionalRecoveryRule()]
#' @export
#' @examples
#' set.seed(1)
#' df <- do.call(rbind, lapply(1:40, function(s) {
#'   fast <- s <= 20; t <- 0:6
#'   y <- (if (fast) 5 * t else 0.5 * t) + rnorm(7, 0, 1.5)
#'   data.frame(subject = s, time = t, y = y)
#' }))
#' latentClassGrowth(df, "subject", "time", "y", n_classes = 1:3)
latentClassGrowth <- function(data, subject, time, outcome, n_classes = 1:3,
                              degree = 1, seed = NULL) {
  if (!requireNamespace("flexmix", quietly = TRUE)) {
    stop("latentClassGrowth() requires the 'flexmix' package.", call. = FALSE)
  }
  for (col in c(subject, time, outcome)) {
    if (!col %in% names(data)) stop("column '", col, "' not in data.",
                                    call. = FALSE)
  }
  d <- data.frame(.subj = factor(data[[subject]]),
                  .t = as.numeric(data[[time]]),
                  .y = as.numeric(data[[outcome]]))
  d <- d[stats::complete.cases(d), , drop = FALSE]
  n_classes <- sort(unique(as.integer(n_classes)))
  if (any(n_classes < 1L)) stop("'n_classes' must be >= 1.", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)

  rhs <- if (degree > 1L) sprintf("poly(.t, %d)", degree) else ".t"
  form <- stats::as.formula(sprintf(".y ~ %s | .subj", rhs))

  fits <- list(); bics <- rep(NA_real_, length(n_classes))
  for (i in seq_along(n_classes)) {
    k <- n_classes[i]
    m <- tryCatch(
      flexmix::flexmix(form, data = d, k = k,
                       control = list(minprior = 0.05)),
      error = function(e) NULL)
    # keep only models that did not collapse to fewer components
    if (!is.null(m) && length(unique(flexmix::clusters(m))) == k) {
      fits[[as.character(k)]] <- m
      bics[i] <- stats::BIC(m)
    }
  }
  if (!length(fits)) stop("no latent-class model converged.", call. = FALSE)

  best_k <- n_classes[which.min(bics)]
  best <- fits[[as.character(best_k)]]
  post <- flexmix::posterior(best)                        # obs x class
  assign_obs <- flexmix::clusters(best)
  # one row per subject (all rows of a subject share a class)
  subj <- d$.subj
  first <- !duplicated(subj)
  assignment <- data.frame(subject = as.character(subj[first]),
                           class = assign_obs[first],
                           row.names = NULL, stringsAsFactors = FALSE)
  post_subj <- post[first, , drop = FALSE]
  colnames(post_subj) <- paste0("class", seq_len(ncol(post_subj)))

  res <- list(
    n_classes = best_k,
    bic = data.frame(k = n_classes, BIC = bics),
    assignment = assignment,
    posterior = as.data.frame(post_subj),
    entropy = .lcga_entropy(post_subj),
    class_sizes = as.integer(table(factor(assignment$class,
                                          levels = seq_len(best_k)))),
    model = best)

  PhysioCore::AnalysisResult(
    type = "latent_class_growth", result = res,
    parameters = list(n_classes = best_k, degree = degree),
    estimate = best_k, method = "flexmix_lcga")
}

# normalised entropy of the posterior classification (1 = perfect separation)
.lcga_entropy <- function(post) {
  K <- ncol(post)
  if (K < 2L) return(1)
  p <- pmax(post, .Machine$double.eps)
  ent <- -sum(post * log(p))
  1 - ent / (nrow(post) * log(K))
}
