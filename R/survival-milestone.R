# Time-to-milestone: turn a longitudinal outcome that must cross a clinical
# threshold into a right-censored survival endpoint (event = first crossing,
# censoring = never attained).

#' Derive a time-to-milestone survival endpoint
#'
#' From a long-format longitudinal series, finds for each subject the first time
#' its outcome crosses a clinical `threshold`, producing a right-censored
#' `(time, event)` pair per subject: `event = 1` at the first crossing, or
#' `event = 0` (right-censored) at the last observation - or at `censor_at` - for
#' subjects who never attain the milestone.
#'
#' @param data A long-format data.frame with one row per subject-visit.
#' @param id,time,value Column names of the subject id, the visit time and the
#'   outcome value.
#' @param threshold The clinical milestone value.
#' @param direction `"increase"` (default; the milestone is reaching a value
#'   `>= threshold`) or `"decrease"` (a value `<= threshold`).
#' @param censor_at Optional common administrative censoring time for
#'   non-attainers; defaults to each subject's last observation time.
#' @return A data.frame with columns `id`, `time`, `event` (1 attained /
#'   0 censored) and `attained`, carrying a `Surv` object as the `"surv"`
#'   attribute.
#' @seealso [survivalFit()], [milestoneHazard()]
#' @export
#' @examples
#' d <- data.frame(
#'   id = rep(c("a", "b"), each = 3), time = rep(1:3, 2),
#'   value = c(10, 30, 55, 12, 20, 28))
#' timeToMilestone(d, "id", "time", "value", threshold = 50)
timeToMilestone <- function(data, id, time, value, threshold,
                            direction = c("increase", "decrease"),
                            censor_at = NULL) {
  direction <- match.arg(direction)
  stopifnot(is.data.frame(data),
            all(c(id, time, value) %in% names(data)),
            is.numeric(threshold), length(threshold) == 1)
  crossed <- if (direction == "increase") {
    function(v) v >= threshold
  } else {
    function(v) v <= threshold
  }

  ids <- unique(data[[id]])
  out <- lapply(ids, function(i) {
    sub <- data[data[[id]] == i, , drop = FALSE]
    sub <- sub[order(sub[[time]]), , drop = FALSE]
    tv <- sub[[time]]; vv <- sub[[value]]
    keep <- !is.na(tv) & !is.na(vv)
    tv <- tv[keep]; vv <- vv[keep]
    if (!length(tv)) {
      return(data.frame(id = i, time = NA_real_, event = 0L,
                        attained = FALSE, stringsAsFactors = FALSE))
    }
    hit <- which(crossed(vv))
    if (length(hit)) {
      data.frame(id = i, time = tv[hit[1]], event = 1L, attained = TRUE,
                 stringsAsFactors = FALSE)
    } else {
      ct <- if (is.null(censor_at)) max(tv) else censor_at
      data.frame(id = i, time = ct, event = 0L, attained = FALSE,
                 stringsAsFactors = FALSE)
    }
  })
  res <- do.call(rbind, out)
  names(res)[1] <- "id"
  if (requireNamespace("survival", quietly = TRUE)) {
    attr(res, "surv") <- survival::Surv(res$time, res$event)
  }
  res
}
