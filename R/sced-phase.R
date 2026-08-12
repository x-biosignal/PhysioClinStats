# ABAB phase-design orchestration: pull the phases out of a PhysioExperiment
# single-channel series or a plain data.frame, then run the non-overlap effect
# sizes, 2-SD flags and celeration lines for each reversal contrast.

# Extract a tidy data.frame(time, value, phase) from the accepted inputs.
.sced_extract <- function(x, value = NULL, phase = NULL, time = NULL,
                          assay = NULL, channel = 1L) {
  if (inherits(x, "PhysioExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE) ||
        !requireNamespace("S4Vectors", quietly = TRUE)) {
      stop("Reading a PhysioExperiment needs the SummarizedExperiment and ",
           "S4Vectors packages.", call. = FALSE)
    }
    a <- if (is.null(assay)) SummarizedExperiment::assayNames(x)[1] else assay
    m <- SummarizedExperiment::assay(x, a)
    v <- as.numeric(m[, channel])
    md <- S4Vectors::metadata(x)
    ph <- phase %||% md$phase %||%
      SummarizedExperiment::rowData(x)[["phase"]]
    if (is.null(ph)) stop("No phase labels found; supply 'phase'.", call. = FALSE)
    tm <- time %||% seq_along(v)
    out <- data.frame(time = .sced_numeric(tm, "time"),
                      value = .sced_numeric(v, "value"),
                      phase = as.character(ph), stringsAsFactors = FALSE)
  } else if (is.data.frame(x)) {
    vcol <- value %||% "value"; pcol <- phase %||% "phase"
    if (!all(c(vcol, pcol) %in% names(x))) {
      stop("data.frame needs 'value' and 'phase' columns (or name them).",
           call. = FALSE)
    }
    tcol <- time %||% (if ("time" %in% names(x)) "time" else NULL)
    tm <- if (is.null(tcol)) seq_len(nrow(x)) else x[[tcol]]
    out <- data.frame(time = .sced_numeric(tm, "time"),
                      value = .sced_numeric(x[[vcol]], "value"),
                      phase = as.character(x[[pcol]]), stringsAsFactors = FALSE)
  } else {
    stop("x must be a PhysioExperiment or a data.frame.", call. = FALSE)
  }
  # Order-dependent statistics (the Tau-U baseline-trend correction and the 2-SD
  # consecutive-run rule) must see observations in time order, which the plot
  # method also uses. A stable sort keeps ties in their given order.
  out[order(out$time), , drop = FALSE]
}
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Analyse an ABAB (reversal) single-case design
#'
#' Orchestrates the SCED analysis of a reversal design: it separates the phases
#' (by default the four phases `A1, B1, A2, B2`), and for each intervention
#' contrast (`B1` vs `A1` and `B2` vs `A2`) computes the non-overlap effect
#' sizes, the 2-SD-band decision and per-phase celeration lines.
#'
#' @param x A `PhysioExperiment` single-channel series, or a
#'   `data.frame(time, value, phase)`. Observations are ordered by `time` before
#'   analysis, so a `time` column need not be pre-sorted; without one, row order
#'   is used. Ordering matters because the Tau-U baseline-trend correction and
#'   the 2-SD consecutive-run rule are sequence-dependent.
#' @param phases Ordered character vector of the phase labels (default
#'   `c("A1","B1","A2","B2")`).
#' @param contrasts A list of length-2 `c(baseline, intervention)` phase-label
#'   pairs (default the two reversal contrasts).
#' @param improvement `"increase"` (default) or `"decrease"`.
#' @param value,phase,time,assay,channel Column/assay selectors passed to the
#'   extractor when `x` is a data.frame or PhysioExperiment.
#' @return An object of class `sced_abab` (a list) with `phase_data`,
#'   `contrasts` (per contrast: NAP/Tau/Tau-U/PND/PEM results, the 2-SD decision
#'   and both phase celeration lines) and a `summary` data.frame, with `print`
#'   and `plot` methods.
#' @seealso [scedNAP()], [scedTwoSDBand()], [scedCelerationLine()]
#' @export
#' @examples
#' df <- data.frame(
#'   value = c(10, 11, 9, 18, 19, 20, 10, 12, 11, 21, 22, 20),
#'   phase = rep(c("A1", "B1", "A2", "B2"), each = 3))
#' res <- scedABAB(df)
#' res
scedABAB <- function(x, phases = c("A1", "B1", "A2", "B2"),
                     contrasts = list(c("A1", "B1"), c("A2", "B2")),
                     improvement = c("increase", "decrease"),
                     value = NULL, phase = NULL, time = NULL,
                     assay = NULL, channel = 1L) {
  improvement <- match.arg(improvement)
  d <- .sced_extract(x, value = value, phase = phase, time = time,
                     assay = assay, channel = channel)
  present <- unique(d$phase)
  if (!all(phases %in% present)) {
    stop(sprintf("Phases missing from the data: %s",
                 paste(setdiff(phases, present), collapse = ", ")),
         call. = FALSE)
  }
  split_ph <- lapply(phases, function(p) d[d$phase == p, , drop = FALSE])
  names(split_ph) <- phases

  con <- lapply(contrasts, function(pair) {
    A <- d$value[d$phase == pair[1]]
    B <- d$value[d$phase == pair[2]]
    cel_A <- scedCelerationLine(A, time = d$time[d$phase == pair[1]])
    cel_B <- scedCelerationLine(B, time = d$time[d$phase == pair[2]])
    list(baseline = pair[1], intervention = pair[2],
         nap = scedNAP(A, B, improvement = improvement),
         tau = scedTau(A, B, improvement = improvement),
         tau_u = scedTauU(A, B, improvement = improvement),
         pnd = scedPND(A, B, improvement = improvement),
         pem = scedPEM(A, B, improvement = improvement),
         band = scedTwoSDBand(A, B, improvement = improvement),
         celeration = list(baseline = cel_A, intervention = cel_B))
  })
  names(con) <- vapply(contrasts, function(p) paste(p[2], "vs", p[1]),
                       character(1))

  rv <- function(r) PhysioCore::resultValue(r)
  summ <- do.call(rbind, lapply(con, function(cc) {
    data.frame(contrast = paste(cc$intervention, "vs", cc$baseline),
               NAP = rv(cc$nap)$estimate, Tau = rv(cc$tau)$estimate,
               Tau_U = rv(cc$tau_u)$estimate, PND = rv(cc$pnd)$estimate,
               PEM = rv(cc$pem)$estimate,
               band_flag = rv(cc$band)$estimate,
               p_value = rv(cc$tau_u)$p_value,
               stringsAsFactors = FALSE)
  }))
  rownames(summ) <- NULL

  structure(list(phase_data = d, phases = phases, split = split_ph,
                 contrasts = con, summary = summ, improvement = improvement),
            class = "sced_abab")
}

#' @export
print.sced_abab <- function(x, ...) {
  cat("Single-case ABAB analysis\n")
  cat(sprintf("  phases: %s (%d observations)\n",
              paste(x$phases, collapse = ", "), nrow(x$phase_data)))
  cat(sprintf("  improvement: %s\n\n", x$improvement))
  s <- x$summary
  num <- vapply(s, is.numeric, logical(1))
  s[num] <- lapply(s[num], function(v) round(v, 4))
  print(s, row.names = FALSE)
  invisible(x)
}

#' @export
plot.sced_abab <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot.sced_abab needs the ggplot2 package.", call. = FALSE)
  }
  d <- x$phase_data
  d$phase <- factor(d$phase, levels = x$phases)
  # phase boundaries for vertical separators
  ord <- d[order(d$time), ]
  bnds <- ord$time[which(diff(as.integer(ord$phase)) != 0)] + 0.5
  bands <- do.call(rbind, lapply(x$contrasts, function(cc) {
    b <- PhysioCore::resultValue(cc$band)
    ph <- cc$intervention
    data.frame(phase = ph, mean = b$mean, upper = b$upper, lower = b$lower)
  }))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$time, y = .data$value)) +
    ggplot2::geom_vline(xintercept = bnds, linetype = "dashed",
                        colour = "grey60") +
    ggplot2::geom_line(ggplot2::aes(group = .data$phase)) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$phase)) +
    ggplot2::labs(x = "Session", y = "Outcome", colour = "Phase",
                  title = "ABAB single-case design") +
    ggplot2::theme_minimal()

  # celeration lines per phase
  for (cc in x$contrasts) {
    for (side in c("baseline", "intervention")) {
      ph <- if (side == "baseline") cc$baseline else cc$intervention
      seg <- d[d$phase == ph, ]
      cel <- PhysioCore::resultValue(cc$celeration[[side]])
      p <- p + ggplot2::geom_line(
        data = data.frame(time = seg$time, value = cel$fitted),
        ggplot2::aes(x = .data$time, y = .data$value),
        colour = "steelblue", linewidth = 0.6)
    }
  }
  p
}
