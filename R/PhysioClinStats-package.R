#' PhysioClinStats: clinical inference engine
#'
#' Mixed-effects / MMRM longitudinal models, single-case (N-of-1) designs,
#' estimands with multiple imputation, causal mediation, declared target-trial
#' emulation, and per-estimate uncertainty. Heavy modelling backends are
#' optional, guarded Suggests.
#'
#' @keywords internal
"_PACKAGE"

# `.data` is the tidy-evaluation pronoun used inside the (Suggests-only) ggplot2
# aes() calls of plot.sced_abab; declare it so R CMD check does not read it as an
# undefined global.
utils::globalVariables(".data")
