# Multiply impute missing data

Wraps [`mice::mice`](https://amices.org/mice/reference/mice.html) for
missing-at-random multiple imputation (the default), or delegates to
rbmi for reference-based (jump-to-reference) imputation. Imputations are
reproducible under a fixed `seed`.

## Usage

``` r
multipleImputation(
  data,
  method = "pmm",
  m = 5,
  predictors = NULL,
  seed = NULL,
  reference_based = FALSE,
  ...
)
```

## Arguments

- data:

  A data frame with missing values (`NA`).

- method:

  Imputation method passed to mice (default `"pmm"`); recycled across
  incomplete columns.

- m:

  Number of imputations (default 5).

- predictors:

  Optional: a character vector of variables to use as predictors for
  every imputed column, or a full mice `predictorMatrix`. `NULL` uses
  the mice default.

- seed:

  Optional integer seed for reproducibility.

- reference_based:

  Logical; if `TRUE`, perform reference-based (jump-to-reference)
  imputation via rbmi (which must be installed).

- ...:

  Further arguments passed to
  [`mice::mice`](https://amices.org/mice/reference/mice.html).

## Value

A `mice` `"mids"` object (MAR path) with `m` and `method` recorded.

## References

van Buuren & Groothuis-Oudshoorn 2011 (mice); Carpenter et al. 2013
(reference-based MI); rbmi.

## See also

[`poolEstimates()`](https://x-biosignal.github.io/PhysioClinStats/reference/poolEstimates.md),
[`analyseEstimand()`](https://x-biosignal.github.io/PhysioClinStats/reference/analyseEstimand.md)

## Examples

``` r
# \donttest{
if (requireNamespace("mice", quietly = TRUE)) {
  imp <- multipleImputation(mice::nhanes, m = 5, seed = 1)
}
# }
```
