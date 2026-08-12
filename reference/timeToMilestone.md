# Derive a time-to-milestone survival endpoint

From a long-format longitudinal series, finds for each subject the first
time its outcome crosses a clinical `threshold`, producing a
right-censored `(time, event)` pair per subject: `event = 1` at the
first crossing, or `event = 0` (right-censored) at the last
observation - or at `censor_at` - for subjects who never attain the
milestone.

## Usage

``` r
timeToMilestone(
  data,
  id,
  time,
  value,
  threshold,
  direction = c("increase", "decrease"),
  censor_at = NULL
)
```

## Arguments

- data:

  A long-format data.frame with one row per subject-visit.

- id, time, value:

  Column names of the subject id, the visit time and the outcome value.

- threshold:

  The clinical milestone value.

- direction:

  `"increase"` (default; the milestone is reaching a value
  `>= threshold`) or `"decrease"` (a value `<= threshold`).

- censor_at:

  Optional common administrative censoring time for non-attainers;
  defaults to each subject's last observation time.

## Value

A data.frame with columns `id`, `time`, `event` (1 attained / 0
censored) and `attained`, carrying a `Surv` object as the `"surv"`
attribute.

## See also

[`survivalFit()`](https://x-biosignal.github.io/PhysioClinStats/reference/survivalFit.md),
[`milestoneHazard()`](https://x-biosignal.github.io/PhysioClinStats/reference/milestoneHazard.md)

## Examples

``` r
d <- data.frame(
  id = rep(c("a", "b"), each = 3), time = rep(1:3, 2),
  value = c(10, 30, 55, 12, 20, 28))
timeToMilestone(d, "id", "time", "value", threshold = 50)
#>   id time event attained
#> 1  a    3     1     TRUE
#> 2  b    3     0    FALSE
```
