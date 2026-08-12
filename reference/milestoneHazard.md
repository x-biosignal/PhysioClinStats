# Hazard ratio for attaining a milestone

Fits a Cox model for the time-to-milestone endpoint against a grouping
variable and returns the group hazard ratio - the relative rate of
attaining the milestone - with its confidence interval.

## Usage

``` r
milestoneHazard(milestone, group, conf_level = 0.95)
```

## Arguments

- milestone:

  A time-to-milestone data.frame from
  [`timeToMilestone()`](https://x-biosignal.github.io/PhysioClinStats/reference/timeToMilestone.md)
  (columns `time`, `event`), joined with the grouping variable.

- group:

  Column name of the grouping variable.

- conf_level:

  Confidence level (default 0.95).

## Value

An `AnalysisResult` (`type = "milestone_hazard"`) as
[`coxModel()`](https://x-biosignal.github.io/PhysioClinStats/reference/coxModel.md).

## See also

[`timeToMilestone()`](https://x-biosignal.github.io/PhysioClinStats/reference/timeToMilestone.md),
[`coxModel()`](https://x-biosignal.github.io/PhysioClinStats/reference/coxModel.md)
