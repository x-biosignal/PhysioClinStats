# Reshape a longitudinal PhysioExperiment to long format

Flattens an epoched or longitudinal `PhysioExperiment` into a long
`data.frame(subject, time, channel, value)` suitable for
[`fitMixedModel()`](https://x-biosignal.github.io/PhysioClinStats/reference/fitMixedModel.md),
[`fitMMRM()`](https://x-biosignal.github.io/PhysioClinStats/reference/fitMMRM.md)
and the marginal-means helpers.

## Usage

``` r
peToLong(pe, assay = NULL, value_name = "value", subjects = NULL)
```

## Arguments

- pe:

  A `PhysioExperiment`. A 2D assay is treated as time x channel for a
  single subject; a 3D assay as time x channel x subject/trial.

- assay:

  Assay name (default the first).

- value_name:

  Name for the value column (default `"value"`).

- subjects:

  Optional character/factor labels for the third dimension
  (subjects/trials); defaults to the `colnames`/index.

## Value

A long `data.frame` with columns `subject`, `time`, `channel` and the
value column.

## See also

[`fitMixedModel()`](https://x-biosignal.github.io/PhysioClinStats/reference/fitMixedModel.md),
[`fitMMRM()`](https://x-biosignal.github.io/PhysioClinStats/reference/fitMMRM.md)
