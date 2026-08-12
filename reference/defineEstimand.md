# Define an ICH E9(R1) estimand

Constructs a validated estimand from its five ICH E9(R1) attributes. The
returned object is a fully-named list (class `"estimand"`) that
round-trips through the `estimand` slot of an
[`AnalysisResult`](https://x-biosignal.r-universe.dev/PhysioCore/reference/AnalysisResult.html).

## Usage

``` r
defineEstimand(
  treatment,
  population,
  endpoint,
  intercurrent = list(event = NA_character_, strategy = "treatment-policy"),
  summary_measure = "difference in means"
)
```

## Arguments

- treatment:

  The treatment condition(s) being compared.

- population:

  The target population.

- endpoint:

  The endpoint / variable of interest.

- intercurrent:

  A named list with `event` (the intercurrent event) and `strategy` (one
  of
  [ESTIMAND_STRATEGIES](https://x-biosignal.github.io/PhysioClinStats/reference/ESTIMAND_STRATEGIES.md)).

- summary_measure:

  The population-level summary (e.g. `"difference in means"`).

## Value

An `"estimand"` object (a named list of the five attributes plus the
strategy's analysis `recipe`).

## References

ICH E9(R1) addendum on estimands and sensitivity analysis (2019).

## See also

[`analyseEstimand()`](https://x-biosignal.github.io/PhysioClinStats/reference/analyseEstimand.md),
[`multipleImputation()`](https://x-biosignal.github.io/PhysioClinStats/reference/multipleImputation.md)

## Examples

``` r
defineEstimand(
  treatment = "rehab protocol A vs B", population = "post-stroke",
  endpoint = "6-month gait speed",
  intercurrent = list(event = "treatment discontinuation",
                      strategy = "treatment-policy"),
  summary_measure = "difference in means")
#> <estimand> ICH E9(R1)
#>   treatment:   rehab protocol A vs B 
#>   population:  post-stroke 
#>   endpoint:    6-month gait speed 
#>   intercurrent event: treatment discontinuation 
#>   strategy:    treatment-policy 
#>   summary:     difference in means 
```
