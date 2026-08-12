# Which optional modelling backends are available

The inference functions load heavy modelling packages only when present
(guarded via `requireNamespace`). This reports availability so callers
can fail early with a helpful message.

## Usage

``` r
clinStatsBackends()
```

## Value

A named logical vector for `mmrm`, `emmeans`, `lme4`, and
`SingleCaseES`.

## Examples

``` r
clinStatsBackends()
#> mmrm() registered as emmeans extension
#>         mmrm      emmeans         lme4 SingleCaseES 
#>         TRUE         TRUE         TRUE         TRUE 
```
