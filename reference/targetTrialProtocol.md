# Declare a target-trial protocol

Defines the seven protocol components that must be fixed before
emulating a target trial. `time_zero` is a separate alignment rule and
is never inferred from the observed data. Static strategies are scalar
sustained treatment values. Dynamic strategy functions receive the
current row as `data` and the participant's rows through the current
time as `history`; they must not access later rows.

## Usage

``` r
targetTrialProtocol(
  eligibility,
  treatment_strategies,
  assignment,
  time_zero,
  follow_up,
  outcome,
  causal_contrast,
  analysis_plan,
  protocol_id = NULL,
  version = "1.0.0"
)
```

## Arguments

- eligibility:

  Eligibility rule. Use a function of baseline data or a one-sided
  formula for an executable rule.

- treatment_strategies:

  Named list of at least two static scalar or dynamic function
  strategies.

- assignment, outcome, causal_contrast, analysis_plan:

  Non-empty, serializable protocol declarations.

- time_zero:

  Explicit scalar time value or alignment rule.

- follow_up:

  Non-empty follow-up declaration. A finite positive scalar is used as
  the fixed horizon by
  [`targetTrialEmulate()`](https://x-biosignal.github.io/PhysioClinStats/reference/targetTrialEmulate.md).

- protocol_id:

  Optional stable identifier. When omitted, one is derived from the
  protocol declarations.

- version:

  Semantic protocol version.

## Value

An object of class `target_trial_protocol`.

## Examples

``` r
protocol <- targetTrialProtocol(
  eligibility = function(data) data$eligible,
  treatment_strategies = list(never = 0, always = 1),
  assignment = "cloning at eligibility",
  time_zero = 0,
  follow_up = 12,
  outcome = "binary recovery by week 12",
  causal_contrast = "always versus never",
  analysis_plan = "clone-censor-weight"
)
protocol
#> Target-trial protocol: ttp-971eddbb8645 (version 1.0.0 )
#> Strategies: never [static], always [static] 
#> Time zero: 0 
#> Eligibility: function (data) 
#> data$eligible 
#> Assignment: "cloning at eligibility" 
#> Follow-up: 12 
#> Outcome: "binary recovery by week 12" 
#> Causal contrast: "always versus never" 
#> Analysis plan: "clone-censor-weight" 
#> Portable: yes 
#> Validation: 10 of 10 components valid
```
