# WS8-11 independent numeric validation

- R: 4.3.3
- PhysioClinStats: 0.2.0
- Seeded datasets: 100 x 500 subjects
- Clone/censor mismatches: 0 (limit 0)
- Baseline weight max relative error: 6.54e-16 (limit 1e-10)
- Longitudinal weight max relative error: 0 (limit 1e-10)
- ITT RD bias: -0.0065; 95% CI coverage: 0.980
- Per-protocol RD bias: 0.0079; 95% CI coverage: 0.950
- Positivity stress surfaced: TRUE
- Future-information mutation detected: TRUE
- Cumulative-product grouping mutation detected: TRUE
- Treatment-coding mutation detected: TRUE

## Gates

- clone_agreement: PASS
- baseline_weight: PASS
- longitudinal_weight: PASS
- itt_bias: PASS
- pp_bias: PASS
- itt_coverage: PASS
- pp_coverage: PASS
- positivity: PASS
- future_information_mutation: PASS
- grouping_mutation: PASS
- coding_mutation: PASS

This simulation checks numerical implementation under its generating mechanism; it does not validate causal identification assumptions in observational rehabilitation data.
