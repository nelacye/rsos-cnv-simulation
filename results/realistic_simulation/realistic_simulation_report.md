# Realistic Simulation Report

## Scope

- Simulated embryos: 1500
- Calibration fraction: 0.50
- Empirical chromosome cutoff for fixed-threshold logic: 0.668
- Latent statuses: normal, ambiguous, abnormal
- Autosomes only (`chr1`-`chr22`)
- Aggregate embryo rules are reported as out-of-sample evaluation

## Interpretation

- The realistic model inserts an explicit measurement layer between latent chromosome state and observed pseudo-signal.
- Aggregate embryo cutoffs are estimated on a calibration split using only latent normal embryos.
- Threshold and Bayesian chromosome logic are evaluated on the held-out evaluation split.
- The idealized-vs-realistic comparison now uses a parallel comparable-idealized mode with matched latent chromosome burden, so it isolates the measurement layer more cleanly.
- The output tables can be compared directly through sensitivity, specificity, accuracy, and confusion counts.

## Key Output Files

- `decision_rule_comparison.csv`
- `evaluation_metrics.csv`
- `calibration_cutoffs.csv`
- `confusion_matrices.csv`
- `idealized_vs_realistic_metrics.csv`
- `observed_chromosome_signal.csv`
- `embryo_decision_calls.csv`
