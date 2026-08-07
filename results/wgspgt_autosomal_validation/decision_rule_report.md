# Decision Rule Comparison Report

## Scope

This comparison uses `wgspgt_log2like.csv`, autosomes only, and strict mode only.

## Shared Inputs

- Chromosome abnormality cutoff for count-based rules: `|log2like| > 0.444`
- Mean absolute log2like threshold: `0.406`
- Sum absolute log2like threshold: `8.923`
- Mean and sum thresholds were set from the 95th percentile of the clear-normal embryos in strict mode.

## Rules Tested

- >=1 abnormal chromosome
- >=2 abnormal chromosomes
- >=3 abnormal chromosomes
- mean absolute log2like > threshold
- sum absolute log2like > threshold

## Baseline

- Baseline rule (`>=1 abnormal chromosome`) sensitivity: 0.632
- Baseline rule (`>=1 abnormal chromosome`) specificity: 0.429
- Baseline rule (`>=1 abnormal chromosome`) accuracy: 0.577

## Candidate Improvements

- No tested rule improved specificity without a substantial sensitivity drop relative to the baseline.

## Interpretation

- Rules with higher chromosome-count requirements usually trade sensitivity for specificity.
- Aggregate signal rules based on embryo-wide mean or sum absolute log2like can be more conservative than the `>=1 chromosome` rule.
- Use `decision_rule_comparison.csv` for the ranked metrics and `decision_rule_confusion_matrices.csv` for the underlying confusion tables.
