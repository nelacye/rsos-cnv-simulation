# EmbryoStatus Mapping Validation Report

## Unique EmbryoStatus Values

- affected
- carrier
- inconclusive
- unaffected

## Clinical Grouping

- clearly abnormal: `affected`
- clearly normal: `unaffected`
- ambiguous / intermediate: `carrier`, `inconclusive`

## Validation Modes

- Shared chromosome abnormality threshold: `|log2like| > 0.444`
- Threshold source: 95th percentile of autosomal absolute log2-like values

### Strict

- Includes only clear normal vs clear abnormal embryos.
- Excludes ambiguous/intermediate statuses.
- Sensitivity: 0.632
- Specificity: 0.429
- Accuracy: 0.577
- TP / TN / FP / FN: 12 / 3 / 4 / 7

### Inclusive

- Uses the current mapping where non-`unaffected` statuses are treated as abnormal.
- Sensitivity: 0.708
- Specificity: 0.429
- Accuracy: 0.645
- TP / TN / FP / FN: 17 / 3 / 4 / 7

## Comparison

- The strict mode is clinically cleaner because it removes ambiguous statuses from the truth set.
- The inclusive mode keeps more embryos but mixes intermediate categories into the abnormal class.
- Use `validation_metrics_by_mode.csv` and `per_embryo_classification_by_mode.csv` for side-by-side comparison.
