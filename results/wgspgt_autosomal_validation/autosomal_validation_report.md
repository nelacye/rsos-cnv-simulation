# Autosomal Validation Report

## Scope

This analysis uses `wgspgt_log2like.csv` and restricts validation to autosomes only (`chr1`-`chr22`).

## Empirical Threshold

- Chromosome abnormality threshold: `|log2like| > 0.444`
- Threshold source: 95th percentile of `|log2like|` across autosomes

## Embryo-Level Rule

- An embryo is called `abnormal` if it has at least one abnormal autosome.
- Embryo-level summaries include number of abnormal autosomes, mean log2-like value, and variance.

## EmbryoStatus Comparison

- Binary truth mapping used for metrics:
  - `unaffected` -> `normal`
  - all other statuses -> `abnormal`
- This is a pragmatic first-pass comparison because no dedicated diagnosis sheet is being used here.

## Performance

- Sensitivity: 0.708
- Specificity: 0.429
- Accuracy: 0.645
- TP / TN / FP / FN: 17 / 3 / 4 / 7

## Interpretation

- These results should be treated as exploratory because the pseudo log2-like transformation is not yet equivalent to original QDNAseq logR.
- The binary truth mapping based on `EmbryoStatus` is also only an approximation of CNV abnormality.
- Even so, this autosomes-only pass is a useful checkpoint for whether the transformed signal carries embryo-level abnormality information.
