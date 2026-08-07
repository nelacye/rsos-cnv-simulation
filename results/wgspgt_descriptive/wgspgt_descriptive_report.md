# WGS-PGT Descriptive Report

## Scope

This report uses the workbook-specific cleaned dataset derived from `Fig1f_WGS` and `Fig1d_Fig1e` only. No Bayesian classifier or threshold classifier was applied in this analysis.

Outputs in this directory:

- `overall_summary.csv`
- `by_chromosome_summary.csv`
- `by_status_method_summary.csv`
- `by_status_chromosome_group_summary.csv`
- `embryo_level_summary.csv`
- `chromosome_scale_comparison.csv`
- plot PNG files

## Dataset Overview

- Chromosome-level rows: 744
- Embryos: 31
- EmbryoStatus levels: affected, carrier, inconclusive, unaffected
- Method levels observed: 1

## Overall Distribution

- Median value: 1.670
- Mean value: 3.443
- Interquartile range: 1.354 to 2.201
- Minimum to maximum: 0.674 to 55.987
- The distribution is strongly right-skewed with a long upper tail.

## By Chromosome

- Autosomes are broadly concentrated around medians of roughly 1.5 to 2.1.
- `chrX` median is 38.194, which is 23.0x the pooled autosome median.
- `chrAut` median is 1.696, which is 1.0x the pooled autosome median.
- `chrX` is clearly on a different scale from autosomes.
- `chrAut` is elevated relative to autosomes but remains much closer to the autosome range than `chrX`.

## By EmbryoStatus and Method

- Only one method (`WGS`) is present, so there is no between-method distribution comparison available in this dataset.
- Distribution differences across `EmbryoStatus` are visible, but they occur on top of the strong chromosome-specific scale effects, especially `chrX`.

## Embryo-Level Summaries

- For each embryo, the report includes:
  - median chromosome value
  - range across chromosomes
  - number of extreme chromosomes
- Extreme chromosomes were defined using embryo-specific Tukey fences: values below Q1 - 1.5*IQR or above Q3 + 1.5*IQR within that embryo.
- The embryo with the most extreme chromosomes in this first pass is embryo `34` with 5 extreme chromosomes and a chromosome-value range of 47.453.

## Interpretation of Fig1f_WGS Scale

- These values do not behave like standard log2 ratios.
  Standard CNV log2 ratios are typically centered near 0 for diploid regions, whereas these values have a pooled median of 1.670 and extend to 55.987.
- The autosomes are also not centered near a simple diploid baseline such as 0 or 2-copy-equivalent log2 space.
- The `chrX` values are dramatically inflated relative to autosomes, which argues against this being a uniformly normalized log2-ratio track.
- The data do show chromosome-structured behavior rather than looking like a pure random QC/error score, so they are unlikely to be a generic sequencing QC metric alone.
- The most plausible interpretation from this workbook alone is that `Fig1f_WGS` behaves more like a copy-number-like or coverage-derived signal on a transformed positive scale, not a standard log2-ratio output.
- Because the scale is positive, shifted, and chromosome-dependent, it should not be fed directly into a Bayesian CNV classifier that expects conventional log2 ratios without a prior scale-reconstruction step.

## Practical Conclusion

- `Fig1f_WGS` is useful for descriptive exploration and embryo/chromosome ranking.
- It is not yet suitable as a drop-in replacement for a standard chromosome-level log2-ratio input.
- Before Bayesian-vs-threshold comparison, the next technical step should be to recover or reconstruct the original QDNAseq `LogR`/`SegLogR` scale, or at minimum derive a validated mapping from this workbook scale to a CNV-relevant measurement.
