# Pseudo log2-like Transformation Report

## Transformation

For each embryo:

1. Compute the median autosomal value across `chr1`-`chr22`.
2. Divide each chromosome value by that embryo-specific autosomal median.
3. Take `log2` of the normalized value.

## Distribution Check

- Overall median transformed value: 0.002
- Overall interquartile range: -0.071 to 0.127
- Autosome group median: -0.000
- Autosome median absolute deviation from zero (using absolute transformed values): 0.087
- `chrAut` median transformed value: -0.001
- `chrX` median transformed value: 3.020

## Interpretation

- Autosomes are centered very close to 0 after transformation, which is consistent with a log2-ratio-like baseline.
- `chrAut` remains close to the autosomal baseline after transformation.
- `chrX` remains far from a conventional CNV-like range after transformation.
- This normalization improves comparability across embryos because it removes embryo-specific scale offsets.
- However, the transformed values should still be treated as pseudo log2-like rather than validated QDNAseq logR values.

## Practical Conclusion

- The transformation partially improves the scale, but it does not fully establish a trustworthy log2-like CNV signal.
- Use `wgspgt_log2like.csv` for exploratory modeling only until the original QDNAseq `LogR` or `SegLogR` outputs can be recovered.
