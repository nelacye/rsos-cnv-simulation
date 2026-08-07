# Model Redesign

## Why the Old Model Is Idealized

The original simulation framework is intentionally simple:

1. latent mosaicism `M`
2. deterministic chromosome-level copy-ratio relationship
3. direct conversion to log2 ratio
4. additive Gaussian noise
5. fixed threshold or Bayesian decision rule

This structure is useful for studying decision boundaries, but it is idealized in
two important ways:

- it treats the observed chromosome-level value as a nearly direct readout of the
  latent biological state
- it does not separate the biological signal layer from the measurement and
  aggregation layer

In the current repository, the idealized signal generation is mainly implemented
by:

- `R/simulate_cnv.R`
- `R/simulate_embryo_uniform.R`
- `R/simulate_embryo_single_affected.R`

Threshold classification is implemented by:

- `R/classify_embryo.R`

Embryo-level summarization is implemented by:

- `scripts/run_scenario_single.R`
- `scripts/run_scenario_A.R`
- `scripts/run_scenario_B.R`
- parts of `R/simulation.R`

Bayesian posterior calculation currently appears in the analysis layer rather
than a reusable module, especially in:

- `fig_same_ranking_different_decisions.R`
- `R/viotti_comparison.R`

## Why Real-Data Exploration Suggests a Measurement Layer Is Needed

The workbook-derived WGSPGT chromosome values did not behave like native
QDNAseq `LogR` or `SegLogR` outputs.

Key empirical observations from the exploration already performed in this
repository:

- the raw workbook chromosome values were not centered like ordinary log2 ratios
- autosomal median normalization produced a pseudo log2-like signal, but
  discriminative performance remained limited
- autosome-only thresholding preserved sensitivity better than specificity
- changing embryo-level decision rules alone did not solve the core problem

This suggests the main limitation is not only the decision rule. The deeper issue
is that the observed chromosome-level quantity is a distorted or aggregated
measurement of the latent chromosome state.

That motivates an explicit measurement layer in the new simulation.

## New Model Architecture

The new realistic aggregated-signal model is designed as a parallel framework
that keeps the idealized simulation intact.

### Layer 1: Latent Embryo State

Embryos are assigned to one of three latent categories:

- `normal`
- `ambiguous`
- `abnormal`

These states govern:

- how many chromosomes are affected
- how large the latent chromosome effects are
- how often embryo-level false positive and false negative behavior can emerge

### Layer 2: Latent Chromosome Signal

For autosomes `chr1`-`chr22`, each embryo receives a latent chromosome-level
signal in a log2-like biological space.

The minimal prototype models:

- no affected chromosomes for most `normal` embryos
- smaller and fewer affected chromosomes for `ambiguous` embryos
- larger and more numerous affected chromosomes for `abnormal` embryos

### Layer 3: Measurement / Aggregation Layer

The realistic model then distorts the latent chromosome signal into an observed
chromosome-level pseudo-signal.

The prototype includes:

- embryo-specific multiplicative baseline scaling
- embryo-specific additive baseline shift
- chromosome-specific multiplicative distortion
- coverage-dependent noise variance
- per-embryo autosomal median normalization
- observed pseudo log2-like signal after aggregation

This layer is meant to reproduce the empirical pattern seen in the workbook-based
exploration:

- autosomes roughly centered after within-embryo normalization
- moderate sensitivity
- weak specificity under simple embryo-level rules

### Layer 4: Embryo-Level Decision Rule

The prototype supports two decision logics:

- fixed threshold logic on observed pseudo-signal
- Bayesian posterior logic on observed pseudo-signal

The embryo-level rules compared are:

- `>=1` abnormal chromosome
- `>=2` abnormal chromosomes
- `>=3` abnormal chromosomes
- mean aggregate evidence
- sum aggregate evidence

## Comparable Outputs Between Old and New Models

The new realistic framework is intended to produce outputs comparable to the old
idealized framework, even though the internal signal path is different.

Comparable outputs include:

- chromosome-level observed signal distributions
- embryo-level calls
- confusion matrices
- sensitivity / specificity / accuracy
- decision-rule comparison tables

In the minimal prototype, the idealized comparison is represented by a simple
reuse of the old scenario-style logic, while the realistic comparison comes from
the new measurement-aware pipeline.

## Minimal Prototype Assumptions

The new framework intentionally starts small.

Assumptions in the prototype:

- autosomes only
- latent embryo state is known in simulation
- ambiguous embryos are biologically intermediate rather than explicitly mosaic
  at the cell level
- the Bayesian layer is a pragmatic mixture model on the observed pseudo-signal,
  not a full reconstruction of QDNAseq internals
- the empirical goal is to recover the broad failure mode seen in the workbook
  exploration, not to exactly replicate the public WGSPGT data source

These assumptions are explicit in code comments and can be relaxed later.
