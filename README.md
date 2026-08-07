# Borderline CNVs in PGT-A: a simulation study

## Overview
This repository contains reproducible simulation code for the paper:

Title: Borderline copy number variations in preimplantation genetic testing: when signal becomes noise
Authors: Anere Cye
Journal: Royal Society Open Science (2026)
DOI: 10.5281/zenodo.19421746

## Repository structure

- R/ - simulation functions
- scripts/ - run_all.R and scenario scripts
- results/ - generated automatically

## How to run

Open R in the project root and run:
source(scripts/run_all.R)

All simulations use set.seed(2025). 5000 replicates per condition.

## WGS-PGT validation

The first-pass real-data workflow uses only one workbook and only two sheets:

- `data/wgspgt/Fig.1d_Fig.1e_Fig.1f_Fig.1g_Fig.2b_Fig5.xlsx`
- `Fig1f_WGS` for chromosome-level input
- `Fig1d_Fig1e` for embryo-level metadata

The loader in `R/load_wgspgt_data.R`:
- treats `EmbryoNumber` as `embryo_id`
- reshapes the WGS sheet from wide chromosome columns to tidy rows
- joins `meanDepth`, `meanBreadth`, `EmbryoStatus`, and `Method`
- skips all other sheets for now

Run the validation analysis from the project root with:
`Rscript analysis_wgspgt_validation.R`

This first pass does not require a diagnosis sheet yet. It saves a cleaned
dataset and descriptive threshold-based outputs to `results/wgspgt_validation/`:
- `wgspgt_cleaned.csv`
- `wgspgt_threshold_calls.csv`
- `wgspgt_chromosome_summary.csv`
- `wgspgt_embryo_summary.csv`
- `wgspgt_method_summary.csv`

## Realistic simulation

The original simulation code is preserved as the idealized model. A new parallel
prototype framework models a more realistic aggregated-signal path:

1. latent embryo state
2. latent chromosome-level signal
3. measurement / distortion / aggregation layer
4. embryo-level decision rule

Core files:

- `MODEL_REDESIGN.md`
- `R/simulate_latent_embryo_state.R`
- `R/simulate_latent_chromosome_signal.R`
- `R/simulate_measurement_layer.R`
- `R/classify_embryo_realistic.R`
- `analysis_realistic_simulation.R`

Run the realistic prototype from the project root with:
`Rscript analysis_realistic_simulation.R`

Outputs are written to `results/realistic_simulation/`, including:

- chromosome-level observed signal tables
- embryo-level summaries
- calibration cutoffs for out-of-sample evaluation
- evaluation metrics for out-of-sample evaluation
- confusion matrices
- sensitivity / specificity / accuracy tables
- decision-rule comparison tables
- idealized vs realistic comparison tables
- signal distribution plots

## Key results (Scenario B, n=5000)

False normal rate M=15%: 34.5% [33.2-35.9%]
False positive rate diploid: 54.4% [53.0-55.8%]
