# Mesa_Project

Reproducible workflow for MESA heart-failure biomarker analyses, with a public **aggregate-only replication layer** and a local **controlled-access analysis layer**.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18987368.svg)](https://doi.org/10.5281/zenodo.18987368)

## Repository layout

- `analysis/` — controlled-access analysis scripts (run locally with BioLINCC-approved data access)
  - `mesa_hf_ethnicity_biomarker.R` (canonical analysis script)
  - `run_main_pipeline.R` (entrypoint)
  - `generate_table_s5_minimal_models.R` (minimally adjusted sensitivity table export)
  - `generate_posthoc_interaction_power.R` (post hoc interaction-power calculation export)
- `replication/` — public replication package using **precomputed aggregate CSVs only**
  - `replicate_tables_figures.R` / `.py`
  - `data/` (aggregate summary inputs)
  - `tables/`, `figures/` (rendered outputs)
- `archive/legacy_scripts/` — historical pipeline versions (`v1`–`v6`) retained for provenance
- `data/rawdata/` — local-only controlled-access files (ignored by Git)

## Current status (2026-09)

The project now supports manuscript-aligned replication and supplemental updates, including:

- minimally adjusted sensitivity results for supplemental Table S5
- post hoc interaction-power calculation artifacts
- updated manuscript-facing aggregate outputs in replication data/tables

## How to run

### Controlled-access analysis (local only)

Run from repo root in R:

```r
source("analysis/run_main_pipeline.R")
```

This requires approved MESA/BioLINCC data files in `data/rawdata/` and must not be pushed publicly.

### Public replication (aggregate only)

From `replication/`:

```bash
Rscript replicate_tables_figures.R
```

This regenerates manuscript tables/figures from precomputed aggregate CSV inputs.

## Data governance

This repository must not publish participant-level MESA data or identifiable records. Only aggregate, non-identifiable outputs are intended for public sharing.

## Citation

If you use this workflow, cite the Zenodo record and the MESA cohort reference in your manuscript.
