# Mesa_Project

Pipeline + modeling scripts for **MESA (Multi-Ethnic Study of Atherosclerosis)** analyses focused on **incident Heart Failure (HF)** and **race/ethnicity × biomarker** associations. The workflow is designed to run locally (in R) using MESA CSV extracts obtained via controlled-access mechanisms (NHLBI BioLINCC).

---

## What’s in this repo

### R scripts (modeling + cohort assembly)
- `mesa_hf_ethnicity_biomarker_v1.r`  
  Early version of the pipeline. Uses a placeholder absolute path (`BASE_DIR <- "/path/to/mesa_data"`).

- `mesa_hf_ethnicity_biomarker_v2.r` (**recommended**)  
  Updated version of the pipeline. Uses repo-local paths (`BASE_DIR <- "data"`), exports analysis-ready datasets, runs interaction scans, and produces publication-ready figures.

### Reference documentation (PDF)
- `MESA_v2023a.pdf`
- `Data Dictionary - Ancillary Studies.pdf`

### Data directory (not tracked with raw data)
- `data/`  
  Intended location for controlled-access CSVs **on your local machine**.
  - `data/derived/` is created automatically for outputs.

---

## Core logic / conventions

### Join key (important)
When merging MESA tables, **always use `MESAID` as the primary join key**. Both scripts enforce this (they stop if `mesaid` is missing after name cleaning).

### Efficiency
The R pipeline is written using vectorized operations (tidyverse). Avoid explicit loops for participant-level or longitudinal operations.

> If you later scale up to very large extracts, prefer **Python (PySpark)** for ingestion/joins/feature engineering and keep **R** for statistical modeling.

### Statistical approach (high-level)
For each candidate biomarker, the scripts fit Cox proportional hazards models and test **biomarker × race/ethnicity** interaction using a likelihood ratio test (reduced vs. interaction model). Results include:
- HR per 1 SD biomarker (with 95% CI)
- main-effect p-value
- interaction p-value
- multiplicity adjustments (Holm and BH/FDR) for interaction p-values

---

## Input data expectations

The scripts expect MESA datasets as **CSV** files on disk. Filenames are hard-coded in the R scripts (see the `paths <- list(...)` section).

Key inputs used (non-exhaustive):
- Events:
  - Preferred: `mesaevthr2020_drepos_20241120.csv` (overall HF endpoint/time)
  - Optional: `mesaevefthru2015_drepos_20200330.csv` (EF subtype fields if present)
- Exam 1 baseline:
  - `mesae1dres20220813.csv`
  - `mesa_site_drepos_20181106.csv`
- Cardiac biomarkers:
  - `mesaas079_drepos_20151118.csv`
  - `mesaas244_drepos_20161011.csv`
- Environmental/geocode (subset creation):
  - `mesa_airexpos_ds_20231211.csv`
  - `mesaaire1_drepos_20240603.csv` (optional)
  - `mesaas023raceseg_ds_20220111.csv`
- Exam 4 + immune:
  - `mesae4dres06222012.csv`
  - `mesaas042_drepos_20150819.csv`

All column names are cleaned to snake_case via `janitor::clean_names()`.

---

## Outputs

All outputs are written under `data/derived/` (v2) or under `${BASE_DIR}/derived` (v1).

### v1 outputs (RDS + results table)
- `mesa_hf_ethnicity_biomarker_interactions.csv`
- `mesa_hf_biomarker_main.rds`
- `mesa_hf_air_inflammation_subset.rds`
- `mesa_hf_exam4_landmark_ready.rds`

### v2 outputs (CSV + results table + figures)
- Results:
  - `mesa_hf_ethnicity_biomarker_interactions.csv`
- Analysis-ready datasets:
  - `mesa_hf_biomarker_main.csv`
  - `mesa_hf_air_inflammation_subset.csv`
  - `mesa_hf_exam4_landmark_ready.csv`
- Figures (PDF and PNG):
  - `fig1_forest_biomarker_hr.(pdf|png)` — forest plot of HRs per SD
  - `fig2_interaction_pvalues.(pdf|png)` — interaction p-values (−log10 p)
  - `fig3_km_hf_free_survival.(pdf|png)` — Kaplan–Meier HF-free survival by race/ethnicity
  - `fig4_biomarker_dist_by_race.(pdf|png)` — biomarker distributions by race/ethnicity

---

## How to run (in R)

### 1) Put controlled-access CSVs in `data/`
From the repo root, create/populate:

```text
data/
  mesae1dres20220813.csv
  mesaevthr2020_drepos_20241120.csv
  ...
```

### 2) Run the recommended pipeline (v2)
In R:

```r
source("mesa_hf_ethnicity_biomarker_v2.r")
```

### 3) (Optional) Run v1
Edit `BASE_DIR` in `mesa_hf_ethnicity_biomarker_v1.r`, then:

```r
source("mesa_hf_ethnicity_biomarker_v1.r")
```

---

## Notes on missing data and standardization

- Certain biomarkers use MESA missingness indicators (e.g., `crp1m`, `fib1m`) to set values to `NA`.
- Biomarkers are transformed to standardized metrics:
  - Z-score: `make_z()`
  - log + Z-score (with small constant): `make_log_z()`

If you plan extended modeling, consider adding:
- Multiple imputation for covariates (e.g., `mice`)
- Systematic z-score normalization for environmental exposures across exam windows

---

## Data Privacy

MESA/NHLBI data are **controlled-access** and must be handled in accordance with **NHLBI BioLINCC data privacy and data use guidelines** (including any applicable Data Use Agreement).  
This repository does not contain any participant-level MESA data.

Users must:

- obtain datasets independently through BioLINCC

- keep raw and derived data local

- avoid committing any identifiable or participant-level data to GitHub

- Only aggregate, non-identifiable outputs should be shared publicly.

---

## Acknowledgements

This project uses data and/or biospecimen-derived measures obtained through the **National Institutes of Health (NIH)** and the **National Heart, Lung, and Blood Institute (NHLBI)**. We acknowledge and thank **NHLBI’s BioLINCC (Biological Specimen and Data Repository Information Coordination Center)** for providing access to MESA resources and supporting data distribution.

---

## Suggested next steps (optional)

- Add `renv` for reproducible R package management.
- Add a `data/README.md` documenting the exact expected input filenames and provenance.
- Add a Python/PySpark ingestion layer for scalable joins/feature engineering (keeping `MESAID` as the primary key).
- For mixed-effects modeling, default to `lme4`/`nlme` with random intercepts for **MESAID** and **SITE** when appropriate.

---

## License / citation

- All analysis code in this repository is released under the MIT License unless otherwise specified.

The license applies only to the analysis scripts and workflow, not to any datasets referenced in this repository.

Copyright (c) [2026] Félix E. Rivera-Mariani



- If you use this workflow or scripts in a publication, please cite:

Rivera-Mariani FE.
MESA Heart Failure Biomarker–Ethnicity Analysis Pipeline.
GitHub repository. 2026.


- Publications using MESA data should also cite the primary cohort description:

Bild DE, Bluemke DA, Burke GL, et al.
Multi-Ethnic Study of Atherosclerosis: objectives and design.
American Journal of Epidemiology. 2002;156(9):871–881.
https://doi.org/10.1093/aje/kwf113
