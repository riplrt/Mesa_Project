# Replication Package — MESA Inflammatory Burden & HF Disparity Study
## Rivera-Mariani Lab | Manuscript Draft v1.0

This package regenerates every table and figure in the manuscript
"Inflammatory Burden, Structural Exposure, and Racial Disparities in Incident
Heart Failure: Evidence from the Multi-Ethnic Study of Atherosclerosis (MESA)."

All statistics were pre-computed under StatisticalSOP_v2. The scripts here
perform NO model fitting and contain NO raw participant-level MESA data —
only aggregate summary results — and are therefore compliant with the MESA
Data Use Agreement for sharing among the authorship team.

## Contents

```
replication/
├── README.md                        ← this file
├── replicate_tables_figures.R       ← primary replication script (R)
├── replicate_tables_figures.py      ← identical pipeline in Python
├── data/                            ← 16 summary CSV inputs
│   ├── mesa_hf_n_by_ethnicity.csv
│   ├── mesa_hf_events_by_ethnicity.csv
│   ├── mesa_hf_biomarker_distribution_by_ethnicity.csv
│   ├── mesa_hf_additive_biomarker_models.csv
│   ├── mesa_hf_additive_biomarker_models_compact.csv
│   ├── mesa_hf_mi_pooled_biomarker_models.csv
│   ├── mesa_hf_ethnicity_biomarker_interactions.csv
│   ├── mesa_hf_biomarker_group_specific_slopes.csv
│   ├── mesa_hf_biomarker_group_specific_slopes_compact.csv
│   ├── mesa_hf_black_white_gap_mediation_compact.csv
│   ├── mesa_hf_mi_vs_complete_case_pooled.csv
│   ├── mesa_hf_incremental_model_performance_compact.csv
│   ├── mesa_hf_environmental_robustness_models_compact.csv
│   ├── mesa_hf_biomarker_missingness_by_ethnicity.csv
│   ├── mesa_hf_composite_inflammatory_burden_additive.csv
│   └── mesa_hf_composite_inflammatory_burden_loadings.csv
├── tables/                          ← generated outputs (pre-run for convenience)
└── figures/                         ← generated outputs (PDF + 300-dpi PNG)
```

## How to run

**R (preferred for the lab pipeline):**
```bash
Rscript replicate_tables_figures.R
# deps: install.packages(c("ggplot2","dplyr","tidyr","scales"))
```

**Python (identical output):**
```bash
python3 replicate_tables_figures.py
# deps: pip install pandas matplotlib
```

Both scripts write to `./tables/` and `./figures/`.

## Output → manuscript mapping

| Output file                                | Manuscript element        |
|--------------------------------------------|---------------------------|
| Table1_baseline_by_ethnicity.csv           | Table 1                   |
| Table2_primary_cox_models.csv              | Table 2                   |
| Table3a_interaction_tests.csv              | Table 3 (panel A)         |
| Table3b_group_specific_HRs.csv             | Table 3 (panel B)         |
| Table4_sequential_mediation.csv            | Table 4                   |
| SuppTableS1_cc_vs_mi.csv                   | Supplementary Table S1    |
| SuppTableS2_incremental_performance.csv    | Supplementary Table S2    |
| SuppTableS3_environmental_robustness.csv   | Supplementary Table S3    |
| SuppTableS4_missingness_by_ethnicity.csv   | Supplementary Table S4    |
| Figure1_biomarker_distributions.pdf/.png   | Figure 1                  |
| Figure2_forest_primary_HRs.pdf/.png        | Figure 2                  |
| Figure3_PCA_loadings.pdf/.png              | Figure 3                  |
| Figure4_group_specific_HRs.pdf/.png        | Figure 4                  |
| Figure5_sequential_mediation.pdf/.png      | Figure 5                  |

## Notes for the team

- Figure styling (fonts, palette) may be adjusted to target-journal specs at
  submission; the underlying values must not change.
- Percent-excess-risk values in the source CSV are stored with negative sign
  (attenuation direction); the scripts plot absolute magnitude, matching the
  manuscript text.
- Table 3 is split into two CSV panels (interaction tests; group-specific HRs)
  to be merged into a single two-panel table at Word formatting stage.
- Questions → PI (Dr. Rivera-Mariani).
