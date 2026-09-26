#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dplyr))
src <- '../data/derived/mesa_hf_additive_biomarker_model_sequence_compact.csv'
out_data <- '../replication_reproduced_2026-08-14/data/mesa_hf_minimally_adjusted_biomarker_models_compact.csv'
out_table <- '../replication_reproduced_2026-08-14/tables/SuppTableS5_minimally_adjusted_models.csv'
x <- read.csv(src, check.names = FALSE)
s5 <- x |>
  filter(model == 'Model 1') |>
  select(biomarker, n, events, `HR (95% CI)`, `Wald P`) |>
  mutate(biomarker = factor(biomarker, levels = c('IL-6 (log-z)', 'D-dimer (log-z)', 'CRP (log-z)', 'Fibrinogen (z)'))) |>
  arrange(biomarker) |>
  mutate(biomarker = as.character(biomarker))
write.csv(s5, out_data, row.names = FALSE)
write.csv(s5, out_table, row.names = FALSE)
message('Wrote: ', out_data)
message('Wrote: ', out_table)
