# ============================================================================
# replicate_tables_figures.R
# MESA Inflammatory Burden & HF Disparity Study — Rivera-Mariani Lab
#
# Regenerates all manuscript tables (CSV) and figures (PDF + PNG) from the
# summary CSV files included in this replication package.
#
# Usage:   Rscript replicate_tables_figures.R
# Inputs:  ./data/*.csv        (summary result files; no raw MESA data)
# Outputs: ./tables/*.csv      Table 1–4, Supplementary S1–S4
#          ./figures/*.pdf/png Figure 1–5
#
# Dependencies: ggplot2, dplyr, tidyr, scales (install.packages if missing)
# All statistics are pre-computed per StatisticalSOP_v2; this script only
# formats and renders — it performs no model fitting and requires no raw
# participant-level data (compliant with the MESA Data Use Agreement).
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(scales)
})

dir.create("tables",  showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

D <- function(f) read.csv(file.path("data", f), check.names = FALSE, stringsAsFactors = FALSE)

save_fig <- function(plot, name, w = 8, h = 5.5) {
  ggsave(file.path("figures", paste0(name, ".pdf")), plot, width = w, height = h, device = "pdf")
  ggsave(file.path("figures", paste0(name, ".png")), plot, width = w, height = h, dpi = 300)
  message("Figure written: ", name)
}

race_levels <- c("White", "Chinese American", "Black", "Hispanic")
race_cols   <- c("White" = "#4C72B0", "Chinese American" = "#55A868",
                 "Black" = "#C44E52", "Hispanic" = "#8172B2")
theme_pub <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12),
        legend.position = "bottom")

# ============================================================================
# TABLE 1 — Baseline characteristics & biomarker distributions by race
# ============================================================================
n_eth   <- D("mesa_hf_n_by_ethnicity.csv")
ev_eth  <- D("mesa_hf_events_by_ethnicity.csv")
bio_eth <- D("mesa_hf_biomarker_distribution_by_ethnicity.csv")

t1_top <- ev_eth %>%
  transmute(Characteristic = "n (HF events, %)",
            race_label,
            value = sprintf("%d (%d, %.2f%%)", n, hf_events, hf_event_pct)) %>%
  pivot_wider(names_from = race_label, values_from = value)

t1_bio <- bio_eth %>%
  transmute(Characteristic = paste0(biomarker_label, ", median (IQR) z"),
            race_label, value = median_iqr) %>%
  pivot_wider(names_from = race_label, values_from = value)

table1 <- bind_rows(t1_top, t1_bio) %>%
  select(Characteristic, all_of(race_levels))
write.csv(table1, "tables/Table1_baseline_by_ethnicity.csv", row.names = FALSE)
message("Table 1 written")

# ============================================================================
# TABLE 2 — Primary Cox models: complete-case + MI pooled
# ============================================================================
cc <- D("mesa_hf_additive_biomarker_models_compact.csv") %>%
  rename(`CC HR (95% CI)` = `HR (95% CI)`, `CC Wald P` = `Wald P`)
mi <- D("mesa_hf_mi_pooled_biomarker_models.csv") %>%
  transmute(biomarker = biomarker_label,
            `MI HR (95% CI)` = hr_ci, `MI P` = p_value_fmt)
table2 <- left_join(cc, mi, by = "biomarker")
write.csv(table2, "tables/Table2_primary_cox_models.csv", row.names = FALSE)
message("Table 2 written")

# ============================================================================
# TABLE 3 — Race × biomarker interactions + group-specific HRs
# ============================================================================
inter <- D("mesa_hf_ethnicity_biomarker_interactions.csv") %>%
  transmute(biomarker,
            `Global interaction P` = signif(p_interaction, 3),
            `Holm P` = p_adj_holm, `BH q` = signif(q_adj_bh, 3))
write.csv(inter, "tables/Table3a_interaction_tests.csv", row.names = FALSE)

grp <- D("mesa_hf_biomarker_group_specific_slopes_compact.csv")
write.csv(grp, "tables/Table3b_group_specific_HRs.csv", row.names = FALSE)
message("Table 3 written (3a interactions, 3b group-specific)")

# ============================================================================
# TABLE 4 — Sequential mediation of Black–White gap (Exam 1)
# ============================================================================
med <- D("mesa_hf_black_white_gap_mediation_compact.csv") %>%
  filter(analysis == "Exam 1")
write.csv(med, "tables/Table4_sequential_mediation.csv", row.names = FALSE)
message("Table 4 written")

# ============================================================================
# SUPPLEMENTARY TABLES
# ============================================================================
write.csv(D("mesa_hf_mi_vs_complete_case_pooled.csv") %>%
            select(biomarker_label, cc_hr_ci, cc_p_value_fmt, mi_hr_ci, mi_p_value_fmt),
          "tables/SuppTableS1_cc_vs_mi.csv", row.names = FALSE)
write.csv(D("mesa_hf_incremental_model_performance_compact.csv"),
          "tables/SuppTableS2_incremental_performance.csv", row.names = FALSE)
write.csv(D("mesa_hf_environmental_robustness_models_compact.csv"),
          "tables/SuppTableS3_environmental_robustness.csv", row.names = FALSE)
write.csv(D("mesa_hf_biomarker_missingness_by_ethnicity.csv"),
          "tables/SuppTableS4_missingness_by_ethnicity.csv", row.names = FALSE)
message("Supplementary tables S1–S4 written")

# ============================================================================
# FIGURE 1 — Biomarker distributions by race/ethnicity (median + IQR)
# ============================================================================
f1 <- bio_eth %>%
  mutate(race_label = factor(race_label, levels = race_levels)) %>%
  ggplot(aes(x = biomarker_label, y = median, fill = race_label)) +
  geom_col(position = position_dodge(width = .82), width = .75) +
  geom_errorbar(aes(ymin = q1, ymax = q3),
                position = position_dodge(width = .82), width = .25, linewidth = .4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = race_cols, name = NULL) +
  labs(title = "Figure 1. Standardized biomarker distributions by race/ethnicity",
       subtitle = "Bars = median z-score; error bars = IQR (Exam 1)",
       x = NULL, y = "Standardized value (z)") +
  theme_pub
save_fig(f1, "Figure1_biomarker_distributions", w = 9, h = 5.5)

# ============================================================================
# FIGURE 2 — Forest plot: primary HRs (CC + MI) + composite index
# ============================================================================
parse_ci <- function(s) {
  m <- regmatches(s, gregexpr("[0-9.]+", s))[[1]]
  as.numeric(m[1:3])
}
cc_raw <- D("mesa_hf_additive_biomarker_models.csv") %>%
  transmute(label = biomarker_label, hr, lcl, ucl, set = "Complete-case")
mi_raw <- D("mesa_hf_mi_pooled_biomarker_models.csv") %>%
  transmute(label = biomarker_label, hr, lcl, ucl, set = "Multiple imputation")
comp <- D("mesa_hf_composite_inflammatory_burden_additive.csv") %>%
  filter(analysis == "Exam 1") %>%
  transmute(label = "Composite burden (PC1)", hr, lcl, ucl, set = "Complete-case")

forest <- bind_rows(cc_raw, mi_raw, comp) %>%
  mutate(label = factor(label, levels = rev(c("IL-6 (log-z)", "D-dimer (log-z)",
                                              "CRP (log-z)", "Fibrinogen (z)",
                                              "Composite burden (PC1)"))))
f2 <- ggplot(forest, aes(x = hr, y = label, color = set, shape = set)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_pointrange(aes(xmin = lcl, xmax = ucl),
                  position = position_dodge(width = .55), size = .55) +
  scale_color_manual(values = c("Complete-case" = "#1F3864",
                                "Multiple imputation" = "#2E75B6"), name = NULL) +
  scale_shape_manual(values = c(16, 1), name = NULL) +
  scale_x_continuous(limits = c(0.9, 1.5)) +
  labs(title = "Figure 2. Adjusted hazard ratios for incident HF per 1-SD biomarker",
       x = "Hazard ratio (95% CI)", y = NULL) +
  theme_pub
save_fig(f2, "Figure2_forest_primary_HRs", w = 8, h = 4.8)

# ============================================================================
# FIGURE 3 — PCA loadings (Exam 1 composite)
# ============================================================================
load1 <- D("mesa_hf_composite_inflammatory_burden_loadings.csv") %>%
  filter(analysis == "Exam 1")
f3 <- ggplot(load1, aes(x = reorder(marker_label, loading_pc1), y = loading_pc1)) +
  geom_col(fill = "#0D6B8A", width = .6) +
  geom_text(aes(label = sprintf("%.3f", loading_pc1)), hjust = -0.15, size = 3.4) +
  coord_flip() +
  scale_y_continuous(limits = c(0, .7)) +
  labs(title = "Figure 3. PC1 loadings — composite inflammatory burden index",
       subtitle = sprintf("Variance explained = %.1f%% (n = %s complete)",
                          100 * load1$variance_explained_pc1[1],
                          format(load1$n_complete[1], big.mark = ",")),
       x = NULL, y = "PC1 loading") +
  theme_pub
save_fig(f3, "Figure3_PCA_loadings", w = 7, h = 4.2)

# ============================================================================
# FIGURE 4 — Group-specific HRs (4 biomarkers × 4 race groups)
# ============================================================================
grp_full <- D("mesa_hf_biomarker_group_specific_slopes.csv") %>%
  mutate(race_label = factor(race_label, levels = race_levels))
f4 <- ggplot(grp_full, aes(x = hr, y = race_label, color = race_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_pointrange(aes(xmin = lcl, xmax = ucl), size = .45) +
  facet_wrap(~ biomarker_label, ncol = 2) +
  scale_color_manual(values = race_cols, guide = "none") +
  coord_cartesian(xlim = c(0.6, 2.0)) +
  labs(title = "Figure 4. Group-specific hazard ratios by race/ethnicity",
       subtitle = "All race × biomarker multiplicative interactions null (Holm-corrected P = 1.0)",
       x = "Hazard ratio (95% CI) per 1-SD", y = NULL) +
  theme_pub
save_fig(f4, "Figure4_group_specific_HRs", w = 9, h = 6.5)

# ============================================================================
# FIGURE 5 — Sequential mediation waterfall (Black–White gap, Exam 1)
# ============================================================================
med_f <- med %>%
  mutate(model = factor(model, levels = model),
         pct = `% excess risk explained vs clinical-only`,
         hr_num = as.numeric(sub(" .*", "", `Black vs White HR (95% CI)`)))
f5 <- ggplot(med_f, aes(x = model, y = abs(pct))) +
  geom_col(fill = "#C44E52", width = .62) +
  geom_text(aes(label = sprintf("%.1f%%", abs(pct))), vjust = -0.4, size = 3.4) +
  scale_y_continuous(limits = c(0, 70)) +
  labs(title = "Figure 5. Sequential mediation of the Black–White HF risk gap",
       subtitle = "Percent excess risk explained vs. clinical-only model (Exam 1 biomarkers)",
       x = NULL, y = "% excess risk explained") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 28, hjust = 1))
save_fig(f5, "Figure5_sequential_mediation", w = 9, h = 5.5)

message("\nAll tables and figures regenerated successfully.")
message("Tables:  ./tables/   |  Figures: ./figures/")
