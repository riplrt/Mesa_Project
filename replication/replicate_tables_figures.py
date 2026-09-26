#!/usr/bin/env python3
# ============================================================================
# replicate_tables_figures.py  (Python mirror of replicate_tables_figures.R)
# MESA Inflammatory Burden & HF Disparity Study — Rivera-Mariani Lab
#
# Usage:   python3 replicate_tables_figures.py
# Inputs:  ./data/*.csv   (summary results; no raw MESA participant data)
# Outputs: ./tables/*.csv, ./figures/*.pdf + *.png
# Deps:    pandas, matplotlib  (pip install pandas matplotlib)
# ============================================================================
import os, re
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.makedirs("tables", exist_ok=True)
os.makedirs("figures", exist_ok=True)

def D(f): return pd.read_csv(os.path.join("data", f))

RACE_LEVELS = ["White", "Chinese American", "Black", "Hispanic"]
RACE_COLS = {"White":"#4C72B0","Chinese American":"#55A868","Black":"#C44E52","Hispanic":"#8172B2"}

def save_fig(fig, name, w=8, h=5.5):
    fig.set_size_inches(w, h)
    fig.savefig(f"figures/{name}.pdf", bbox_inches="tight")
    fig.savefig(f"figures/{name}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Figure written:", name)

# ── TABLE 1 ─────────────────────────────────────────────────────────────────
ev  = D("mesa_hf_events_by_ethnicity.csv")
bio = D("mesa_hf_biomarker_distribution_by_ethnicity.csv")

row_top = {"Characteristic": "n (HF events, %)"}
for _, r in ev.iterrows():
    row_top[r["race_label"]] = f'{r["n"]} ({r["hf_events"]}, {r["hf_event_pct"]:.2f}%)'
rows = [row_top]
for lab, sub in bio.groupby("biomarker_label"):
    row = {"Characteristic": f"{lab}, median (IQR) z"}
    for _, r in sub.iterrows():
        row[r["race_label"]] = r["median_iqr"]
    rows.append(row)
pd.DataFrame(rows)[["Characteristic"] + RACE_LEVELS].to_csv(
    "tables/Table1_baseline_by_ethnicity.csv", index=False)
print("Table 1 written")

# ── TABLE 2 ─────────────────────────────────────────────────────────────────
cc = D("mesa_hf_additive_biomarker_models_compact.csv").rename(
    columns={"HR (95% CI)": "CC HR (95% CI)", "Wald P": "CC Wald P"})
mi = D("mesa_hf_mi_pooled_biomarker_models.csv")[["biomarker_label","hr_ci","p_value_fmt"]]
mi.columns = ["biomarker", "MI HR (95% CI)", "MI P"]
cc.merge(mi, on="biomarker", how="left").to_csv("tables/Table2_primary_cox_models.csv", index=False)
print("Table 2 written")

# ── TABLE 3 ─────────────────────────────────────────────────────────────────
inter = D("mesa_hf_ethnicity_biomarker_interactions.csv")
inter_out = inter[["biomarker","p_interaction","p_adj_holm","q_adj_bh"]].copy()
inter_out.columns = ["biomarker","Global interaction P","Holm P","BH q"]
inter_out["Global interaction P"] = inter_out["Global interaction P"].round(3)
inter_out["BH q"] = inter_out["BH q"].round(3)
inter_out.to_csv("tables/Table3a_interaction_tests.csv", index=False)
D("mesa_hf_biomarker_group_specific_slopes_compact.csv").to_csv(
    "tables/Table3b_group_specific_HRs.csv", index=False)
print("Table 3 written (3a, 3b)")

# ── TABLE 4 ─────────────────────────────────────────────────────────────────
med = D("mesa_hf_black_white_gap_mediation_compact.csv")
med_e1 = med[med["analysis"] == "Exam 1"].copy()
med_e1.to_csv("tables/Table4_sequential_mediation.csv", index=False)
print("Table 4 written")

# ── SUPPLEMENTARY S1–S4 ─────────────────────────────────────────────────────
D("mesa_hf_mi_vs_complete_case_pooled.csv")[
    ["biomarker_label","cc_hr_ci","cc_p_value_fmt","mi_hr_ci","mi_p_value_fmt"]
].to_csv("tables/SuppTableS1_cc_vs_mi.csv", index=False)
D("mesa_hf_incremental_model_performance_compact.csv").to_csv(
    "tables/SuppTableS2_incremental_performance.csv", index=False)
D("mesa_hf_environmental_robustness_models_compact.csv").to_csv(
    "tables/SuppTableS3_environmental_robustness.csv", index=False)
D("mesa_hf_biomarker_missingness_by_ethnicity.csv").to_csv(
    "tables/SuppTableS4_missingness_by_ethnicity.csv", index=False)
print("Supplementary tables S1–S4 written")

# ── FIGURE 1 ────────────────────────────────────────────────────────────────
fig, ax = plt.subplots()
markers = sorted(bio["biomarker_label"].unique())
x = range(len(markers)); bw = 0.19
for i, race in enumerate(RACE_LEVELS):
    sub = bio[bio["race_label"] == race].set_index("biomarker_label").loc[markers]
    pos = [xi + (i - 1.5) * bw for xi in x]
    ax.bar(pos, sub["median"], width=bw, color=RACE_COLS[race], label=race)
    ax.errorbar(pos, sub["median"],
                yerr=[sub["median"] - sub["q1"], sub["q3"] - sub["median"]],
                fmt="none", ecolor="black", elinewidth=0.7, capsize=2)
ax.axhline(0, ls="--", c="grey", lw=0.8)
ax.set_xticks(list(x)); ax.set_xticklabels(markers, fontsize=9)
ax.set_ylabel("Standardized value (z)")
ax.set_title("Figure 1. Standardized biomarker distributions by race/ethnicity\n(bars = median z; error bars = IQR, Exam 1)", fontsize=11)
ax.legend(fontsize=8, ncol=4, loc="upper left")
save_fig(fig, "Figure1_biomarker_distributions", 9, 5.5)

# ── FIGURE 2 ────────────────────────────────────────────────────────────────
cc_raw = D("mesa_hf_additive_biomarker_models.csv")[["biomarker_label","hr","lcl","ucl"]]
cc_raw["set"] = "Complete-case"
mi_raw = D("mesa_hf_mi_pooled_biomarker_models.csv")[["biomarker_label","hr","lcl","ucl"]]
mi_raw["set"] = "Multiple imputation"
comp = D("mesa_hf_composite_inflammatory_burden_additive.csv")
comp = comp[comp["analysis"] == "Exam 1"][["hr","lcl","ucl"]]
comp["biomarker_label"] = "Composite burden (PC1)"; comp["set"] = "Complete-case"
forest = pd.concat([cc_raw, mi_raw, comp], ignore_index=True)
order = ["Composite burden (PC1)","Fibrinogen (z)","CRP (log-z)","D-dimer (log-z)","IL-6 (log-z)"]
ymap = {l: i for i, l in enumerate(order)}
fig, ax = plt.subplots()
for s, off, col, mk in [("Complete-case", -0.13, "#1F3864", "o"),
                        ("Multiple imputation", 0.13, "#2E75B6", "o")]:
    sub = forest[forest["set"] == s]
    ys = [ymap[l] + off for l in sub["biomarker_label"]]
    ax.errorbar(sub["hr"], ys, xerr=[sub["hr"]-sub["lcl"], sub["ucl"]-sub["hr"]],
                fmt=mk, color=col, label=s, capsize=3,
                markerfacecolor=col if s == "Complete-case" else "white")
ax.axvline(1, ls="--", c="grey", lw=0.8)
ax.set_yticks(range(len(order))); ax.set_yticklabels(order, fontsize=9)
ax.set_xlabel("Hazard ratio (95% CI) per 1-SD")
ax.set_title("Figure 2. Adjusted hazard ratios for incident heart failure", fontsize=11)
ax.set_xlim(0.9, 1.5); ax.legend(fontsize=8)
save_fig(fig, "Figure2_forest_primary_HRs", 8, 4.8)

# ── FIGURE 3 ────────────────────────────────────────────────────────────────
load1 = D("mesa_hf_composite_inflammatory_burden_loadings.csv")
load1 = load1[load1["analysis"] == "Exam 1"].sort_values("loading_pc1")
fig, ax = plt.subplots()
ax.barh(load1["marker_label"], load1["loading_pc1"], color="#0D6B8A", height=0.55)
for y, v in enumerate(load1["loading_pc1"]):
    ax.text(v + 0.008, y, f"{v:.3f}", va="center", fontsize=9)
ax.set_xlim(0, 0.7); ax.set_xlabel("PC1 loading")
ve = load1["variance_explained_pc1"].iloc[0] * 100
nc = int(load1["n_complete"].iloc[0])
ax.set_title(f"Figure 3. PC1 loadings — composite inflammatory burden index\nVariance explained = {ve:.1f}% (n = {nc:,} complete)", fontsize=11)
save_fig(fig, "Figure3_PCA_loadings", 7, 4.2)

# ── FIGURE 4 ────────────────────────────────────────────────────────────────
grp = D("mesa_hf_biomarker_group_specific_slopes.csv")
bio_labels = sorted(grp["biomarker_label"].unique())
fig, axes = plt.subplots(2, 2, sharex=True)
for ax, lab in zip(axes.flat, bio_labels):
    sub = grp[grp["biomarker_label"] == lab]
    for i, race in enumerate(RACE_LEVELS):
        r = sub[sub["race_label"] == race].iloc[0]
        ax.errorbar(r["hr"], i, xerr=[[r["hr"]-r["lcl"]],[r["ucl"]-r["hr"]]],
                    fmt="o", color=RACE_COLS[race], capsize=3)
    ax.axvline(1, ls="--", c="grey", lw=0.8)
    ax.set_yticks(range(4)); ax.set_yticklabels(RACE_LEVELS, fontsize=8)
    ax.set_title(lab, fontsize=10); ax.set_xlim(0.6, 2.0)
fig.suptitle("Figure 4. Group-specific hazard ratios by race/ethnicity\nAll multiplicative interactions null (Holm-corrected P = 1.0)", fontsize=11)
fig.supxlabel("Hazard ratio (95% CI) per 1-SD", fontsize=10)
fig.tight_layout(rect=[0, 0.02, 1, 0.93])
save_fig(fig, "Figure4_group_specific_HRs", 9, 6.5)

# ── FIGURE 5 ────────────────────────────────────────────────────────────────
pct_col = "% excess risk explained vs clinical-only"
med_f = med_e1.copy()
med_f["pct_abs"] = med_f[pct_col].abs()
fig, ax = plt.subplots()
ax.bar(range(len(med_f)), med_f["pct_abs"], color="#C44E52", width=0.6)
for i, v in enumerate(med_f["pct_abs"]):
    ax.text(i, v + 1.2, f"{v:.1f}%", ha="center", fontsize=9)
ax.set_xticks(range(len(med_f)))
ax.set_xticklabels(med_f["model"], rotation=28, ha="right", fontsize=8)
ax.set_ylim(0, 70); ax.set_ylabel("% excess risk explained")
ax.set_title("Figure 5. Sequential mediation of the Black–White HF risk gap\n(vs. clinical-only model, Exam 1 biomarkers)", fontsize=11)
save_fig(fig, "Figure5_sequential_mediation", 9, 5.5)

print("\nAll tables and figures regenerated successfully.")
print("Tables: ./tables/  |  Figures: ./figures/")
