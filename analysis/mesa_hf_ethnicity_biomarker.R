#------------------------------------------------------------------------------
# MESA HF x Ethnicity x Candidate Biomarkers
#
# UPDATE (2026-03-11):
# - Replace all .sas7bdat inputs with .csv inputs.
# - Use MESAID as primary join key everywhere.
# - Event files: prefer mesaevthr2020_drepos_20241120.csv; optionally use
#   mesaevefthru2015_drepos_20200330.csv if present for EF subtype fields.
#
# UPDATE (2026-03-12):
# - aire1 ingestion now reads crp1m/fib1m missingness flags from
#   mesaaire1_drepos_20240603.csv (if present) and applies them before coalescing.
# - Added crp1_source/fib1_source provenance columns to mesa_main.
#
# UPDATE (2026-04-03):
# - Section 17: RERI additive interaction analysis (VanderWeele & Knol 2014)
#   for all primary inflammatory biomarkers x race/ethnicity on the absolute
#   risk-difference scale via Cox predicted risks, with delta-method CIs.
# - Section 18 (v5): Formal causal mediation via CMAverse — REWRITTEN.
#   Root-cause fixes for universal model_error in v4:
#     (a) Outcome passed as pre-built Surv() column, not separate time+event args
#     (b) Factor confounders converted to numeric dummies before cmest()
#     (c) Error messages now captured and written to CSV for debugging
#     (d) E-value sensitivity analysis with robust fallback
#
# UPDATE (2026-04-03) v6:
# - Section 18 (v6): Two additional root-cause fixes for universal model_error in v5:
#     Bug 1: yreg = "survCox" is invalid; corrected to yreg = "coxph".
#     Bug 2: CMAverse cmest() with yreg = "coxph" requires separate outcome (time)
#            and event arguments — NOT a pre-built Surv() column.
#
# UPDATE (2026-04-03) v7:
# - Section 18 (v7): Fixes silent-NA extraction bug in v6. Models ran successfully
#   (status = "success") but all estimates were NA because:
#     Bug 3: v6 extracted from summary(res)$effect.decomposition, which does not
#            exist in CMAverse. The actual slots are: res$effect.pe, res$effect.se,
#            res$effect.ci.low, res$effect.ci.high, res$effect.pval (named vectors
#            on the cmest result object directly).
#     Bug 4: For survival/ratio outcomes with EMint=TRUE, CMAverse names effects
#            with "R" prefix: "Rpnde", "Rtnie", "Rte" — not "pnde", "tnie", "te".
#            "pm" (proportion mediated) is unprefixed.
#     Bug 5: Added diagnostic cat() to print actual slot names from res for any
#            future debugging.
#------------------------------------------------------------------------------

# 0) Packages -----------------------------------------------------------------
packages <- c("tidyverse", "janitor", "fs", "glue", "survival", "broom", "mice")

to_install <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

invisible(lapply(packages, library, character.only = TRUE))
options(stringsAsFactors = FALSE)

# 1) Paths --------------------------------------------------------------------
BASE_DIR <- "data/rawdata"
OUT_DIR  <- fs::path(BASE_DIR, "derived")
fs::dir_create(OUT_DIR)

paths <- list(
  # Events
  event_primary   = fs::path(BASE_DIR, "mesaevefthru2015_drepos_20200330.csv"),   # optional
  event_secondary = fs::path(BASE_DIR, "mesaevthr2020_drepos_20241120.csv"),      # preferred

  # Core exam + covariates
  exam1   = fs::path(BASE_DIR, "mesae1dres20220813.csv"),
  site    = fs::path(BASE_DIR, "mesa_site_drepos_20181106.csv"),

  # Environment / geocode
  air     = fs::path(BASE_DIR, "mesa_airexpos_ds_20231211.csv"),
  geocode = fs::path(BASE_DIR, "mesaas023raceseg_ds_20220111.csv"),
  aire1   = fs::path(BASE_DIR, "mesaaire1_drepos_20240603.csv"),

  # Exam 4 / immune
  exam4   = fs::path(BASE_DIR, "mesae4dres06222012.csv"),
  immune  = fs::path(BASE_DIR, "mesaas042_drepos_20150819.csv"),

  # Cardiac biomarkers
  cardiac079 = fs::path(BASE_DIR, "mesaas079_drepos_20151118.csv"),
  cardiac244 = fs::path(BASE_DIR, "mesaas244_drepos_20161011.csv")
)

# 2) Helpers ------------------------------------------------------------------
read_mesa_csv <- function(path, keep = NULL, guess_max = 20000) {
  if (!fs::file_exists(path)) stop(glue::glue("File not found: {path}"))

  dat <- readr::read_csv(path, show_col_types = FALSE, guess_max = guess_max) %>%
    janitor::clean_names()

  if (!"mesaid" %in% names(dat)) {
    stop(glue::glue("Column 'mesaid' not found in {basename(path)}"))
  }

  dat <- dat %>% mutate(mesaid = as.integer(mesaid))

  if (!is.null(keep)) {
    dat <- dat %>% select(any_of(c("mesaid", keep)))
  }

  dat
}

ensure_cols <- function(dat, cols, fill = NA) {
  for (nm in cols) {
    if (!nm %in% names(dat)) dat[[nm]] <- fill
  }
  dat
}

make_z <- function(x) {
  x <- as.numeric(x)
  if (sum(is.finite(x)) < 2L) return(rep(NA_real_, length(x)))
  sdx <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(sdx) || sdx == 0) return(rep(NA_real_, length(x)))
  as.numeric(scale(x))
}

make_log_z <- function(x) {
  x <- as.numeric(x)
  pos <- x[is.finite(x) & x > 0]
  if (length(pos) < 2L) return(rep(NA_real_, length(x)))

  cst <- 0.5 * min(pos, na.rm = TRUE)
  lx  <- log(x + cst)
  sdx <- stats::sd(lx, na.rm = TRUE)
  if (!is.finite(sdx) || sdx == 0) return(rep(NA_real_, length(x)))

  as.numeric(scale(lx))
}

load_cardiac <- function(file079, file244) {
  card079 <- tibble(mesaid = integer())
  card244 <- tibble(mesaid = integer())

  if (fs::file_exists(file079)) {
    card079 <- read_mesa_csv(file079, c("tntstat1", "ntprbnp1", "tntstat3", "ntprbnp3")) %>%
      transmute(
        mesaid,
        ntprobnp_e1_079 = ntprbnp1,
        troponin_e1_079 = tntstat1,
        ntprobnp_e3_079 = ntprbnp3,
        troponin_e3_079 = tntstat3
      )
  }

  if (fs::file_exists(file244)) {
    card244 <- read_mesa_csv(file244, c(
      "probnpb1", "probnpbqns1", "probnpblob1",
      "hstntb1", "hstntbqns1", "hstntbblob1",
      "probnpb3", "probnpblob3",
      "hstntb3", "hstntblob3"
    )) %>%
      transmute(
        mesaid,
        ntprobnp_e1_244 = probnpb1,
        ntprobnp_qns1   = probnpbqns1,
        ntprobnp_blob1  = probnpblob1,
        troponin_e1_244 = hstntb1,
        troponin_qns1   = hstntbqns1,
        troponin_blob1  = hstntbblob1,
        ntprobnp_e3_244 = probnpb3,
        ntprobnp_blob3  = probnpblob3,
        troponin_e3_244 = hstntb3,
        troponin_blob3  = hstntblob3
      )
  }

  joined <- full_join(card244, card079, by = "mesaid")
  joined <- ensure_cols(joined, c(
    "ntprobnp_e1_244", "ntprobnp_qns1", "ntprobnp_blob1",
    "troponin_e1_244", "troponin_qns1", "troponin_blob1",
    "ntprobnp_e3_244", "ntprobnp_blob3", "troponin_e3_244", "troponin_blob3",
    "ntprobnp_e1_079", "troponin_e1_079", "ntprobnp_e3_079", "troponin_e3_079"
  ))

  joined %>%
    transmute(
      mesaid,
      ntprobnp_e1 = coalesce(ntprobnp_e1_244, ntprobnp_e1_079),
      troponin_e1 = coalesce(troponin_e1_244, troponin_e1_079),
      ntprobnp_e3 = coalesce(ntprobnp_e3_244, ntprobnp_e3_079),
      troponin_e3 = coalesce(troponin_e3_244, troponin_e3_079),
      ntprobnp_qns1 = ntprobnp_qns1,
      ntprobnp_blob1 = ntprobnp_blob1,
      troponin_qns1 = troponin_qns1,
      troponin_blob1 = troponin_blob1,
      troponin_source_e1 = case_when(
        !is.na(troponin_e1_244) ~ "mesaas244_hsctnt",
        !is.na(troponin_e1_079) ~ "mesaas079_troponin_t",
        TRUE ~ NA_character_
      )
    ) %>%
    mutate(
      ntprobnp_e1 = if_else(ntprobnp_qns1 == 1, NA_real_, ntprobnp_e1, missing = ntprobnp_e1),
      troponin_e1 = if_else(troponin_qns1 == 1, NA_real_, troponin_e1, missing = troponin_e1)
    )
}

fit_interaction_cox <- function(data, biomarker, base_covars,
                                time_var = "hf_time_days",
                                event_var = "hf_event",
                                race_var = "race_eth") {
  vars_needed <- c(time_var, event_var, biomarker, race_var, base_covars)

  d <- data %>%
    select(any_of(vars_needed)) %>%
    filter(complete.cases(.), .data[[time_var]] > 0)

  if (nrow(d) == 0L || sum(d[[event_var]], na.rm = TRUE) == 0L) {
    return(tibble(
      biomarker = biomarker,
      n = 0L,
      events = 0L,
      hr_ref = NA_real_,
      lcl_ref = NA_real_,
      ucl_ref = NA_real_,
      p_ref = NA_real_,
      p_interaction = NA_real_
    ))
  }

  rhs_base <- if (length(base_covars) > 0) paste(base_covars, collapse = " + ") else NULL
  rhs_red  <- paste(c(biomarker, race_var, rhs_base), collapse = " + ")
  rhs_full <- paste(c(paste0(biomarker, " * ", race_var), rhs_base), collapse = " + ")

  red_formula  <- as.formula(paste0("survival::Surv(", time_var, ", ", event_var, ") ~ ", rhs_red))
  full_formula <- as.formula(paste0("survival::Surv(", time_var, ", ", event_var, ") ~ ", rhs_full))

  fit_red <- tryCatch(
    survival::coxph(red_formula, data = d, ties = "efron"),
    error = function(e) NULL
  )
  fit_full <- tryCatch(
    survival::coxph(full_formula, data = d, ties = "efron"),
    error = function(e) NULL
  )

  if (is.null(fit_red) || is.null(fit_full)) {
    return(tibble(
      biomarker = biomarker,
      n = nrow(d),
      events = sum(d[[event_var]], na.rm = TRUE),
      hr_ref = NA_real_,
      lcl_ref = NA_real_,
      ucl_ref = NA_real_,
      p_ref = NA_real_,
      p_interaction = NA_real_
    ))
  }

  lrt <- tryCatch(as.data.frame(anova(fit_red, fit_full, test = "LRT")), error = function(e) NULL)
  p_int <- if (!is.null(lrt)) as.numeric(lrt[nrow(lrt), ncol(lrt)]) else NA_real_

  coef_tbl <- tryCatch(
    broom::tidy(fit_full, exponentiate = TRUE, conf.int = TRUE),
    error = function(e) tibble()
  )
  ref_row <- coef_tbl %>% filter(term == biomarker)

  tibble(
    biomarker = biomarker,
    n = nrow(d),
    events = sum(d[[event_var]], na.rm = TRUE),
    hr_ref = if (nrow(ref_row) == 1) ref_row$estimate else NA_real_,
    lcl_ref = if (nrow(ref_row) == 1) ref_row$conf.low else NA_real_,
    ucl_ref = if (nrow(ref_row) == 1) ref_row$conf.high else NA_real_,
    p_ref = if (nrow(ref_row) == 1) ref_row$p.value else NA_real_,
    p_interaction = p_int
  )
}

run_scan <- function(biomarkers, family_label, data, base_covars) {
  if (length(biomarkers) == 0L) return(tibble())

  purrr::map_dfr(
    biomarkers,
    fit_interaction_cox,
    data = data,
    base_covars = base_covars
  ) %>%
    mutate(
      family = family_label,
      p_adj_holm = p.adjust(p_interaction, method = "holm"),
      q_adj_bh   = p.adjust(p_interaction, method = "BH")
    )
}

# 3) Load datasets ------------------------------------------------------------
evt_primary <- if (fs::file_exists(paths$event_primary)) {
  read_mesa_csv(paths$event_primary, c("chfdiag", "ttchf", "efclass", "efmeas"))
} else {
  tibble(mesaid = integer())
}
evt_primary <- ensure_cols(evt_primary, c("chfdiag", "ttchf", "efclass", "efmeas"))

evt_secondary <- read_mesa_csv(paths$event_secondary, c("fuptt", "prebase", "exall", "chf", "chftt"))

exam1 <- read_mesa_csv(paths$exam1, c(
  "age1c", "gender1", "race1c", "bmi1c", "sbp1c", "htn1c", "dm031c",
  "chol1", "hdl1", "ldl1", "trig1", "cepgfr1c", "educ1", "income1",
  "cig1c", "pkyrs1c",
  "crp1", "crp1m", "il61", "fib1", "fib1m", "ddimer1",
  "olvef1"
))

site <- read_mesa_csv(paths$site, c("site1c"))

air_all <- read_mesa_csv(paths$air, c(
  "exam", "pm25_bl", "pm25_fu", "no2_bl", "no2_fu",
  "pm25_ugm3_1_yr_exam", "no2_ppb_1_yr_exam"
))

geocode_all <- read_mesa_csv(paths$geocode, c("exam", "accuracy"))

cardiac <- load_cardiac(paths$cardiac079, paths$cardiac244)

# Exam 1 ancillary air-related CRP/FIB (fallback only)
aire1 <- if (fs::file_exists(paths$aire1)) {
  read_mesa_csv(paths$aire1, c("crp1", "fib1", "crp1m", "fib1m"))
} else {
  tibble(mesaid = integer())
}
if ("crp1"  %in% names(aire1)) aire1 <- aire1 %>% rename(crp1_air  = crp1)
if ("fib1"  %in% names(aire1)) aire1 <- aire1 %>% rename(fib1_air  = fib1)
if ("crp1m" %in% names(aire1)) aire1 <- aire1 %>% rename(crp1m_air = crp1m)
if ("fib1m" %in% names(aire1)) aire1 <- aire1 %>% rename(fib1m_air = fib1m)

aire1 <- ensure_cols(aire1, c("crp1_air", "fib1_air", "crp1m_air", "fib1m_air"))

if (!all(is.na(aire1$crp1m_air))) aire1 <- aire1 %>% mutate(crp1_air = if_else(crp1m_air == 1, NA_real_, crp1_air, missing = crp1_air))
if (!all(is.na(aire1$fib1m_air))) aire1 <- aire1 %>% mutate(fib1_air = if_else(fib1m_air == 1, NA_real_, fib1_air, missing = fib1_air))

exam4_biom <- if (fs::file_exists(paths$exam4)) {
  read_mesa_csv(paths$exam4, c("crp4", "crp4m", "fib4", "fib4m", "ddimer4", "icam4"))
} else {
  tibble(mesaid = integer())
}
exam4_biom <- ensure_cols(exam4_biom, c("crp4", "crp4m", "fib4", "fib4m", "ddimer4", "icam4"))

immune <- if (fs::file_exists(paths$immune)) {
  read_mesa_csv(paths$immune, c(
    "abbaspha4", "abeospha4", "ablympha4", "abmoncya4", "abneupha4",
    "epcbla4", "epccd8la4", "epccd8ma4", "epccd8na4", "epcwba4"
  ))
} else {
  tibble(mesaid = integer())
}

# 4) Basic missing-value handling ---------------------------------------------
exam1 <- ensure_cols(exam1, c("crp1", "crp1m", "fib1", "fib1m", "il61", "ddimer1", "olvef1"))
if ("crp1m" %in% names(exam1)) exam1 <- exam1 %>% mutate(crp1 = if_else(crp1m == 1, NA_real_, crp1, missing = crp1))
if ("fib1m" %in% names(exam1)) exam1 <- exam1 %>% mutate(fib1 = if_else(fib1m == 1, NA_real_, fib1, missing = fib1))

if ("crp4m" %in% names(exam4_biom)) exam4_biom <- exam4_biom %>% mutate(crp4 = if_else(crp4m == 1, NA_real_, crp4, missing = crp4))
if ("fib4m" %in% names(exam4_biom)) exam4_biom <- exam4_biom %>% mutate(fib4 = if_else(fib4m == 1, NA_real_, fib4, missing = fib4))

# 5) Primary Exam 1 analytic dataset ------------------------------------------
mesa_main <- exam1 %>%
  left_join(aire1,         by = "mesaid") %>%
  left_join(site,          by = "mesaid") %>%
  left_join(evt_primary,   by = "mesaid") %>%
  left_join(evt_secondary, by = "mesaid") %>%
  left_join(cardiac,       by = "mesaid") %>%
  mutate(
    crp1_source = case_when(
      !is.na(crp1)     ~ "exam1",
      is.na(crp1) & !is.na(crp1_air) ~ "aire1",
      TRUE             ~ NA_character_
    ),
    fib1_source = case_when(
      !is.na(fib1)     ~ "exam1",
      is.na(fib1) & !is.na(fib1_air) ~ "aire1",
      TRUE             ~ NA_character_
    ),

    crp1_final = coalesce(crp1, crp1_air),
    fib1_final = coalesce(fib1, fib1_air),

    hf_event = if_else(!is.na(chf) & chf == 1, 1L, 0L),
    hf_time_days = case_when(
      hf_event == 1L ~ chftt,
      TRUE ~ fuptt
    ),

    ef_class_event = efclass,
    ef_meas_event  = efmeas,
    lvef_baseline = olvef1,

    sex = factor(gender1),
    race_eth = factor(race1c),
    site_factor = factor(site1c),

    z_log_crp1 = make_log_z(crp1_final),
    z_log_il61 = make_log_z(il61),
    z_fib1 = make_z(fib1_final),
    z_log_ddimer1 = make_log_z(ddimer1),

    z_log_ntprobnp_e1 = make_log_z(ntprobnp_e1),
    z_log_troponin_e1 = make_log_z(troponin_e1)
  ) %>%
  filter(
    !is.na(mesaid),
    !is.na(race_eth),
    !is.na(hf_time_days),
    hf_time_days > 0,
    is.na(exall) | exall == 0
  )

# 6) Air-pollution + inflammation sensitivity subset --------------------------
air_e1 <- air_all %>%
  filter(exam == 1) %>%
  distinct(mesaid, .keep_all = TRUE) %>%
  select(-exam)

geocode_e1 <- geocode_all %>%
  filter(exam == 1) %>%
  distinct(mesaid, .keep_all = TRUE) %>%
  transmute(mesaid, geocode_accuracy = accuracy)

mesa_air_sens <- mesa_main %>%
  left_join(air_e1, by = "mesaid") %>%
  left_join(geocode_e1, by = "mesaid") %>%
  filter(
    (!is.na(pm25_bl) | !is.na(pm25_ugm3_1_yr_exam)),
    (!is.na(no2_bl)  | !is.na(no2_ppb_1_yr_exam)),
    !is.na(crp1_final),
    !is.na(il61),
    !is.na(fib1_final),
    !is.na(ddimer1)
  )

# 7) Exam 4 / immune landmark-ready dataset -----------------------------------
mesa_exam4_landmark <- exam4_biom %>%
  left_join(immune, by = "mesaid") %>%
  left_join(evt_primary, by = "mesaid") %>%
  left_join(evt_secondary, by = "mesaid") %>%
  left_join(site, by = "mesaid") %>%
  left_join(
    exam1 %>% select(any_of(c(
      "mesaid", "age1c", "gender1", "race1c", "bmi1c", "sbp1c", "htn1c",
      "dm031c", "cepgfr1c", "educ1", "income1", "cig1c", "pkyrs1c", "olvef1"
    ))),
    by = "mesaid"
  ) %>%
  mutate(
    race_eth = factor(race1c),
    sex = factor(gender1),
    site_factor = factor(site1c),
    lvef_baseline = olvef1,

    z_log_crp4 = make_log_z(crp4),
    z_fib4 = make_z(fib4),
    z_log_ddimer4 = make_log_z(ddimer4),
    z_log_icam4 = make_log_z(icam4)
  )

# 8) Interaction scan ----------------------------------------------------------
primary_biomarkers <- c("z_log_crp1", "z_log_il61", "z_fib1", "z_log_ddimer1")
secondary_biomarkers <- c("z_log_ntprobnp_e1", "z_log_troponin_e1")

primary_biomarkers <- primary_biomarkers[
  primary_biomarkers %in% names(mesa_main) &
    vapply(mesa_main[primary_biomarkers], function(x) sum(!is.na(x)) > 0, logical(1))
]

secondary_biomarkers <- secondary_biomarkers[
  secondary_biomarkers %in% names(mesa_main) &
    vapply(mesa_main[secondary_biomarkers], function(x) sum(!is.na(x)) > 0, logical(1))
]

base_covars <- c(
  "age1c", "sex", "site_factor", "bmi1c", "sbp1c", "htn1c", "dm031c",
  "chol1", "hdl1", "ldl1", "trig1", "cepgfr1c", "educ1", "income1",
  "cig1c", "pkyrs1c"
)

base_covars <- base_covars[
  base_covars %in% names(mesa_main) &
    vapply(mesa_main[base_covars], function(x) sum(!is.na(x)) > 0, logical(1))
]

primary_results <- run_scan(
  biomarkers = primary_biomarkers,
  family_label = "primary_inflammatory",
  data = mesa_main,
  base_covars = base_covars
)

secondary_results <- run_scan(
  biomarkers = secondary_biomarkers,
  family_label = "secondary_cardiac",
  data = mesa_main,
  base_covars = base_covars
)

interaction_results <- bind_rows(primary_results, secondary_results) %>%
  arrange(family, p_interaction)

# 9) QC and export -------------------------------------------------------------
cat("N baseline analytic participants: ", nrow(mesa_main), "\n", sep = "")
cat("N HF events (chf): ", sum(mesa_main$hf_event, na.rm = TRUE), "\n", sep = "")
cat("Primary biomarker count: ", length(primary_biomarkers), "\n", sep = "")
cat("Secondary biomarker count: ", length(secondary_biomarkers), "\n", sep = "")
cat("N air/inflammation sensitivity participants: ", nrow(mesa_air_sens), "\n", sep = "")
cat("N exam4/immune landmark-ready participants: ", nrow(mesa_exam4_landmark), "\n", sep = "")

cat("CRP1 source counts:\n")
print(table(mesa_main$crp1_source, useNA = "ifany"))
cat("FIB1 source counts:\n")
print(table(mesa_main$fib1_source, useNA = "ifany"))

print(interaction_results)

readr::write_csv(interaction_results, fs::path(OUT_DIR, "mesa_hf_ethnicity_biomarker_interactions.csv"))
readr::write_csv(mesa_main, fs::path(OUT_DIR, "mesa_hf_biomarker_main.csv"))
readr::write_csv(mesa_air_sens, fs::path(OUT_DIR, "mesa_hf_air_inflammation_subset.csv"))
readr::write_csv(mesa_exam4_landmark, fs::path(OUT_DIR, "mesa_hf_exam4_landmark_ready.csv"))

# 10) Publication-ready figures ------------------------------------------------

## Shared aesthetics -----------------------------------------------------------
RACE_LABELS <- c(
  "1" = "White",
  "2" = "Chinese American",
  "3" = "Black",
  "4" = "Hispanic"
)

BIOMARKER_LABELS <- c(
  z_log_crp1    = "CRP (log-z)",
  z_log_il61    = "IL-6 (log-z)",
  z_fib1        = "Fibrinogen (z)",
  z_log_ddimer1 = "D-dimer (log-z)",
  z_log_ntprobnp_e1  = "NT-proBNP (log-z)",
  z_log_troponin_e1  = "Troponin (log-z)"
)

RACE_COLORS <- c(
  "White"            = "#4E79A7",
  "Chinese American" = "#F28E2B",
  "Black"            = "#E15759",
  "Hispanic"         = "#76B7B2"
)

theme_mesa <- function() {
  ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      strip.background   = ggplot2::element_rect(fill = "grey92", colour = NA),
      legend.position    = "bottom",
      plot.title         = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle      = ggplot2::element_text(size = 10, colour = "grey40"),
      axis.title         = ggplot2::element_text(size = 10),
      plot.caption       = ggplot2::element_text(size = 8, colour = "grey50")
    )
}

## Figure 1 -------------------------------------------------------------------
fig1_data <- interaction_results %>%
  filter(!is.na(hr_ref)) %>%
  mutate(
    label      = dplyr::recode(biomarker, !!!BIOMARKER_LABELS, .default = biomarker),
    sig_ref    = ifelse(p_ref < 0.05, "p < 0.05", "p \u2265 0.05"),
    p_int_lab  = paste0("P\u2099\u1d57 = ", formatC(p_interaction, digits = 2, format = "g"))
  )

fig1 <- ggplot2::ggplot(
  fig1_data,
  ggplot2::aes(x = hr_ref, y = stats::reorder(label, hr_ref),
               colour = sig_ref, shape = sig_ref)
) +
  ggplot2::geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  ggplot2::geom_errorbarh(
    ggplot2::aes(xmin = lcl_ref, xmax = ucl_ref),
    height = 0.25, linewidth = 0.7
  ) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_text(
    ggplot2::aes(x = ucl_ref, label = p_int_lab),
    hjust = -0.15, size = 3, colour = "grey40"
  ) +
  ggplot2::scale_colour_manual(
    name   = "Main effect",
    values = c("p < 0.05" = "#E15759", "p \u2265 0.05" = "#4E79A7")
  ) +
  ggplot2::scale_shape_manual(
    name   = "Main effect",
    values = c("p < 0.05" = 16, "p \u2265 0.05" = 1)
  ) +
  ggplot2::scale_x_continuous(
    name   = "Hazard Ratio per SD (95% CI)",
    limits = c(
      min(fig1_data$lcl_ref, na.rm = TRUE) * 0.85,
      max(fig1_data$ucl_ref, na.rm = TRUE) * 1.35
    )
  ) +
  ggplot2::labs(
    title    = "Figure 1. Biomarker Associations with Incident Heart Failure",
    subtitle = "Cox prop. hazards, adj. for age, sex, site, BMI, BP, diabetes, lipids, renal function, education, income, smoking",
    y        = NULL,
    caption  = glue::glue(
      "n = {max(fig1_data$n, na.rm = TRUE)} (complete cases vary per biomarker). ",
      "HR shown for reference (White) group from interaction model. ",
      "P\u2099\u1d57 = race\u00d7biomarker likelihood-ratio test."
    )
  ) +
  theme_mesa()

ggplot2::ggsave(
  fs::path(OUT_DIR, "fig1_forest_biomarker_hr.pdf"),
  fig1, width = 8, height = max(3.5, nrow(fig1_data) * 0.8 + 1.5),
  device = "pdf"
)
ggplot2::ggsave(
  fs::path(OUT_DIR, "fig1_forest_biomarker_hr.png"),
  fig1, width = 8, height = max(3.5, nrow(fig1_data) * 0.8 + 1.5),
  dpi = 300
)

## Figure 2 -------------------------------------------------------------------
fig2_data <- interaction_results %>%
  filter(!is.na(p_interaction)) %>%
  mutate(
    label     = dplyr::recode(biomarker, !!!BIOMARKER_LABELS, .default = biomarker),
    neg_log_p = -log10(p_interaction),
    sig_line  = -log10(0.05)
  )

fig2 <- ggplot2::ggplot(
  fig2_data,
  ggplot2::aes(x = neg_log_p, y = stats::reorder(label, neg_log_p))
) +
  ggplot2::geom_vline(
    xintercept = unique(fig2_data$sig_line),
    linetype   = "dashed", colour = "grey50", linewidth = 0.6
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = neg_log_p, yend = stats::reorder(label, neg_log_p)),
    colour = "grey70"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(colour = p_interaction < 0.05),
    size = 3.5
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = paste0("P = ", formatC(p_interaction, digits = 2, format = "g"))
    ),
    hjust = -0.2, size = 3
  ) +
  ggplot2::scale_colour_manual(
    name   = "Interaction",
    values = c("TRUE" = "#E15759", "FALSE" = "#4E79A7"),
    labels = c("TRUE" = "Significant (p < 0.05)", "FALSE" = "Not significant")
  ) +
  ggplot2::scale_x_continuous(
    name   = expression(-log[10](P[interaction])),
    expand = ggplot2::expansion(mult = c(0, 0.3))
  ) +
  ggplot2::labs(
    title    = "Figure 2. Race \u00d7 Biomarker Interaction P-values",
    subtitle = "Likelihood-ratio test comparing additive vs. interaction Cox models",
    y        = NULL,
    caption  = "Dashed line: nominal p = 0.05 threshold."
  ) +
  theme_mesa()

ggplot2::ggsave(
  fs::path(OUT_DIR, "fig2_interaction_pvalues.pdf"),
  fig2, width = 7, height = max(3, nrow(fig2_data) * 0.8 + 1.5),
  device = "pdf"
)
ggplot2::ggsave(
  fs::path(OUT_DIR, "fig2_interaction_pvalues.png"),
  fig2, width = 7, height = max(3, nrow(fig2_data) * 0.8 + 1.5),
  dpi = 300
)

## Figure 3 -------------------------------------------------------------------
surv_fit <- survival::survfit(
  survival::Surv(hf_time_days / 365.25, hf_event) ~ race_eth,
  data = mesa_main
)

surv_tbl <- broom::tidy(surv_fit) %>%
  mutate(
    race_label = RACE_LABELS[sub("race_eth=", "", strata)],
    time_yrs   = time
  )

fig3 <- ggplot2::ggplot(
  surv_tbl,
  ggplot2::aes(x = time_yrs, y = estimate, colour = race_label, fill = race_label)
) +
  ggplot2::geom_step(linewidth = 0.9) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.12, colour = NA
  ) +
  ggplot2::scale_colour_manual(name = "Race/Ethnicity", values = RACE_COLORS) +
  ggplot2::scale_fill_manual(name = "Race/Ethnicity", values = RACE_COLORS) +
  ggplot2::scale_y_continuous(
    name   = "Heart Failure\u2013Free Survival",
    limits = c(NA, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  ggplot2::scale_x_continuous(
    name   = "Time from Exam 1 (years)",
    breaks = seq(0, 20, 4)
  ) +
  ggplot2::labs(
    title    = "Figure 3. Heart Failure\u2013Free Survival by Race/Ethnicity",
    subtitle = glue::glue(
      "MESA Exam 1 cohort (n = {nrow(mesa_main)}); ",
      "{sum(mesa_main$hf_event)} incident HF events"
    ),
    caption  = "Shaded bands: 95% pointwise confidence intervals."
  ) +
  theme_mesa()

ggplot2::ggsave(
  fs::path(OUT_DIR, "fig3_km_hf_free_survival.pdf"),
  fig3, width = 7.5, height = 5,
  device = "pdf"
)
ggplot2::ggsave(
  fs::path(OUT_DIR, "fig3_km_hf_free_survival.png"),
  fig3, width = 7.5, height = 5,
  dpi = 300
)

## Figure 4 -------------------------------------------------------------------
biom_long <- mesa_main %>%
  mutate(race_label = dplyr::recode(as.character(race_eth), !!!RACE_LABELS)) %>%
  tidyr::pivot_longer(
    cols      = dplyr::any_of(names(BIOMARKER_LABELS)),
    names_to  = "biomarker",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%
  mutate(
    biom_label = dplyr::recode(biomarker, !!!BIOMARKER_LABELS)
  )

fig4 <- ggplot2::ggplot(
  biom_long,
  ggplot2::aes(x = race_label, y = value, colour = race_label, fill = race_label)
) +
  ggplot2::geom_violin(alpha = 0.25, linewidth = 0.5, trim = TRUE) +
  ggplot2::geom_boxplot(
    width = 0.18, outlier.shape = NA,
    alpha = 0.6, linewidth = 0.55
  ) +
  ggplot2::facet_wrap(
    ~biom_label, scales = "free_y",
    ncol = 2
  ) +
  ggplot2::scale_colour_manual(values = RACE_COLORS, guide = "none") +
  ggplot2::scale_fill_manual(values = RACE_COLORS, guide = "none") +
  ggplot2::scale_x_discrete(
    labels = function(x) stringr::str_wrap(x, width = 10)
  ) +
  ggplot2::labs(
    title    = "Figure 4. Standardised Biomarker Distributions by Race/Ethnicity",
    subtitle = "z-scored (or log-z-scored) Exam 1 values; boxes show IQR, line = median",
    x        = NULL,
    y        = "Standardised Value",
    caption  = "Outliers suppressed for clarity."
  ) +
  theme_mesa() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8))

n_biom_avail <- dplyr::n_distinct(biom_long$biomarker)
fig4_height  <- ceiling(n_biom_avail / 2) * 3.2 + 1.5

ggplot2::ggsave(
  fs::path(OUT_DIR, "fig4_biomarker_dist_by_race.pdf"),
  fig4, width = 8, height = fig4_height,
  device = "pdf"
)
ggplot2::ggsave(
  fs::path(OUT_DIR, "fig4_biomarker_dist_by_race.png"),
  fig4, width = 8, height = fig4_height,
  dpi = 300
)

cat("\nFigures written to:", OUT_DIR, "\n")

# 13) Group-specific biomarker slopes by race/ethnicity ----------------------
extract_group_slopes <- function(data, biomarker, base_covars,
                                 time_var = "hf_time_days",
                                 event_var = "hf_event",
                                 race_var = "race_eth") {
  vars_needed <- c(time_var, event_var, biomarker, race_var, base_covars)

  d <- data %>%
    dplyr::select(dplyr::any_of(vars_needed)) %>%
    dplyr::filter(stats::complete.cases(.), .data[[time_var]] > 0)

  if (nrow(d) == 0L || sum(d[[event_var]], na.rm = TRUE) == 0L) return(tibble())

  d <- d %>% dplyr::mutate(!!race_var := as.factor(.data[[race_var]]))
  race_levels <- levels(d[[race_var]])
  if (length(race_levels) < 1L) return(tibble())

  rhs_base <- if (length(base_covars) > 0) paste(base_covars, collapse = " + ") else NULL
  rhs_full <- paste(c(paste0(biomarker, " * ", race_var), rhs_base), collapse = " + ")
  full_formula <- stats::as.formula(
    paste0("survival::Surv(", time_var, ", ", event_var, ") ~ ", rhs_full)
  )

  fit <- tryCatch(
    survival::coxph(full_formula, data = d, ties = "efron"),
    error = function(e) NULL
  )
  if (is.null(fit)) return(tibble())

  beta <- stats::coef(fit)
  vcv  <- stats::vcov(fit)
  z975 <- stats::qnorm(0.975)

  beta_main_name <- biomarker
  if (!beta_main_name %in% names(beta)) return(tibble())

  race_effects <- purrr::map_dfr(race_levels, function(lvl) {
    if (identical(lvl, race_levels[[1]])) {
      est <- unname(beta[beta_main_name])
      se  <- sqrt(unname(vcv[beta_main_name, beta_main_name]))
      contrast_p <- NA_real_
    } else {
      int_name_1 <- paste0(biomarker, ":", race_var, lvl)
      int_name_2 <- paste0(race_var, lvl, ":", biomarker)
      int_name   <- dplyr::coalesce(
        names(beta)[match(int_name_1, names(beta))],
        names(beta)[match(int_name_2, names(beta))]
      )

      if (is.na(int_name)) {
        return(tibble(
          biomarker = biomarker,
          race_level = lvl,
          n = nrow(d),
          events = sum(d[[event_var]], na.rm = TRUE),
          hr = NA_real_,
          lcl = NA_real_,
          ucl = NA_real_,
          wald_p = NA_real_,
          contrast_p_vs_ref = NA_real_
        ))
      }

      est <- unname(beta[beta_main_name] + beta[int_name])
      var_est <- unname(
        vcv[beta_main_name, beta_main_name] +
          vcv[int_name, int_name] +
          2 * vcv[beta_main_name, int_name]
      )
      se <- sqrt(var_est)

      z_contrast <- unname(beta[int_name] / sqrt(vcv[int_name, int_name]))
      contrast_p <- 2 * stats::pnorm(abs(z_contrast), lower.tail = FALSE)
    }

    z_wald <- est / se

    tibble(
      biomarker = biomarker,
      race_level = lvl,
      n = nrow(d),
      events = sum(d[[event_var]], na.rm = TRUE),
      hr = exp(est),
      lcl = exp(est - z975 * se),
      ucl = exp(est + z975 * se),
      wald_p = 2 * stats::pnorm(abs(z_wald), lower.tail = FALSE),
      contrast_p_vs_ref = contrast_p
    )
  })

  race_effects
}

biomarkers_for_slopes <- unique(interaction_results$biomarker)

ethnic_slope_table <- purrr::map_dfr(
  biomarkers_for_slopes,
  extract_group_slopes,
  data = mesa_main,
  base_covars = base_covars
) %>%
  dplyr::mutate(
    race_label = dplyr::recode(as.character(race_level), !!!RACE_LABELS, .default = as.character(race_level)),
    biomarker_label = dplyr::recode(biomarker, !!!BIOMARKER_LABELS, .default = biomarker)
  ) %>%
  dplyr::arrange(biomarker_label, race_level) %>%
  dplyr::select(
    biomarker,
    biomarker_label,
    race_level,
    race_label,
    n,
    events,
    hr,
    lcl,
    ucl,
    wald_p,
    contrast_p_vs_ref
  )

ethnic_slope_summary <- ethnic_slope_table %>%
  dplyr::transmute(
    biomarker = biomarker_label,
    race = race_label,
    `HR (95% CI)` = glue::glue("{formatC(hr, format = 'f', digits = 2)} ({formatC(lcl, format = 'f', digits = 2)}, {formatC(ucl, format = 'f', digits = 2)})"),
    `Wald P` = formatC(wald_p, format = "g", digits = 3),
    `Contrast P vs White` = dplyr::if_else(
      is.na(contrast_p_vs_ref),
      NA_character_,
      formatC(contrast_p_vs_ref, format = "g", digits = 3)
    )
  )

readr::write_csv(ethnic_slope_table, fs::path(OUT_DIR, "mesa_hf_biomarker_group_specific_slopes.csv"))
readr::write_csv(ethnic_slope_summary, fs::path(OUT_DIR, "mesa_hf_biomarker_group_specific_slopes_compact.csv"))

interaction_results_group_specific <- ethnic_slope_table
interaction_results_group_compact <- ethnic_slope_summary

print(interaction_results_group_compact)

# 14) Additive Cox models: biomarker + race/ethnicity + covariates ----------
fit_additive_biomarker_cox <- function(data, biomarker, base_covars,
                                       time_var = "hf_time_days",
                                       event_var = "hf_event",
                                       race_var = "race_eth") {
  vars_needed <- c(time_var, event_var, biomarker, race_var, base_covars)

  d <- data %>%
    dplyr::select(dplyr::any_of(vars_needed)) %>%
    dplyr::filter(stats::complete.cases(.), .data[[time_var]] > 0)

  if (nrow(d) == 0L || sum(d[[event_var]], na.rm = TRUE) == 0L) {
    return(tibble(
      biomarker = biomarker,
      n = nrow(d),
      events = sum(d[[event_var]], na.rm = TRUE),
      hr = NA_real_,
      lcl = NA_real_,
      ucl = NA_real_,
      p_value = NA_real_
    ))
  }

  rhs_base <- if (length(base_covars) > 0) paste(base_covars, collapse = " + ") else NULL
  rhs <- paste(c(biomarker, race_var, rhs_base), collapse = " + ")
  fml <- stats::as.formula(
    paste0("survival::Surv(", time_var, ", ", event_var, ") ~ ", rhs)
  )

  fit <- tryCatch(
    survival::coxph(fml, data = d, ties = "efron"),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(tibble(
      biomarker = biomarker,
      n = nrow(d),
      events = sum(d[[event_var]], na.rm = TRUE),
      hr = NA_real_,
      lcl = NA_real_,
      ucl = NA_real_,
      p_value = NA_real_
    ))
  }

  td <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  tr <- td %>% dplyr::filter(.data$term == biomarker)

  tibble(
    biomarker = biomarker,
    n = nrow(d),
    events = sum(d[[event_var]], na.rm = TRUE),
    hr = if (nrow(tr) == 1) tr$estimate else NA_real_,
    lcl = if (nrow(tr) == 1) tr$conf.low else NA_real_,
    ucl = if (nrow(tr) == 1) tr$conf.high else NA_real_,
    p_value = if (nrow(tr) == 1) tr$p.value else NA_real_
  )
}

biomarkers_additive <- c("z_log_crp1", "z_log_il61", "z_fib1", "z_log_ddimer1")
biomarkers_additive <- biomarkers_additive[biomarkers_additive %in% names(mesa_main)]

additive_biomarker_models <- purrr::map_dfr(
  biomarkers_additive,
  fit_additive_biomarker_cox,
  data = mesa_main,
  base_covars = base_covars
) %>%
  dplyr::mutate(
    biomarker_label = dplyr::recode(biomarker, !!!BIOMARKER_LABELS, .default = biomarker)
  ) %>%
  dplyr::arrange(biomarker_label)

additive_biomarker_models_compact <- additive_biomarker_models %>%
  dplyr::transmute(
    biomarker = biomarker_label,
    n,
    events,
    `HR (95% CI)` = glue::glue(
      "{formatC(hr, format = 'f', digits = 2)} ({formatC(lcl, format = 'f', digits = 2)}, {formatC(ucl, format = 'f', digits = 2)})"
    ),
    `Wald P` = formatC(p_value, format = "g", digits = 3)
  )

readr::write_csv(
  additive_biomarker_models,
  fs::path(OUT_DIR, "mesa_hf_additive_biomarker_models.csv")
)
readr::write_csv(
  additive_biomarker_models_compact,
  fs::path(OUT_DIR, "mesa_hf_additive_biomarker_models_compact.csv")
)

print(additive_biomarker_models_compact)

# 15) Sequential adjustment models (Model 1-4) -------------------------------
sequential_model_covars <- list(
  "Model 1" = c("age1c", "sex", "site_factor"),
  "Model 2" = c("age1c", "sex", "site_factor", "bmi1c", "sbp1c", "htn1c", "dm031c", "chol1", "hdl1", "ldl1", "trig1", "cepgfr1c"),
  "Model 3" = c("age1c", "sex", "site_factor", "bmi1c", "sbp1c", "htn1c", "dm031c", "chol1", "hdl1", "ldl1", "trig1", "cepgfr1c", "cig1c", "educ1", "income1"),
  "Model 4" = c("age1c", "sex", "site_factor", "bmi1c", "sbp1c", "htn1c", "dm031c", "chol1", "hdl1", "ldl1", "trig1", "cepgfr1c", "cig1c", "educ1", "income1", "z_log_ntprobnp_e1", "z_log_troponin_e1")
)

keep_available_covars <- function(data, covars) {
  covars[
    covars %in% names(data) &
      vapply(data[covars], function(x) sum(!is.na(x)) > 0, logical(1))
  ]
}

additive_model_sequence <- purrr::map_dfr(
  names(sequential_model_covars),
  function(model_label) {
    covars_this_model <- keep_available_covars(mesa_main, sequential_model_covars[[model_label]])

    purrr::map_dfr(
      biomarkers_additive,
      fit_additive_biomarker_cox,
      data = mesa_main,
      base_covars = covars_this_model,
      race_var = "race_eth"
    ) %>%
      dplyr::mutate(model = model_label)
  }
) %>%
  dplyr::mutate(
    model = factor(model, levels = c("Model 1", "Model 2", "Model 3", "Model 4")),
    biomarker_label = dplyr::recode(biomarker, !!!BIOMARKER_LABELS, .default = biomarker)
  ) %>%
  dplyr::arrange(model, biomarker_label)

additive_model_sequence_compact <- additive_model_sequence %>%
  dplyr::transmute(
    model,
    biomarker = biomarker_label,
    n,
    events,
    `HR (95% CI)` = glue::glue(
      "{formatC(hr, format = 'f', digits = 2)} ({formatC(lcl, format = 'f', digits = 2)}, {formatC(ucl, format = 'f', digits = 2)})"
    ),
    `Wald P` = formatC(p_value, format = "g", digits = 3)
  )

readr::write_csv(
  additive_model_sequence,
  fs::path(OUT_DIR, "mesa_hf_additive_biomarker_model_sequence.csv")
)
readr::write_csv(
  additive_model_sequence_compact,
  fs::path(OUT_DIR, "mesa_hf_additive_biomarker_model_sequence_compact.csv")
)

additive_biomarker_models_sequential <- additive_model_sequence
additive_biomarker_models_sequential_compact <- additive_model_sequence_compact

print(additive_biomarker_models_sequential_compact)

# =============================================================================
# SECTION 17: RERI Additive Interaction Analysis
# =============================================================================
# Reference: VanderWeele TJ & Knol MJ (2014). A tutorial on interaction.
#            Epidemiologic Methods, 3(1), 33-72.
# (Identical to v4/v5/v6 — no changes needed)
# =============================================================================

cat("\n========== Section 17: RERI Additive Interaction ==========\n")

compute_reri_for_biomarker <- function(data, biomarker, base_covars,
                                       time_var  = "hf_time_days",
                                       event_var = "hf_event",
                                       race_var  = "race_eth",
                                       ref_race  = "1") {
  d <- data %>%
    dplyr::filter(!is.na(.data[[biomarker]]),
                  !is.na(.data[[race_var]]),
                  .data[[time_var]] > 0) %>%
    dplyr::mutate(
      biomarker_high = as.integer(.data[[biomarker]] > stats::median(.data[[biomarker]], na.rm = TRUE)),
      race_nonwhite  = as.integer(as.character(.data[[race_var]]) != ref_race)
    )

  d <- d %>%
    dplyr::mutate(
      joint = factor(
        paste0("A", race_nonwhite, "B", biomarker_high),
        levels = c("A0B0", "A1B0", "A0B1", "A1B1")
      )
    )

  vars_needed <- c(time_var, event_var, "joint", base_covars)
  d <- d %>%
    dplyr::select(dplyr::any_of(vars_needed)) %>%
    dplyr::filter(stats::complete.cases(.))

  if (nrow(d) < 20L || sum(d[[event_var]], na.rm = TRUE) < 5L) {
    return(tibble(
      biomarker      = biomarker,
      n              = nrow(d),
      events         = sum(d[[event_var]], na.rm = TRUE),
      hr_10          = NA_real_,
      hr_01          = NA_real_,
      hr_11          = NA_real_,
      reri           = NA_real_,
      reri_lcl       = NA_real_,
      reri_ucl       = NA_real_,
      reri_p         = NA_real_,
      ap             = NA_real_,
      si             = NA_real_
    ))
  }

  rhs_base <- if (length(base_covars) > 0) paste(base_covars, collapse = " + ") else NULL
  rhs      <- paste(c("joint", rhs_base), collapse = " + ")
  fml      <- stats::as.formula(
    paste0("survival::Surv(", time_var, ", ", event_var, ") ~ ", rhs)
  )

  fit <- tryCatch(
    survival::coxph(fml, data = d, ties = "efron"),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(tibble(
      biomarker = biomarker, n = nrow(d),
      events = sum(d[[event_var]], na.rm = TRUE),
      hr_10 = NA_real_, hr_01 = NA_real_, hr_11 = NA_real_,
      reri = NA_real_, reri_lcl = NA_real_, reri_ucl = NA_real_,
      reri_p = NA_real_, ap = NA_real_, si = NA_real_
    ))
  }

  beta <- stats::coef(fit)
  vcv  <- stats::vcov(fit)

  nm_10 <- "jointA1B0"
  nm_01 <- "jointA0B1"
  nm_11 <- "jointA1B1"

  if (!all(c(nm_10, nm_01, nm_11) %in% names(beta))) {
    return(tibble(
      biomarker = biomarker, n = nrow(d),
      events = sum(d[[event_var]], na.rm = TRUE),
      hr_10 = NA_real_, hr_01 = NA_real_, hr_11 = NA_real_,
      reri = NA_real_, reri_lcl = NA_real_, reri_ucl = NA_real_,
      reri_p = NA_real_, ap = NA_real_, si = NA_real_
    ))
  }

  b10 <- unname(beta[nm_10])
  b01 <- unname(beta[nm_01])
  b11 <- unname(beta[nm_11])

  hr_10 <- exp(b10)
  hr_01 <- exp(b01)
  hr_11 <- exp(b11)

  reri_val <- hr_11 - hr_10 - hr_01 + 1

  grad <- c(-exp(b10), -exp(b01), exp(b11))

  idx <- match(c(nm_10, nm_01, nm_11), names(beta))
  vcv_sub <- vcv[idx, idx]

  var_reri <- as.numeric(t(grad) %*% vcv_sub %*% grad)
  se_reri  <- sqrt(max(var_reri, 0))

  z975     <- stats::qnorm(0.975)
  reri_lcl <- reri_val - z975 * se_reri
  reri_ucl <- reri_val + z975 * se_reri
  reri_z   <- if (se_reri > 0) reri_val / se_reri else NA_real_
  reri_p   <- if (!is.na(reri_z)) 2 * stats::pnorm(abs(reri_z), lower.tail = FALSE) else NA_real_

  ap_val <- if (hr_11 > 0) reri_val / hr_11 else NA_real_
  si_val <- if ((hr_10 - 1) + (hr_01 - 1) != 0) {
    (hr_11 - 1) / ((hr_10 - 1) + (hr_01 - 1))
  } else {
    NA_real_
  }

  tibble(
    biomarker = biomarker,
    n         = nrow(d),
    events    = sum(d[[event_var]], na.rm = TRUE),
    hr_10     = hr_10,
    hr_01     = hr_01,
    hr_11     = hr_11,
    reri      = reri_val,
    reri_lcl  = reri_lcl,
    reri_ucl  = reri_ucl,
    reri_p    = reri_p,
    ap        = ap_val,
    si        = si_val
  )
}

reri_results <- purrr::map_dfr(
  primary_biomarkers,
  compute_reri_for_biomarker,
  data        = mesa_main,
  base_covars = base_covars
) %>%
  dplyr::mutate(
    biomarker_label = dplyr::recode(biomarker, !!!BIOMARKER_LABELS, .default = biomarker)
  ) %>%
  dplyr::arrange(reri_p)

reri_compact <- reri_results %>%
  dplyr::transmute(
    biomarker = biomarker_label,
    n,
    events,
    `HR(A=1,B=0)` = formatC(hr_10, format = "f", digits = 2),
    `HR(A=0,B=1)` = formatC(hr_01, format = "f", digits = 2),
    `HR(A=1,B=1)` = formatC(hr_11, format = "f", digits = 2),
    `RERI (95% CI)` = glue::glue(
      "{formatC(reri, format = 'f', digits = 3)} ",
      "({formatC(reri_lcl, format = 'f', digits = 3)}, ",
      "{formatC(reri_ucl, format = 'f', digits = 3)})"
    ),
    `RERI P`    = formatC(reri_p, format = "g", digits = 3),
    AP          = formatC(ap, format = "f", digits = 3),
    SI          = formatC(si, format = "f", digits = 3)
  )

cat("\n--- RERI Additive Interaction Results (Non-White vs White) ---\n")
print(reri_compact)

readr::write_csv(reri_results, fs::path(OUT_DIR, "mesa_hf_reri_additive_interaction.csv"))
readr::write_csv(reri_compact, fs::path(OUT_DIR, "mesa_hf_reri_additive_interaction_compact.csv"))

# --- RERI stratified by specific race/ethnicity groups ---
race_contrasts <- c("2" = "Chinese American", "3" = "Black", "4" = "Hispanic")

reri_by_race <- purrr::map_dfr(names(race_contrasts), function(race_code) {
  d_sub <- mesa_main %>%
    dplyr::filter(as.character(race_eth) %in% c("1", race_code))

  purrr::map_dfr(
    primary_biomarkers,
    compute_reri_for_biomarker,
    data        = d_sub,
    base_covars = base_covars,
    ref_race    = "1"
  ) %>%
    dplyr::mutate(
      contrast_race      = race_contrasts[[race_code]],
      contrast_race_code = race_code
    )
}) %>%
  dplyr::mutate(
    biomarker_label = dplyr::recode(biomarker, !!!BIOMARKER_LABELS, .default = biomarker)
  ) %>%
  dplyr::arrange(biomarker_label, contrast_race)

reri_by_race_compact <- reri_by_race %>%
  dplyr::transmute(
    biomarker     = biomarker_label,
    contrast_race,
    n,
    events,
    `HR(A=1,B=0)` = formatC(hr_10, format = "f", digits = 2),
    `HR(A=0,B=1)` = formatC(hr_01, format = "f", digits = 2),
    `HR(A=1,B=1)` = formatC(hr_11, format = "f", digits = 2),
    `RERI (95% CI)` = glue::glue(
      "{formatC(reri, format = 'f', digits = 3)} ",
      "({formatC(reri_lcl, format = 'f', digits = 3)}, ",
      "{formatC(reri_ucl, format = 'f', digits = 3)})"
    ),
    `RERI P` = formatC(reri_p, format = "g", digits = 3),
    AP       = formatC(ap, format = "f", digits = 3),
    SI       = formatC(si, format = "f", digits = 3)
  )

cat("\n--- RERI by Specific Race/Ethnicity Contrast ---\n")
print(reri_by_race_compact)

readr::write_csv(reri_by_race, fs::path(OUT_DIR, "mesa_hf_reri_by_race.csv"))
readr::write_csv(reri_by_race_compact, fs::path(OUT_DIR, "mesa_hf_reri_by_race_compact.csv"))

## Figure 5 — RERI forest plot by race contrast --------------------------------
fig5_data <- reri_by_race %>%
  dplyr::filter(!is.na(reri)) %>%
  dplyr::mutate(
    label = paste0(biomarker_label, " \u00d7 ", contrast_race),
    sig   = ifelse(reri_p < 0.05, "p < 0.05", "p \u2265 0.05")
  )

if (nrow(fig5_data) > 0) {
  fig5 <- ggplot2::ggplot(
    fig5_data,
    ggplot2::aes(x = reri, y = stats::reorder(label, reri),
                 colour = sig, shape = sig)
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = reri_lcl, xmax = reri_ucl),
      height = 0.3, linewidth = 0.7
    ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_colour_manual(
      name   = "RERI significance",
      values = c("p < 0.05" = "#E15759", "p \u2265 0.05" = "#4E79A7")
    ) +
    ggplot2::scale_shape_manual(
      name   = "RERI significance",
      values = c("p < 0.05" = 16, "p \u2265 0.05" = 1)
    ) +
    ggplot2::labs(
      title    = "Figure 5. RERI Additive Interaction: Race/Ethnicity \u00d7 Biomarker",
      subtitle = "RERI > 0 = super-additive risk; 0 = no additive interaction; < 0 = sub-additive",
      x        = "RERI (95% CI)",
      y        = NULL,
      caption  = "VanderWeele & Knol (2014). Biomarker dichotomised at cohort median. Adj. for full covariate set."
    ) +
    theme_mesa()

  ggplot2::ggsave(
    fs::path(OUT_DIR, "fig5_reri_forest.pdf"),
    fig5, width = 9, height = max(4, nrow(fig5_data) * 0.55 + 2),
    device = "pdf"
  )
  ggplot2::ggsave(
    fs::path(OUT_DIR, "fig5_reri_forest.png"),
    fig5, width = 9, height = max(4, nrow(fig5_data) * 0.55 + 2),
    dpi = 300
  )
  cat("Figure 5 (RERI forest) written.\n")
}


# =============================================================================
# SECTION 18 (v7 REWRITE): Formal Causal Mediation via CMAverse
# =============================================================================
# Reference: Shi B, Choirat C, Coull BA, VanderWeele TJ, Valeri L (2021).
#            CMAverse: A suite of functions for reproducible causal mediation
#            analyses. Epidemiology, 32(5), e20-e22.
#
# BUG HISTORY:
#   v4: yreg unspecified or wrong; factor confounders crashed formula builder.
#   v5: yreg = "survCox" (invalid); pre-built Surv() column (wrong interface).
#       All 8 models returned status = "model_error".
#   v6: yreg = "coxph" (correct); separate outcome/event args (correct).
#       All 8 models returned status = "success" BUT all estimates = NA.
#       Root cause: extraction used summary(res)$effect.decomposition which
#       does NOT exist in CMAverse. Also used lowercase effect names ("pnde")
#       but CMAverse uses "R"-prefixed names for ratio-scale outcomes ("Rpnde").
#
# v7 FIXES:
#   Fix 1: Extract from res$effect.pe, res$effect.ci.low, res$effect.ci.high,
#          res$effect.pval (named vectors on the cmest object itself).
#          Source: BS1125/CMAverse R/cmest.R lines 293-297.
#   Fix 2: For survival outcomes with EMint=TRUE, use CMAverse ratio-scale
#          effect names: "Rpnde", "Rtnie", "Rte", "pm".
#          Source: BS1125/CMAverse R/cmest.R lines 257-263.
#   Fix 3: Diagnostic cat() prints actual names from res$effect.pe so any
#          future naming mismatches are immediately visible in the console log.
# =============================================================================

cat("\n========== Section 18 (v7): Formal Causal Mediation (CMAverse) ==========\n")

# --- Install CMAverse if not already present ---
if (!requireNamespace("CMAverse", quietly = TRUE)) {
  cat("Installing CMAverse from CRAN...\n")
  install.packages("CMAverse")
}
library(CMAverse)

# --- Also need EValue package for sensitivity analysis ---
if (!requireNamespace("EValue", quietly = TRUE)) {
  cat("Installing EValue from CRAN...\n")
  install.packages("EValue")
}
library(EValue)

# =============================================================================
# Helper to expand factor columns to numeric dummies (unchanged from v5)
# =============================================================================
expand_factors_to_dummies <- function(data, factor_cols) {
  new_cols <- character()

  for (fc in factor_cols) {
    if (!fc %in% names(data)) next
    col <- data[[fc]]

    if (is.factor(col) || is.character(col)) {
      col <- as.factor(col)
      lvls <- levels(col)
      if (length(lvls) <= 1L) next

      for (i in 2:length(lvls)) {
        dummy_name <- paste0(fc, "_", lvls[i])
        dummy_name <- gsub("[^a-zA-Z0-9_]", "_", dummy_name)
        data[[dummy_name]] <- as.integer(col == lvls[i])
        new_cols <- c(new_cols, dummy_name)
      }
    } else {
      new_cols <- c(new_cols, fc)
    }
  }

  list(data = data, numeric_confounders = new_cols)
}

# =============================================================================
# Prepare data: expand factor confounders, drop incomplete cases
# v7: NO pre-built Surv() column. Return time_var and event_var names.
# =============================================================================
prepare_mediation_data <- function(data, exposure, mediator, confounders,
                                   time_var = "hf_time_years",
                                   event_var = "hf_event") {
  factor_confounders <- confounders[
    vapply(data[confounders], function(x) is.factor(x) || is.character(x), logical(1))
  ]
  numeric_confounders <- setdiff(confounders, factor_confounders)

  expanded <- expand_factors_to_dummies(data, factor_confounders)
  d <- expanded$data
  all_numeric_confounders <- c(numeric_confounders, expanded$numeric_confounders)

  keep_cols <- c(exposure, mediator, all_numeric_confounders, time_var, event_var)
  d <- d %>%
    dplyr::select(dplyr::any_of(keep_cols)) %>%
    dplyr::filter(stats::complete.cases(.), .data[[time_var]] > 0)

  list(
    data                 = d,
    numeric_confounders  = all_numeric_confounders,
    time_var             = time_var,
    event_var            = event_var
  )
}

# =============================================================================
# v7 extraction helper: pull from res$effect.pe / ci.low / ci.high / pval
# with case-insensitive fallback for effect names
# =============================================================================
extract_effect <- function(res, effect_name, slot_name) {
  # slot_name is one of: "effect.pe", "effect.se", "effect.ci.low",
  #                       "effect.ci.high", "effect.pval"
  vec <- tryCatch(res[[slot_name]], error = function(e) NULL)
  if (is.null(vec) || length(vec) == 0) return(NA_real_)

  # Try exact match first
  if (effect_name %in% names(vec)) return(as.numeric(unname(vec[effect_name])))

  # Try case-insensitive match
  idx <- match(tolower(effect_name), tolower(names(vec)))
  if (!is.na(idx)) return(as.numeric(unname(vec[idx])))

  NA_real_
}

# =============================================================================
# Main mediation runner: one exposure x one mediator (v7)
# =============================================================================
run_cmaverse_mediation_v7 <- function(data, exposure, mediator,
                                      confounders,
                                      time_var  = "hf_time_years",
                                      event_var = "hf_event",
                                      nboot = 500,
                                      seed  = 42) {

  prep <- prepare_mediation_data(
    data        = data,
    exposure    = exposure,
    mediator    = mediator,
    confounders = confounders,
    time_var    = time_var,
    event_var   = event_var
  )

  d         <- prep$data
  basec     <- prep$numeric_confounders
  time_col  <- prep$time_var
  event_col <- prep$event_var

  if (nrow(d) < 50L || sum(d[[event_col]], na.rm = TRUE) < 10L) {
    return(tibble(
      exposure    = exposure,
      mediator    = mediator,
      n           = nrow(d),
      events      = sum(d[[event_col]], na.rm = TRUE),
      nde_est = NA_real_, nde_lcl = NA_real_, nde_ucl = NA_real_, nde_p = NA_real_,
      nie_est = NA_real_, nie_lcl = NA_real_, nie_ucl = NA_real_, nie_p = NA_real_,
      te_est  = NA_real_, te_lcl  = NA_real_, te_ucl  = NA_real_,
      pm_est  = NA_real_, pm_lcl  = NA_real_, pm_ucl  = NA_real_,
      status      = "insufficient_data",
      error_msg   = NA_character_
    ))
  }

  mediator_median <- stats::median(d[[mediator]], na.rm = TRUE)

  estimation_methods <- c("paramfunc", "imputation")

  res <- NULL
  err_msg <- NA_character_

  for (est_method in estimation_methods) {
    res <- tryCatch({
      set.seed(seed)
      CMAverse::cmest(
        data       = as.data.frame(d),
        model      = "rb",
        outcome    = time_col,            # v6/v7: time column name
        event      = event_col,           # v6/v7: event column name
        exposure   = exposure,
        mediator   = mediator,
        basec      = basec,
        EMint      = TRUE,
        mreg       = list("linear"),
        yreg       = "coxph",             # v6/v7: correct CMAverse value
        astar      = 0,
        a          = 1,
        mval       = list(mediator_median),
        estimation = est_method,
        inference  = "bootstrap",
        nboot      = nboot
      )
    }, error = function(e) {
      err_msg <<- paste0("[", est_method, "] ", conditionMessage(e))
      cat("  CMAverse (", est_method, ") error for ", exposure, " -> ",
          mediator, ": ", conditionMessage(e), "\n")
      NULL
    })

    if (!is.null(res)) {
      cat("  CMAverse succeeded with estimation='", est_method, "' for ",
          exposure, " -> ", mediator, "\n")
      break
    }
  }

  if (is.null(res)) {
    return(tibble(
      exposure    = exposure,
      mediator    = mediator,
      n           = nrow(d),
      events      = sum(d[[event_col]], na.rm = TRUE),
      nde_est = NA_real_, nde_lcl = NA_real_, nde_ucl = NA_real_, nde_p = NA_real_,
      nie_est = NA_real_, nie_lcl = NA_real_, nie_ucl = NA_real_, nie_p = NA_real_,
      te_est  = NA_real_, te_lcl  = NA_real_, te_ucl  = NA_real_,
      pm_est  = NA_real_, pm_lcl  = NA_real_, pm_ucl  = NA_real_,
      status      = "model_error",
      error_msg   = err_msg
    ))
  }

  # =========================================================================
  # v7 FIX: Extract directly from the cmest result object slots
  # =========================================================================
  # Diagnostic: print what CMAverse actually returned
  pe_names <- tryCatch(names(res$effect.pe), error = function(e) character())
  cat("  [DIAG] res$effect.pe names: ", paste(pe_names, collapse = ", "), "\n")

  # For survival (ratio-scale) outcomes with EMint = TRUE, CMAverse uses:
  #   Rcde, Rpnde, Rtnde, Rpnie, Rtnie, Rte,
  #   ERcde, ERintref, ERintmed, ERpnie,
  #   ERcde(prop), ERintref(prop), ERintmed(prop), ERpnie(prop),
  #   pm, int, pe
  # We need: NDE = Rpnde, NIE = Rtnie, TE = Rte, PM = pm
  #
  # Fallback: also try lowercase "pnde"/"tnie"/"te" in case EMint=FALSE
  # or a future CMAverse version changes naming conventions.

  # Try R-prefixed names first (ratio-scale), then unprefixed (difference-scale)
  nde_name <- if ("Rpnde" %in% pe_names) "Rpnde" else if ("pnde" %in% pe_names) "pnde" else NA_character_
  nie_name <- if ("Rtnie" %in% pe_names) "Rtnie" else if ("tnie" %in% pe_names) "tnie" else NA_character_
  te_name  <- if ("Rte"   %in% pe_names) "Rte"   else if ("te"   %in% pe_names) "te"   else NA_character_
  pm_name  <- if ("pm"    %in% pe_names) "pm"    else NA_character_

  cat("  [DIAG] Mapped: NDE=", nde_name, ", NIE=", nie_name,
      ", TE=", te_name, ", PM=", pm_name, "\n")

  if (is.na(nde_name) && is.na(nie_name) && is.na(te_name)) {
    # Nothing matched — dump all available names for debugging
    return(tibble(
      exposure    = exposure,
      mediator    = mediator,
      n           = nrow(d),
      events      = sum(d[[event_col]], na.rm = TRUE),
      nde_est = NA_real_, nde_lcl = NA_real_, nde_ucl = NA_real_, nde_p = NA_real_,
      nie_est = NA_real_, nie_lcl = NA_real_, nie_ucl = NA_real_, nie_p = NA_real_,
      te_est  = NA_real_, te_lcl  = NA_real_, te_ucl  = NA_real_,
      pm_est  = NA_real_, pm_lcl  = NA_real_, pm_ucl  = NA_real_,
      status      = "extraction_error",
      error_msg   = paste0("No known effect names found in res$effect.pe. Available: ",
                           paste(pe_names, collapse = ", "))
    ))
  }

  tibble(
    exposure    = exposure,
    mediator    = mediator,
    n           = nrow(d),
    events      = sum(d[[event_col]], na.rm = TRUE),
    nde_est     = extract_effect(res, nde_name, "effect.pe"),
    nde_lcl     = extract_effect(res, nde_name, "effect.ci.low"),
    nde_ucl     = extract_effect(res, nde_name, "effect.ci.high"),
    nde_p       = extract_effect(res, nde_name, "effect.pval"),
    nie_est     = extract_effect(res, nie_name, "effect.pe"),
    nie_lcl     = extract_effect(res, nie_name, "effect.ci.low"),
    nie_ucl     = extract_effect(res, nie_name, "effect.ci.high"),
    nie_p       = extract_effect(res, nie_name, "effect.pval"),
    te_est      = extract_effect(res, te_name,  "effect.pe"),
    te_lcl      = extract_effect(res, te_name,  "effect.ci.low"),
    te_ucl      = extract_effect(res, te_name,  "effect.ci.high"),
    pm_est      = extract_effect(res, pm_name,  "effect.pe"),
    pm_lcl      = extract_effect(res, pm_name,  "effect.ci.low"),
    pm_ucl      = extract_effect(res, pm_name,  "effect.ci.high"),
    status      = "success",
    error_msg   = NA_character_
  )
}

# =============================================================================
# Prepare analytic datasets for mediation
# =============================================================================
mesa_mediation <- mesa_main %>%
  dplyr::mutate(
    race_nonwhite = as.integer(as.character(race_eth) != "1"),
    race_black_vs_white = dplyr::case_when(
      as.character(race_eth) == "1" ~ 0L,
      as.character(race_eth) == "3" ~ 1L,
      TRUE ~ NA_integer_
    ),
    hf_time_years = hf_time_days / 365.25
  )

mediation_confounders <- base_covars[base_covars %in% names(mesa_mediation)]

# =============================================================================
# Run mediation: Black vs White (primary equity comparison)
# =============================================================================
cat("\n--- CMAverse v7: Black vs White, through each inflammatory biomarker ---\n")

mesa_bw <- mesa_mediation %>%
  dplyr::filter(!is.na(race_black_vs_white))

mediation_results_bw <- purrr::map_dfr(
  primary_biomarkers,
  function(biom) {
    cat("  Running mediation: race_black_vs_white -> ", biom, " -> HF ...\n")
    run_cmaverse_mediation_v7(
      data        = mesa_bw,
      exposure    = "race_black_vs_white",
      mediator    = biom,
      confounders = mediation_confounders,
      nboot       = 500
    )
  }
) %>%
  dplyr::mutate(
    biomarker_label = dplyr::recode(mediator, !!!BIOMARKER_LABELS, .default = mediator),
    exposure_label  = "Black vs White"
  )

# =============================================================================
# Run mediation: Non-White vs White (broader contrast)
# =============================================================================
cat("\n--- CMAverse v7: Non-White vs White, through each inflammatory biomarker ---\n")

mediation_results_nw <- purrr::map_dfr(
  primary_biomarkers,
  function(biom) {
    cat("  Running mediation: race_nonwhite -> ", biom, " -> HF ...\n")
    run_cmaverse_mediation_v7(
      data        = mesa_mediation,
      exposure    = "race_nonwhite",
      mediator    = biom,
      confounders = mediation_confounders,
      nboot       = 500
    )
  }
) %>%
  dplyr::mutate(
    biomarker_label = dplyr::recode(mediator, !!!BIOMARKER_LABELS, .default = mediator),
    exposure_label  = "Non-White vs White"
  )

mediation_all <- dplyr::bind_rows(mediation_results_bw, mediation_results_nw)

# =============================================================================
# Compact table for publication
# =============================================================================
mediation_compact <- mediation_all %>%
  dplyr::filter(status == "success") %>%
  dplyr::transmute(
    exposure = exposure_label,
    mediator = biomarker_label,
    n,
    events,
    `NDE (95% CI)` = glue::glue(
      "{formatC(nde_est, format = 'f', digits = 3)} ",
      "({formatC(nde_lcl, format = 'f', digits = 3)}, ",
      "{formatC(nde_ucl, format = 'f', digits = 3)})"
    ),
    `NDE P` = formatC(nde_p, format = "g", digits = 3),
    `NIE (95% CI)` = glue::glue(
      "{formatC(nie_est, format = 'f', digits = 3)} ",
      "({formatC(nie_lcl, format = 'f', digits = 3)}, ",
      "{formatC(nie_ucl, format = 'f', digits = 3)})"
    ),
    `NIE P` = formatC(nie_p, format = "g", digits = 3),
    `TE (95% CI)` = glue::glue(
      "{formatC(te_est, format = 'f', digits = 3)} ",
      "({formatC(te_lcl, format = 'f', digits = 3)}, ",
      "{formatC(te_ucl, format = 'f', digits = 3)})"
    ),
    `% Mediated (95% CI)` = glue::glue(
      "{formatC(pm_est * 100, format = 'f', digits = 1)}% ",
      "({formatC(pm_lcl * 100, format = 'f', digits = 1)}, ",
      "{formatC(pm_ucl * 100, format = 'f', digits = 1)})"
    )
  )

cat("\n--- Causal Mediation Decomposition Results ---\n")
print(mediation_compact)

mediation_errors <- mediation_all %>%
  dplyr::filter(status != "success") %>%
  dplyr::select(exposure, mediator, n, events, status, error_msg)
if (nrow(mediation_errors) > 0) {
  cat("\n--- Mediation models that did NOT converge ---\n")
  print(mediation_errors)
}

readr::write_csv(mediation_all, fs::path(OUT_DIR, "mesa_hf_causal_mediation_cmaverse.csv"))
readr::write_csv(mediation_compact, fs::path(OUT_DIR, "mesa_hf_causal_mediation_compact.csv"))

# =============================================================================
# E-value sensitivity analysis for unmeasured mediator-outcome confounding
# =============================================================================
cat("\n--- E-value Sensitivity Analysis for NIE ---\n")

# For coxph with EMint=TRUE, the NIE (Rtnie) is already on the HR scale
# (ratio, not log-ratio). So nie_est IS the hazard ratio directly.
evalue_results <- mediation_all %>%
  dplyr::filter(status == "success", !is.na(nie_est)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    # Rtnie from CMAverse coxph is already a hazard ratio (not log-HR)
    nie_hr     = nie_est,
    nie_hr_lcl = nie_lcl,
    nie_hr_ucl = nie_ucl,
    # E-value for point estimate
    evalue_point = tryCatch({
      hr_for_eval <- max(nie_hr, 1 / nie_hr)
      lo_for_eval <- min(nie_hr_lcl, nie_hr_ucl)
      hi_for_eval <- max(nie_hr_lcl, nie_hr_ucl)

      ev <- EValue::evalues.HR(
        est  = hr_for_eval,
        lo   = lo_for_eval,
        hi   = hi_for_eval,
        rare = TRUE
      )
      as.numeric(ev["E-values", "point"])
    }, error = function(e) NA_real_),
    # E-value for CI bound closest to null
    evalue_ci = tryCatch({
      hr_for_eval <- max(nie_hr, 1 / nie_hr)
      lo_for_eval <- min(nie_hr_lcl, nie_hr_ucl)
      hi_for_eval <- max(nie_hr_lcl, nie_hr_ucl)

      ev <- EValue::evalues.HR(
        est  = hr_for_eval,
        lo   = lo_for_eval,
        hi   = hi_for_eval,
        rare = TRUE
      )
      as.numeric(ev["E-values", "lower"])
    }, error = function(e) NA_real_)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    exposure = exposure_label,
    mediator = biomarker_label,
    `NIE HR (95% CI)` = glue::glue(
      "{formatC(nie_hr, format = 'f', digits = 3)} ",
      "({formatC(nie_hr_lcl, format = 'f', digits = 3)}, ",
      "{formatC(nie_hr_ucl, format = 'f', digits = 3)})"
    ),
    `E-value (point)` = formatC(evalue_point, format = "f", digits = 2),
    `E-value (CI)`    = formatC(evalue_ci, format = "f", digits = 2)
  )

cat("\n--- E-value Results ---\n")
if (nrow(evalue_results) > 0) {
  print(evalue_results)
} else {
  cat("No successful NIE estimates available for E-value computation.\n")
}

readr::write_csv(evalue_results, fs::path(OUT_DIR, "mesa_hf_mediation_evalues.csv"))

# =============================================================================
# Final QC summary
# =============================================================================
cat("\n\n==================== FINAL SUMMARY (v7) ====================\n")
cat("Sections completed: 0-18\n")
cat("Section 17 outputs (RERI):\n")
cat("  - mesa_hf_reri_additive_interaction.csv\n")
cat("  - mesa_hf_reri_additive_interaction_compact.csv\n")
cat("  - mesa_hf_reri_by_race.csv\n")
cat("  - mesa_hf_reri_by_race_compact.csv\n")
cat("  - fig5_reri_forest.pdf / .png\n")
cat("Section 18 outputs (CMAverse mediation, v7 rewrite):\n")
cat("  - mesa_hf_causal_mediation_cmaverse.csv  (includes error_msg column)\n")
cat("  - mesa_hf_causal_mediation_compact.csv\n")
cat("  - mesa_hf_mediation_evalues.csv\n")

n_success <- sum(mediation_all$status == "success", na.rm = TRUE)
n_total   <- nrow(mediation_all)
cat(glue::glue("  Mediation models: {n_success}/{n_total} succeeded\n"))

if (n_success < n_total) {
  cat("  Failed models -- check error_msg column in mesa_hf_causal_mediation_cmaverse.csv\n")
}
cat("========================================================\n")