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
      ntprobnp_e1 = if_else(ntprobnp_qns1 == 1, NA_real_, ntprobnp_e1),
      troponin_e1 = if_else(troponin_qns1 == 1, NA_real_, troponin_e1)
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
  # Read values + (optional) missingness flags if present in this ancillary extract
  read_mesa_csv(paths$aire1, c("crp1", "fib1", "crp1m", "fib1m"))
} else {
  tibble(mesaid = integer())
}
if ("crp1"  %in% names(aire1)) aire1 <- aire1 %>% rename(crp1_air  = crp1)
if ("fib1"  %in% names(aire1)) aire1 <- aire1 %>% rename(fib1_air  = fib1)
if ("crp1m" %in% names(aire1)) aire1 <- aire1 %>% rename(crp1m_air = crp1m)
if ("fib1m" %in% names(aire1)) aire1 <- aire1 %>% rename(fib1m_air = fib1m)

aire1 <- ensure_cols(aire1, c("crp1_air", "fib1_air", "crp1m_air", "fib1m_air"))

# Apply ancillary missingness flags if provided (robust if all NA / not present)
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
# Exam 1: apply missingness flags from main Exam 1 file
exam1 <- ensure_cols(exam1, c("crp1", "crp1m", "fib1", "fib1m", "il61", "ddimer1", "olvef1"))
if ("crp1m" %in% names(exam1)) exam1 <- exam1 %>% mutate(crp1 = if_else(crp1m == 1, NA_real_, crp1, missing = crp1))
if ("fib1m" %in% names(exam1)) exam1 <- exam1 %>% mutate(fib1 = if_else(fib1m == 1, NA_real_, fib1, missing = fib1))

# Exam 4: apply missingness flags
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
    # provenance BEFORE coalesce (so it reflects the actual chosen value)
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

    # Main endpoint: overall HF from secondary events CSV
    hf_event = if_else(!is.na(chf) & chf == 1, 1L, 0L),
    hf_time_days = case_when(
      hf_event == 1L ~ chftt,
      TRUE ~ fuptt
    ),

    # EF subtype vars only if evt_primary exists (otherwise NA)
    ef_class_event = efclass,
    ef_meas_event  = efmeas,

    # Optional baseline LVEF (exam-based)
    lvef_baseline = olvef1,

    sex = factor(gender1),
    race_eth = factor(race1c),
    site_factor = factor(site1c),

    # Candidate biomarkers (Exam 1)
    z_log_crp1 = make_log_z(crp1_final),
    z_log_il61 = make_log_z(il61),
    z_fib1 = make_z(fib1_final),
    z_log_ddimer1 = make_log_z(ddimer1),

    # Cardiac biomarkers (Exam 1)
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

## Figure 1 — Forest plot: overall HR per SD biomarker ------------------------
fig1_data <- interaction_results %>%
  filter(!is.na(hr_ref)) %>%
  mutate(
    label      = dplyr::recode(biomarker, !!!BIOMARKER_LABELS, .default = biomarker),
    sig_ref    = ifelse(p_ref < 0.05, "p < 0.05", "p ≥ 0.05"),
    p_int_lab  = paste0("Pₙᵗ = ", formatC(p_interaction, digits = 2, format = "g"))
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
    values = c("p < 0.05" = "#E15759", "p ≥ 0.05" = "#4E79A7")
  ) +
  ggplot2::scale_shape_manual(
    name   = "Main effect",
    values = c("p < 0.05" = 16, "p ≥ 0.05" = 1)
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
      "Pₙᵗ = race×biomarker likelihood-ratio test."
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

## Figure 2 — Interaction p-value plot ----------------------------------------
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
    title    = "Figure 2. Race × Biomarker Interaction P-values",
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

## Figure 3 — Kaplan–Meier HF-free survival by race/ethnicity -----------------
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
    name   = "Heart Failure–Free Survival",
    limits = c(NA, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  ggplot2::scale_x_continuous(
    name   = "Time from Exam 1 (years)",
    breaks = seq(0, 20, 4)
  ) +
  ggplot2::labs(
    title    = "Figure 3. Heart Failure–Free Survival by Race/Ethnicity",
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

## Figure 4 — Biomarker distribution by race/ethnicity ------------------------
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

# 16) Incremental prognostic value: IL-6 and D-dimer -------------------------
fit_incremental_cox <- function(data, rhs,
                                time_var = "hf_time_days",
                                event_var = "hf_event") {
  needed_vars <- c(
    time_var,
    event_var,
    all.vars(stats::as.formula(paste0("~", rhs)))
  )

  dat <- data |>
    dplyr::select(dplyr::any_of(needed_vars)) |>
    dplyr::filter(
      dplyr::if_all(dplyr::all_of(needed_vars), ~ !is.na(.x)),
      .data[[time_var]] > 0,
      .data[[event_var]] %in% c(0, 1)
    )

  if (nrow(dat) == 0L) {
    rlang::abort(paste0("No complete observations for model: ", rhs))
  }

  fml <- stats::as.formula(
    paste0("survival::Surv(", time_var, ", ", event_var, ") ~ ", rhs)
  )

  survival::coxph(fml, data = dat, ties = "efron", x = TRUE, y = TRUE)
}

make_step_surv_fun <- function(survfit_obj) {
  if (length(survfit_obj$time) == 0L) {
    return(function(x) rep(1, length(x)))
  }

  stats::stepfun(
    x = survfit_obj$time,
    y = c(1, survfit_obj$surv),
    right = TRUE
  )
}

get_baseline_cumhaz_at <- function(fit, t) {
  bh <- survival::basehaz(fit, centered = FALSE)
  if (nrow(bh) == 0L) return(0)
  idx <- max(which(bh$time <= t))
  if (!is.finite(idx)) return(0)
  as.numeric(bh$hazard[idx])
}

predict_risk_at_time <- function(fit, data, t) {
  lp <- stats::predict(fit, newdata = data, type = "lp")
  h0_t <- get_baseline_cumhaz_at(fit, t)
  1 - exp(-h0_t * exp(lp))
}

calc_td_auc_ipcw <- function(score, time, event, t, g_hat_fun) {
  is_case <- time <= t & event == 1
  is_ctrl <- time > t

  if (sum(is_case) == 0L || sum(is_ctrl) == 0L) return(NA_real_)

  g_case <- pmax(g_hat_fun(pmax(time[is_case] - 1e-8, 0)), 1e-6)
  g_ctrl <- pmax(g_hat_fun(rep(t, sum(is_ctrl))), 1e-6)

  w_case <- 1 / g_case
  w_ctrl <- 1 / g_ctrl

  s_case <- score[is_case]
  s_ctrl <- score[is_ctrl]

  cmp_gt <- outer(s_case, s_ctrl, `>`)
  cmp_eq <- outer(s_case, s_ctrl, `==`)
  w_mat <- outer(w_case, w_ctrl)

  num <- sum((cmp_gt + 0.5 * cmp_eq) * w_mat, na.rm = TRUE)
  den <- sum(w_case) * sum(w_ctrl)

  if (!is.finite(den) || den <= 0) return(NA_real_)
  num / den
}

calc_brier_ipcw <- function(risk, time, event, t, g_hat_fun) {
  y_t <- as.integer(time <= t & event == 1)

  g_t <- pmax(g_hat_fun(rep(t, length(time))), 1e-6)
  g_time <- pmax(g_hat_fun(pmax(time - 1e-8, 0)), 1e-6)

  w <- dplyr::case_when(
    time <= t & event == 1 ~ 1 / g_time,
    time > t ~ 1 / g_t,
    TRUE ~ 0
  )

  denom <- sum(w, na.rm = TRUE)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)

  sum(w * (y_t - risk)^2, na.rm = TRUE) / denom
}

calc_ibs <- function(time, brier) {
  ok <- is.finite(time) & is.finite(brier)
  time <- as.numeric(time[ok])
  brier <- as.numeric(brier[ok])
  if (length(time) == 0L) return(NA_real_)

  ord <- order(time)
  time <- time[ord]
  brier <- brier[ord]

  if (length(unique(time)) < 2L) return(mean(brier, na.rm = TRUE))

  area <- sum(diff(time) * (head(brier, -1) + tail(brier, -1)) / 2)
  area / (max(time) - min(time))
}

compare_lrt <- function(fit_small, fit_large, model_small, model_large) {
  an <- tryCatch(as.data.frame(stats::anova(fit_small, fit_large, test = "LRT")), error = function(e) NULL)
  if (is.null(an) || nrow(an) == 0) {
    return(tibble::tibble(
      model_small = model_small,
      model_large = model_large,
      df_diff = NA_real_,
      chisq = NA_real_,
      p_value = NA_real_
    ))
  }

  last <- an[nrow(an), , drop = FALSE]
  nm <- names(last)
  nm_lower <- tolower(nm)

  df_col <- nm[which(nm_lower %in% c("df", "d.f."))[1]]
  chisq_col <- nm[which(nm_lower %in% c("chisq", "chi sq", "chi-square"))[1]]
  p_col <- nm[which(grepl("pr|p", nm_lower))[1]]

  tibble::tibble(
    model_small = model_small,
    model_large = model_large,
    df_diff = if (!is.na(df_col)) as.numeric(last[[df_col]]) else NA_real_,
    chisq = if (!is.na(chisq_col)) as.numeric(last[[chisq_col]]) else NA_real_,
    p_value = if (!is.na(p_col)) as.numeric(last[[p_col]]) else NA_real_
  )
}

incremental_covars <- c(
  "hf_time_days", "hf_event",
  "age1c", "sex", "race_eth", "site_factor",
  "bmi1c", "sbp1c", "htn1c", "dm031c", "chol1", "hdl1", "ldl1", "trig1", "cepgfr1c",
  "cig1c", "educ1", "income1",
  "z_log_ntprobnp_e1", "z_log_troponin_e1",
  "z_log_il61", "z_log_ddimer1"
)

incremental_covars <- incremental_covars[incremental_covars %in% names(mesa_main)]

incremental_data <- mesa_main |>
  dplyr::select(dplyr::any_of(incremental_covars)) |>
  dplyr::filter(
    !is.na(hf_time_days),
    !is.na(hf_event),
    hf_time_days > 0,
    hf_event %in% c(0, 1)
  )

rhs_M0 <- paste(c(
  "age1c", "sex", "race_eth", "site_factor",
  "bmi1c", "sbp1c", "htn1c", "dm031c", "chol1", "hdl1", "ldl1", "trig1", "cepgfr1c",
  "cig1c", "educ1", "income1"
), collapse = " + ")

rhs_M1 <- paste(rhs_M0, "z_log_il61", sep = " + ")
rhs_M2 <- paste(rhs_M0, "z_log_ddimer1", sep = " + ")
rhs_M3 <- paste(rhs_M0, "z_log_il61", "z_log_ddimer1", sep = " + ")
rhs_M4 <- paste(rhs_M0, "z_log_ntprobnp_e1", "z_log_troponin_e1", sep = " + ")
rhs_M5_il6 <- paste(rhs_M4, "z_log_il61", sep = " + ")
rhs_M5_ddimer <- paste(rhs_M4, "z_log_ddimer1", sep = " + ")
rhs_M5_both <- paste(rhs_M4, "z_log_il61", "z_log_ddimer1", sep = " + ")

n_complete_for_rhs <- function(data, rhs, time_var = "hf_time_days", event_var = "hf_event") {
  needed_vars <- c(time_var, event_var, all.vars(stats::as.formula(paste0("~", rhs))))
  data |>
    dplyr::select(dplyr::any_of(needed_vars)) |>
    dplyr::filter(
      dplyr::if_all(dplyr::all_of(needed_vars), ~ !is.na(.x)),
      .data[[time_var]] > 0,
      .data[[event_var]] %in% c(0, 1)
    ) |>
    nrow()
}

fit_incremental_cox_safe <- function(data, rhs, model_name,
                                     time_var = "hf_time_days",
                                     event_var = "hf_event") {
  n_cc <- n_complete_for_rhs(data, rhs, time_var, event_var)
  if (n_cc == 0L) return(NULL)
  fit_incremental_cox(data, rhs, time_var, event_var)
}

rhs_list <- list(
  M0_clinical = rhs_M0,
  M1_M0_plus_IL6 = rhs_M1,
  M2_M0_plus_Ddimer = rhs_M2,
  M3_M0_plus_IL6_Ddimer = rhs_M3,
  M4_M0_plus_NTproBNP_troponin = rhs_M4,
  M5_M4_plus_IL6 = rhs_M5_il6,
  M5_M4_plus_Ddimer = rhs_M5_ddimer,
  M5_M4_plus_IL6_Ddimer = rhs_M5_both
)

incremental_model_availability <- purrr::imap_dfr(
  rhs_list,
  \(rhs, nm) {
    tibble::tibble(
      model = nm,
      rhs = rhs,
      n_complete = n_complete_for_rhs(incremental_data, rhs),
      fit_ok = n_complete > 0
    )
  }
)

incremental_models <- purrr::imap(
  rhs_list,
  \(rhs, nm) fit_incremental_cox_safe(incremental_data, rhs, nm)
)

# keep only fitted models
incremental_models <- purrr::compact(incremental_models)

event_times <- incremental_data$hf_time_days[incremental_data$hf_event == 1]
eval_times <- stats::quantile(event_times, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
eval_times <- as.numeric(unique(eval_times[is.finite(eval_times) & eval_times > 0]))

g_fit <- survival::survfit(survival::Surv(hf_time_days, 1 - hf_event) ~ 1, data = incremental_data)
g_hat <- make_step_surv_fun(g_fit)

model_time_metrics <- purrr::imap_dfr(
  incremental_models,
  function(fit, model_name) {
    purrr::map_dfr(eval_times, function(t_now) {
      risk_t <- predict_risk_at_time(fit, incremental_data, t_now)

      tibble::tibble(
        model = model_name,
        time = t_now,
        td_auc = calc_td_auc_ipcw(
          score = risk_t,
          time = incremental_data$hf_time_days,
          event = incremental_data$hf_event,
          t = t_now,
          g_hat_fun = g_hat
        ),
        brier = calc_brier_ipcw(
          risk = risk_t,
          time = incremental_data$hf_time_days,
          event = incremental_data$hf_event,
          t = t_now,
          g_hat_fun = g_hat
        )
      )
    })
  }
)

auc_by_time <- model_time_metrics %>%
  dplyr::select(model, time, value = td_auc)

brier_by_time <- model_time_metrics %>%
  dplyr::select(model, time, value = brier)

auc_summary <- auc_by_time %>%
  dplyr::summarize(td_auc_mean = mean(value, na.rm = TRUE), .by = model)

ibs_summary <- brier_by_time %>%
  dplyr::summarize(integrated_brier = calc_ibs(time, value), .by = model)

cindex_summary <- purrr::imap_dfr(
  incremental_models,
  ~ tibble::tibble(
    model = .y,
    c_index = unname(summary(.x)$concordance[1])
  )
)

incremental_performance <- cindex_summary %>%
  dplyr::left_join(auc_summary, by = "model") %>%
  dplyr::left_join(ibs_summary, by = "model") %>%
  dplyr::mutate(
    n = nrow(incremental_data),
    events = sum(incremental_data$hf_event, na.rm = TRUE)
  ) %>%
  dplyr::arrange(model)

incremental_lrt <- dplyr::bind_rows(
  compare_lrt(incremental_models$M0_clinical, incremental_models$M1_M0_plus_IL6, "M0", "M1"),
  compare_lrt(incremental_models$M0_clinical, incremental_models$M2_M0_plus_Ddimer, "M0", "M2"),
  compare_lrt(incremental_models$M0_clinical, incremental_models$M3_M0_plus_IL6_Ddimer, "M0", "M3"),
  compare_lrt(incremental_models$M0_clinical, incremental_models$M4_M0_plus_NTproBNP_troponin, "M0", "M4"),
  compare_lrt(incremental_models$M4_M0_plus_NTproBNP_troponin, incremental_models$M5_M4_plus_IL6, "M4", "M5 (IL-6)"),
  compare_lrt(incremental_models$M4_M0_plus_NTproBNP_troponin, incremental_models$M5_M4_plus_Ddimer, "M4", "M5 (D-dimer)"),
  compare_lrt(incremental_models$M4_M0_plus_NTproBNP_troponin, incremental_models$M5_M4_plus_IL6_Ddimer, "M4", "M5 (IL-6 + D-dimer)")
) %>%
  dplyr::mutate(
    p_value_fmt = formatC(p_value, format = "g", digits = 3)
  )

incremental_performance_compact <- incremental_performance %>%
  dplyr::transmute(
    model,
    n,
    events,
    `C-index` = formatC(c_index, format = "f", digits = 3),
    `Time-dependent AUC (mean)` = formatC(td_auc_mean, format = "f", digits = 3),
    `Integrated Brier Score` = formatC(integrated_brier, format = "f", digits = 3)
  )

readr::write_csv(
  incremental_performance,
  fs::path(OUT_DIR, "mesa_hf_incremental_model_performance.csv")
)
readr::write_csv(
  incremental_performance_compact,
  fs::path(OUT_DIR, "mesa_hf_incremental_model_performance_compact.csv")
)
readr::write_csv(
  incremental_lrt,
  fs::path(OUT_DIR, "mesa_hf_incremental_model_lrt_improvement.csv")
)
readr::write_csv(
  model_time_metrics,
  fs::path(OUT_DIR, "mesa_hf_incremental_model_time_metrics.csv")
)

incremental_models_performance <- incremental_performance
incremental_models_performance_compact <- incremental_performance_compact
incremental_models_lrt <- incremental_lrt

print(incremental_models_performance_compact)
print(incremental_models_lrt)

# 17) Overlap diagnostics: inflammation/hemostatic markers -------------------
overlap_markers <- c("z_log_crp1", "z_log_il61", "z_fib1", "z_log_ddimer1")
overlap_markers <- overlap_markers[overlap_markers %in% names(mesa_main)]

overlap_data <- mesa_main %>%
  dplyr::select(dplyr::any_of(overlap_markers))

marker_pairs <- utils::combn(overlap_markers, 2, simplify = FALSE)

marker_pairwise_cor <- purrr::map_dfr(marker_pairs, function(vv) {
  x <- overlap_data[[vv[[1]]]]
  y <- overlap_data[[vv[[2]]]]

  keep <- is.finite(x) & is.finite(y)
  n_pair <- sum(keep)

  pearson_r <- if (n_pair >= 3L) stats::cor(x[keep], y[keep], method = "pearson") else NA_real_
  spearman_rho <- if (n_pair >= 3L) stats::cor(x[keep], y[keep], method = "spearman") else NA_real_

  tibble::tibble(
    marker_1 = vv[[1]],
    marker_2 = vv[[2]],
    n = n_pair,
    pearson_r = as.numeric(pearson_r),
    spearman_rho = as.numeric(spearman_rho)
  )
}) %>%
  dplyr::mutate(
    marker_1_label = dplyr::recode(marker_1, !!!BIOMARKER_LABELS, .default = marker_1),
    marker_2_label = dplyr::recode(marker_2, !!!BIOMARKER_LABELS, .default = marker_2)
  ) %>%
  dplyr::arrange(marker_1, marker_2)

calc_vif_from_design <- function(data, terms_vec) {
  terms_vec <- terms_vec[terms_vec %in% names(data)]
  if (length(terms_vec) == 0L) {
    return(tibble::tibble(term = character(), vif = numeric(), n_complete = integer()))
  }

  mm <- stats::model.matrix(
    stats::as.formula(paste0("~", paste(terms_vec, collapse = " + "))),
    data = data,
    na.action = stats::na.pass
  )

  mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  if (ncol(mm) == 0L) {
    return(tibble::tibble(term = character(), vif = numeric(), n_complete = integer()))
  }

  mm_df <- as.data.frame(mm)

  purrr::map_dfr(colnames(mm_df), function(target) {
    others <- setdiff(colnames(mm_df), target)

    cc <- stats::complete.cases(mm_df[, c(target, others), drop = FALSE])
    n_cc <- sum(cc)

    if (n_cc < 3L) {
      return(tibble::tibble(term = target, vif = NA_real_, n_complete = n_cc))
    }

    y <- mm_df[[target]][cc]
    if (!any(is.finite(y)) || stats::var(y, na.rm = TRUE) == 0) {
      return(tibble::tibble(term = target, vif = NA_real_, n_complete = n_cc))
    }

    if (length(others) == 0L) {
      return(tibble::tibble(term = target, vif = 1, n_complete = n_cc))
    }

    fit <- stats::lm(
      stats::reformulate(termlabels = others, response = target),
      data = mm_df[cc, , drop = FALSE]
    )

    r2 <- summary(fit)$r.squared
    vif_val <- if (is.finite(r2) && r2 < 1) 1 / (1 - r2) else if (r2 >= 1) Inf else NA_real_

    tibble::tibble(term = target, vif = as.numeric(vif_val), n_complete = n_cc)
  })
}

vif_terms <- c(
  "z_log_il61", "z_log_ddimer1",
  "age1c", "sex", "race_eth", "site_factor",
  "bmi1c", "sbp1c", "htn1c", "dm031c", "chol1", "hdl1", "ldl1", "trig1", "cepgfr1c",
  "cig1c", "educ1", "income1",
  "z_log_ntprobnp_e1", "z_log_troponin_e1"
)
vif_terms <- vif_terms[vif_terms %in% names(incremental_data)]

vif_table <- calc_vif_from_design(incremental_data, vif_terms) %>%
  dplyr::mutate(
    base_term = stringr::str_replace(term, "(sex|race_eth|site_factor).*$", "\\1"),
    marker = dplyr::if_else(base_term %in% c("z_log_il61", "z_log_ddimer1"), base_term, NA_character_),
    marker_label = dplyr::recode(marker, !!!BIOMARKER_LABELS, .default = marker)
  ) %>%
  dplyr::arrange(dplyr::desc(vif), term)

extract_biomarker_term <- function(fit, term_name) {
  td <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  tr <- td %>% dplyr::filter(.data$term == term_name)

  tibble::tibble(
    hr = if (nrow(tr) == 1) tr$estimate else NA_real_,
    lcl = if (nrow(tr) == 1) tr$conf.low else NA_real_,
    ucl = if (nrow(tr) == 1) tr$conf.high else NA_real_,
    p_value = if (nrow(tr) == 1) tr$p.value else NA_real_
  )
}

joint_signal_table <- dplyr::bind_rows(
  extract_biomarker_term(incremental_models$M1_M0_plus_IL6, "z_log_il61") %>% dplyr::mutate(model = "M1", biomarker = "z_log_il61"),
  extract_biomarker_term(incremental_models$M3_M0_plus_IL6_Ddimer, "z_log_il61") %>% dplyr::mutate(model = "M3", biomarker = "z_log_il61"),
  extract_biomarker_term(incremental_models$M2_M0_plus_Ddimer, "z_log_ddimer1") %>% dplyr::mutate(model = "M2", biomarker = "z_log_ddimer1"),
  extract_biomarker_term(incremental_models$M3_M0_plus_IL6_Ddimer, "z_log_ddimer1") %>% dplyr::mutate(model = "M3", biomarker = "z_log_ddimer1")
) %>%
  dplyr::mutate(
    biomarker_label = dplyr::recode(biomarker, !!!BIOMARKER_LABELS, .default = biomarker),
    hr_ci = glue::glue("{formatC(hr, format = 'f', digits = 2)} ({formatC(lcl, format = 'f', digits = 2)}, {formatC(ucl, format = 'f', digits = 2)})"),
    p_value_fmt = formatC(p_value, format = "g", digits = 3)
  ) %>%
  dplyr::select(model, biomarker, biomarker_label, hr, lcl, ucl, p_value, hr_ci, p_value_fmt)

joint_signal_lrt <- dplyr::bind_rows(
  compare_lrt(incremental_models$M1_M0_plus_IL6, incremental_models$M3_M0_plus_IL6_Ddimer, "M1", "M3"),
  compare_lrt(incremental_models$M2_M0_plus_Ddimer, incremental_models$M3_M0_plus_IL6_Ddimer, "M2", "M3")
) %>%
  dplyr::mutate(p_value_fmt = formatC(p_value, format = "g", digits = 3))

joint_signal_compact <- joint_signal_table %>%
  dplyr::transmute(
    model,
    biomarker = biomarker_label,
    `HR (95% CI)` = hr_ci,
    `Wald P` = p_value_fmt
  )

readr::write_csv(
  marker_pairwise_cor,
  fs::path(OUT_DIR, "mesa_hf_inflammation_hemostatic_pairwise_correlations.csv")
)
readr::write_csv(
  vif_table,
  fs::path(OUT_DIR, "mesa_hf_inflammation_hemostatic_vif.csv")
)
readr::write_csv(
  joint_signal_table,
  fs::path(OUT_DIR, "mesa_hf_il6_ddimer_joint_signal.csv")
)
readr::write_csv(
  joint_signal_lrt,
  fs::path(OUT_DIR, "mesa_hf_il6_ddimer_joint_signal_lrt.csv")
)

inflammation_hemostatic_pairwise_correlations <- marker_pairwise_cor
inflammation_hemostatic_vif <- vif_table
il6_ddimer_joint_signal <- joint_signal_table
il6_ddimer_joint_signal_lrt <- joint_signal_lrt

print(joint_signal_compact)
print(il6_ddimer_joint_signal_lrt)

# 18) Subgroup size, events, missingness, and biomarker distributions --------
biomarkers_power_check <- c("z_log_crp1", "z_log_il61", "z_fib1", "z_log_ddimer1")
biomarkers_power_check <- biomarkers_power_check[biomarkers_power_check %in% names(mesa_main)]

n_by_ethnicity <- mesa_main %>%
  dplyr::mutate(
    race_label = dplyr::recode(as.character(race_eth), !!!RACE_LABELS, .default = as.character(race_eth))
  ) %>%
  dplyr::summarize(
    n = dplyr::n(),
    .by = c(race_eth, race_label)
  ) %>%
  dplyr::arrange(race_eth)

hf_events_by_ethnicity <- mesa_main %>%
  dplyr::mutate(
    race_label = dplyr::recode(as.character(race_eth), !!!RACE_LABELS, .default = as.character(race_eth))
  ) %>%
  dplyr::summarize(
    n = dplyr::n(),
    hf_events = sum(hf_event, na.rm = TRUE),
    hf_event_pct = 100 * hf_events / n,
    .by = c(race_eth, race_label)
  ) %>%
  dplyr::arrange(race_eth)

missingness_by_biomarker_ethnicity <- mesa_main %>%
  dplyr::mutate(
    race_label = dplyr::recode(as.character(race_eth), !!!RACE_LABELS, .default = as.character(race_eth))
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::any_of(biomarkers_power_check),
    names_to = "biomarker",
    values_to = "value"
  ) %>%
  dplyr::summarize(
    n = dplyr::n(),
    n_missing = sum(is.na(value)),
    pct_missing = 100 * n_missing / n,
    n_non_missing = sum(!is.na(value)),
    .by = c(race_eth, race_label, biomarker)
  ) %>%
  dplyr::mutate(
    biomarker_label = dplyr::recode(biomarker, !!!BIOMARKER_LABELS, .default = biomarker)
  ) %>%
  dplyr::arrange(biomarker, race_eth)

biomarker_distribution_by_ethnicity <- mesa_main %>%
  dplyr::mutate(
    race_label = dplyr::recode(as.character(race_eth), !!!RACE_LABELS, .default = as.character(race_eth))
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::any_of(biomarkers_power_check),
    names_to = "biomarker",
    values_to = "value"
  ) %>%
  dplyr::summarize(
    n_non_missing = sum(!is.na(value)),
    median = stats::median(value, na.rm = TRUE),
    q1 = stats::quantile(value, probs = 0.25, na.rm = TRUE, names = FALSE),
    q3 = stats::quantile(value, probs = 0.75, na.rm = TRUE, names = FALSE),
    .by = c(race_eth, race_label, biomarker)
  ) %>%
  dplyr::mutate(
    biomarker_label = dplyr::recode(biomarker, !!!BIOMARKER_LABELS, .default = biomarker),
    median_iqr = glue::glue(
      "{formatC(median, format = 'f', digits = 2)} ({formatC(q1, format = 'f', digits = 2)}, {formatC(q3, format = 'f', digits = 2)})"
    )
  ) %>%
  dplyr::arrange(biomarker, race_eth)

power_summary_compact <- n_by_ethnicity %>%
  dplyr::left_join(
    hf_events_by_ethnicity %>% dplyr::select(race_eth, hf_events, hf_event_pct),
    by = "race_eth"
  )

readr::write_csv(
  n_by_ethnicity,
  fs::path(OUT_DIR, "mesa_hf_n_by_ethnicity.csv")
)
readr::write_csv(
  hf_events_by_ethnicity,
  fs::path(OUT_DIR, "mesa_hf_events_by_ethnicity.csv")
)
readr::write_csv(
  missingness_by_biomarker_ethnicity,
  fs::path(OUT_DIR, "mesa_hf_biomarker_missingness_by_ethnicity.csv")
)
readr::write_csv(
  biomarker_distribution_by_ethnicity,
  fs::path(OUT_DIR, "mesa_hf_biomarker_distribution_by_ethnicity.csv")
)

subgroup_n_by_ethnicity <- n_by_ethnicity
subgroup_hf_events_by_ethnicity <- hf_events_by_ethnicity
subgroup_missingness_by_biomarker_ethnicity <- missingness_by_biomarker_ethnicity
subgroup_biomarker_distribution_by_ethnicity <- biomarker_distribution_by_ethnicity

print(power_summary_compact)

# 19) Multiple imputation of covariates with biomarkers left observed --------
extract_pool_p_value <- function(test_obj) {
  if (is.null(test_obj)) return(NA_real_)

  test_df <- tryCatch(as.data.frame(test_obj), error = function(e) NULL)
  if (!is.null(test_df) && nrow(test_df) > 0) {
    p_col <- names(test_df)[grepl("pr|p", names(test_df), ignore.case = TRUE)][1]
    if (!is.na(p_col)) return(as.numeric(test_df[[p_col]][1]))
  }

  if (is.list(test_obj)) {
    for (nm in c("p.value", "pval", "p")) {
      if (!is.null(test_obj[[nm]])) return(as.numeric(test_obj[[nm]][1]))
    }
  }

  NA_real_
}

fit_mi_additive_cox <- function(mids, biomarker, base_covars,
                                time_var = "hf_time_days",
                                event_var = "hf_event",
                                race_var = "race_eth") {
  vars_needed <- unique(c(time_var, event_var, biomarker, race_var, base_covars))
  rhs <- paste(c(biomarker, race_var, base_covars), collapse = " + ")
  fml_txt <- paste0("survival::Surv(", time_var, ", ", event_var, ") ~ ", rhs)

  mira_fit <- with(
    mids, {
      dat <- as.data.frame(mget(vars_needed, inherits = TRUE))
      survival::coxph(stats::as.formula(fml_txt), data = dat, ties = "efron")
    }
  )

  pooled <- mice::pool(mira_fit)
  pooled_sum <- summary(pooled, conf.int = TRUE, exponentiate = TRUE)
  tr <- pooled_sum |> dplyr::filter(.data$term == biomarker)

  tibble::tibble(
    biomarker = biomarker,
    hr = if (nrow(tr) == 1) tr$estimate else NA_real_,
    lcl = if (nrow(tr) == 1) tr$`2.5 %` else NA_real_,
    ucl = if (nrow(tr) == 1) tr$`97.5 %` else NA_real_,
    p_value = if (nrow(tr) == 1) tr$p.value else NA_real_
  )
}

fit_mi_interaction_cox <- function(mids, biomarker, base_covars,
                                   time_var = "hf_time_days",
                                   event_var = "hf_event",
                                   race_var = "race_eth") {
  vars_needed <- unique(c(time_var, event_var, biomarker, race_var, base_covars))

  rhs_red <- paste(c(biomarker, race_var, base_covars), collapse = " + ")
  rhs_full <- paste(c(paste0(biomarker, " * ", race_var), base_covars), collapse = " + ")

  fml_red_txt <- paste0("survival::Surv(", time_var, ", ", event_var, ") ~ ", rhs_red)
  fml_full_txt <- paste0("survival::Surv(", time_var, ", ", event_var, ") ~ ", rhs_full)

  mira_red <- with(
    mids, {
      dat <- as.data.frame(mget(vars_needed, inherits = TRUE))
      survival::coxph(stats::as.formula(fml_red_txt), data = dat, ties = "efron")
    }
  )

  mira_full <- with(
    mids, {
      dat <- as.data.frame(mget(vars_needed, inherits = TRUE))
      survival::coxph(stats::as.formula(fml_full_txt), data = dat, ties = "efron")
    }
  )

  pooled_full <- mice::pool(mira_full)
  pooled_sum <- summary(pooled_full, conf.int = TRUE, exponentiate = TRUE)
  tr <- pooled_sum |> dplyr::filter(.data$term == biomarker)

  d1_test <- tryCatch(mice::D1(mira_full, mira_red), error = function(e) NULL)

  tibble::tibble(
    biomarker = biomarker,
    hr_ref = if (nrow(tr) == 1) tr$estimate else NA_real_,
    lcl_ref = if (nrow(tr) == 1) tr$`2.5 %` else NA_real_,
    ucl_ref = if (nrow(tr) == 1) tr$`97.5 %` else NA_real_,
    p_ref = if (nrow(tr) == 1) tr$p.value else NA_real_,
    p_interaction = extract_pool_p_value(d1_test)
  )
}

biomarkers_mi <- primary_biomarkers[primary_biomarkers %in% names(mesa_main)]

mi_covars <- unique(c(
  "hf_time_days", "hf_event", "race_eth",
  biomarkers_mi,
  base_covars
))
mi_covars <- mi_covars[mi_covars %in% names(mesa_main)]

mi_data <- mesa_main %>%
  dplyr::select(dplyr::any_of(mi_covars)) %>%
  dplyr::filter(
    !is.na(hf_time_days),
    hf_time_days > 0,
    !is.na(hf_event),
    !is.na(race_eth)
  )

mi_init <- mice::mice(mi_data, maxit = 0, printFlag = FALSE)
mi_method <- mi_init$method
mi_predictor_matrix <- mi_init$predictorMatrix

vars_not_imputed <- intersect(c("hf_time_days", "hf_event", "race_eth", biomarkers_mi), names(mi_data))
mi_method[vars_not_imputed] <- ""

mi_predictor_matrix[, c("hf_time_days", "hf_event")] <- 1
diag(mi_predictor_matrix) <- 0

mi_fit <- mice::mice(
  mi_data,
  m = 20,
  maxit = 10,
  method = mi_method,
  predictorMatrix = mi_predictor_matrix,
  seed = 20260313,
  printFlag = FALSE
)

mi_pooled_results <- purrr::map_dfr(
  biomarkers_mi,
  fit_mi_additive_cox,
  mids = mi_fit,
  base_covars = base_covars
) %>%
  dplyr::mutate(
    biomarker_label = dplyr::recode(biomarker, !!!BIOMARKER_LABELS, .default = biomarker),
    hr_ci = glue::glue(
      "{formatC(hr, format = 'f', digits = 2)} ({formatC(lcl, format = 'f', digits = 2)}, {formatC(ucl, format = 'f', digits = 2)})"
    ),
    p_value_fmt = formatC(p_value, format = "g", digits = 3)
  )

mi_interaction_results <- purrr::map_dfr(
  biomarkers_mi,
  fit_mi_interaction_cox,
  mids = mi_fit,
  base_covars = base_covars
) %>%
  dplyr::mutate(
    family = "primary_inflammatory",
    biomarker_label = dplyr::recode(biomarker, !!!BIOMARKER_LABELS, .default = biomarker),
    hr_ref_ci = glue::glue(
      "{formatC(hr_ref, format = 'f', digits = 2)} ({formatC(lcl_ref, format = 'f', digits = 2)}, {formatC(ucl_ref, format = 'f', digits = 2)})"
    ),
    p_ref_fmt = formatC(p_ref, format = "g", digits = 3),
    p_interaction_fmt = formatC(p_interaction, format = "g", digits = 3)
  )

mi_vs_cc_pooled <- additive_biomarker_models %>%
  dplyr::filter(.data$biomarker %in% biomarkers_mi) %>%
  dplyr::select(
    biomarker,
    biomarker_label,
    cc_hr = hr,
    cc_lcl = lcl,
    cc_ucl = ucl,
    cc_p_value = p_value
  ) %>%
  dplyr::full_join(
    mi_pooled_results %>%
      dplyr::select(
        biomarker,
        biomarker_label,
        mi_hr = hr,
        mi_lcl = lcl,
        mi_ucl = ucl,
        mi_p_value = p_value
      ),
    by = c("biomarker", "biomarker_label")
  ) %>%
  dplyr::mutate(
    cc_hr_ci = glue::glue(
      "{formatC(cc_hr, format = 'f', digits = 2)} ({formatC(cc_lcl, format = 'f', digits = 2)}, {formatC(cc_ucl, format = 'f', digits = 2)})"
    ),
    mi_hr_ci = glue::glue(
      "{formatC(mi_hr, format = 'f', digits = 2)} ({formatC(mi_lcl, format = 'f', digits = 2)}, {formatC(mi_ucl, format = 'f', digits = 2)})"
    ),
    cc_p_value_fmt = formatC(cc_p_value, format = "g", digits = 3),
    mi_p_value_fmt = formatC(mi_p_value, format = "g", digits = 3)
  )

mi_vs_cc_interaction <- interaction_results %>%
  dplyr::filter(.data$biomarker %in% biomarkers_mi) %>%
  dplyr::select(
    biomarker,
    cc_hr_ref = hr_ref,
    cc_lcl_ref = lcl_ref,
    cc_ucl_ref = ucl_ref,
    cc_p_ref = p_ref,
    cc_p_interaction = p_interaction
  ) %>%
  dplyr::full_join(
    mi_interaction_results %>%
      dplyr::select(
        biomarker,
        biomarker_label,
        mi_hr_ref = hr_ref,
        mi_lcl_ref = lcl_ref,
        mi_ucl_ref = ucl_ref,
        mi_p_ref = p_ref,
        mi_p_interaction = p_interaction
      ),
    by = "biomarker"
  ) %>%
  dplyr::mutate(
    biomarker_label = dplyr::coalesce(biomarker_label, dplyr::recode(biomarker, !!!BIOMARKER_LABELS, .default = biomarker)),
    cc_hr_ref_ci = glue::glue(
      "{formatC(cc_hr_ref, format = 'f', digits = 2)} ({formatC(cc_lcl_ref, format = 'f', digits = 2)}, {formatC(cc_ucl_ref, format = 'f', digits = 2)})"
    ),
    mi_hr_ref_ci = glue::glue(
      "{formatC(mi_hr_ref, format = 'f', digits = 2)} ({formatC(mi_lcl_ref, format = 'f', digits = 2)}, {formatC(mi_ucl_ref, format = 'f', digits = 2)})"
    ),
    cc_p_ref_fmt = formatC(cc_p_ref, format = "g", digits = 3),
    mi_p_ref_fmt = formatC(mi_p_ref, format = "g", digits = 3),
    cc_p_interaction_fmt = formatC(cc_p_interaction, format = "g", digits = 3),
    mi_p_interaction_fmt = formatC(mi_p_interaction, format = "g", digits = 3)
  )

mi_vs_cc_pooled_compact <- mi_vs_cc_pooled %>%
  dplyr::transmute(
    biomarker = biomarker_label,
    `Complete-case HR (95% CI)` = cc_hr_ci,
    `Complete-case Wald P` = cc_p_value_fmt,
    `MI HR (95% CI)` = mi_hr_ci,
    `MI Wald P` = mi_p_value_fmt
  )

mi_vs_cc_interaction_compact <- mi_vs_cc_interaction %>%
  dplyr::transmute(
    biomarker = biomarker_label,
    `Complete-case ref HR (95% CI)` = cc_hr_ref_ci,
    `Complete-case ref P` = cc_p_ref_fmt,
    `Complete-case interaction P` = cc_p_interaction_fmt,
    `MI ref HR (95% CI)` = mi_hr_ref_ci,
    `MI ref P` = mi_p_ref_fmt,
    `MI interaction P` = mi_p_interaction_fmt
  )

readr::write_csv(
  mi_pooled_results,
  fs::path(OUT_DIR, "mesa_hf_mi_pooled_biomarker_models.csv")
)
readr::write_csv(
  mi_interaction_results,
  fs::path(OUT_DIR, "mesa_hf_mi_interaction_biomarker_models.csv")
)
readr::write_csv(
  mi_vs_cc_pooled,
  fs::path(OUT_DIR, "mesa_hf_mi_vs_complete_case_pooled.csv")
)
readr::write_csv(
  mi_vs_cc_interaction,
  fs::path(OUT_DIR, "mesa_hf_mi_vs_complete_case_interaction.csv")
)

mi_covariate_imputation <- mi_fit
mi_pooled_biomarker_models <- mi_pooled_results
mi_interaction_biomarker_models <- mi_interaction_results
mi_vs_complete_case_pooled <- mi_vs_cc_pooled
mi_vs_complete_case_interaction <- mi_vs_cc_interaction

print(mi_vs_cc_pooled_compact)
print(mi_vs_cc_interaction_compact)

# 20) Environmental robustness: PM2.5, NO2, geocode accuracy ----------------
infer_geocode_restriction <- function(x) {
  if (is.numeric(x)) {
    xx <- x[is.finite(x)]
    if (length(xx) == 0L) return(rep(FALSE, length(x)))

    u <- sort(unique(xx))
    if (all(abs(u - round(u)) < 1e-8) && min(u) >= 1 && max(u) <= 10) {
      return(!is.na(x) & x <= 2)
    }

    med <- stats::median(xx, na.rm = TRUE)
    return(!is.na(x) & x <= med)
  }

  if (is.character(x) || is.factor(x)) {
    xc <- tolower(as.character(x))
    keep <- grepl("rooftop|parcel|exact|address", xc)
    return(!is.na(xc) & keep)
  }

  !is.na(x)
}

inflammatory_markers_env <- c("z_log_crp1", "z_log_il61", "z_fib1", "z_log_ddimer1")
inflammatory_markers_env <- inflammatory_markers_env[inflammatory_markers_env %in% names(mesa_air_sens)]

env_immune_data <- mesa_air_sens %>%
  dplyr::mutate(
    pm25_e1 = dplyr::coalesce(pm25_bl, pm25_ugm3_1_yr_exam),
    no2_e1 = dplyr::coalesce(no2_bl, no2_ppb_1_yr_exam)
  )

env_base_covars <- base_covars[base_covars %in% names(env_immune_data)]
env_covars <- c(env_base_covars, "pm25_e1", "no2_e1")

env_additive_models <- purrr::map_dfr(
  inflammatory_markers_env,
  fit_additive_biomarker_cox,
  data = env_immune_data,
  base_covars = env_covars
) %>%
  dplyr::mutate(
    model = "Air-adjusted (PM2.5 + NO2)",
    biomarker_label = dplyr::recode(biomarker, !!!BIOMARKER_LABELS, .default = biomarker)
  )

geo_keep <- infer_geocode_restriction(env_immune_data$geocode_accuracy)
env_geo_restricted_data <- env_immune_data %>%
  dplyr::filter(geo_keep)

env_geo_additive_models <- purrr::map_dfr(
  inflammatory_markers_env,
  fit_additive_biomarker_cox,
  data = env_geo_restricted_data,
  base_covars = env_covars
) %>%
  dplyr::mutate(
    model = "Air-adjusted + geocode restriction",
    biomarker_label = dplyr::recode(biomarker, !!!BIOMARKER_LABELS, .default = biomarker)
  )

cc_reference_models <- additive_biomarker_models %>%
  dplyr::filter(.data$biomarker %in% inflammatory_markers_env) %>%
  dplyr::mutate(
    model = "Main complete-case",
    biomarker_label = dplyr::recode(biomarker, !!!BIOMARKER_LABELS, .default = biomarker)
  )

environmental_robustness_models <- dplyr::bind_rows(
  cc_reference_models,
  env_additive_models,
  env_geo_additive_models
) %>%
  dplyr::select(model, biomarker, biomarker_label, n, events, hr, lcl, ucl, p_value) %>%
  dplyr::arrange(biomarker_label, model)

environmental_robustness_compact <- environmental_robustness_models %>%
  dplyr::transmute(
    biomarker = biomarker_label,
    model,
    n,
    events,
    `HR (95% CI)` = glue::glue(
      "{formatC(hr, format = 'f', digits = 2)} ({formatC(lcl, format = 'f', digits = 2)}, {formatC(ucl, format = 'f', digits = 2)})"
    ),
    `Wald P` = formatC(p_value, format = "g", digits = 3)
  )

geocode_restriction_summary <- tibble::tibble(
  n_air_subset = nrow(env_immune_data),
  n_geocode_restricted = nrow(env_geo_restricted_data),
  pct_kept = 100 * n_geocode_restricted / n_air_subset,
  rule = dplyr::case_when(
    is.numeric(env_immune_data$geocode_accuracy) ~ "Numeric geocode accuracy restricted to high-accuracy range (<=2 for coded scales, else <= median)",
    TRUE ~ "Text geocode accuracy restricted to labels containing rooftop/parcel/exact/address"
  )
)

readr::write_csv(
  environmental_robustness_models,
  fs::path(OUT_DIR, "mesa_hf_environmental_robustness_models.csv")
)
readr::write_csv(
  environmental_robustness_compact,
  fs::path(OUT_DIR, "mesa_hf_environmental_robustness_models_compact.csv")
)
readr::write_csv(
  geocode_restriction_summary,
  fs::path(OUT_DIR, "mesa_hf_geocode_restriction_summary.csv")
)

environmental_robustness <- environmental_robustness_models
environmental_robustness_compact_table <- environmental_robustness_compact

print(geocode_restriction_summary)
print(environmental_robustness_compact)

# 21) Separate Exam 4-forward analysis (biomarkers + immune phenotypes) -----
immune_cell_vars <- c(
  "abbaspha4", "abeospha4", "ablympha4", "abmoncya4", "abneupha4",
  "epcbla4", "epccd8la4", "epccd8ma4", "epccd8na4", "epcwba4"
)
immune_cell_vars <- immune_cell_vars[immune_cell_vars %in% names(mesa_exam4_landmark)]

exam4_core_markers <- c("z_log_crp4", "z_fib4", "z_log_ddimer4", "z_log_icam4")
exam4_core_markers <- exam4_core_markers[exam4_core_markers %in% names(mesa_exam4_landmark)]

exam4_immune_data <- mesa_exam4_landmark %>%
  dplyr::mutate(
    hf_event_e4 = if_else(!is.na(chf) & chf == 1, 1L, 0L),
    hf_time_days_e4 = dplyr::case_when(
      hf_event_e4 == 1L ~ chftt,
      TRUE ~ fuptt
    )
  ) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::any_of(immune_cell_vars),
      make_z,
      .names = "z_{.col}"
    )
  ) %>%
  dplyr::filter(
    !is.na(race_eth),
    !is.na(hf_time_days_e4),
    hf_time_days_e4 > 0,
    is.na(exall) | exall == 0
  )

immune_marker_z <- paste0("z_", immune_cell_vars)
immune_marker_z <- immune_marker_z[immune_marker_z %in% names(exam4_immune_data)]

exam4_all_markers <- c(exam4_core_markers, immune_marker_z)

exam4_base_covars <- c(
  "age1c", "sex", "race_eth", "site_factor",
  "bmi1c", "sbp1c", "htn1c", "dm031c", "cepgfr1c",
  "educ1", "income1", "cig1c", "pkyrs1c", "lvef_baseline"
)
exam4_base_covars <- exam4_base_covars[
  exam4_base_covars %in% names(exam4_immune_data) &
    vapply(exam4_immune_data[exam4_base_covars], function(x) sum(!is.na(x)) > 0, logical(1))
]

immune_labels <- setNames(
  paste0(stringr::str_to_upper(immune_cell_vars), " (z)"),
  paste0("z_", immune_cell_vars)
)

exam4_marker_labels <- c(
  "z_log_crp4" = "CRP4 (log-z)",
  "z_fib4" = "FIB4 (z)",
  "z_log_ddimer4" = "DDIMER4 (log-z)",
  "z_log_icam4" = "ICAM4 (log-z)",
  immune_labels
)

exam4_additive_models <- purrr::map_dfr(
  exam4_all_markers,
  fit_additive_biomarker_cox,
  data = exam4_immune_data,
  base_covars = exam4_base_covars,
  time_var = "hf_time_days_e4",
  event_var = "hf_event_e4",
  race_var = "race_eth"
) %>%
  dplyr::mutate(
    biomarker_label = dplyr::recode(biomarker, !!!exam4_marker_labels, .default = biomarker),
    analysis = "Exam4-forward separate analysis"
  ) %>%
  dplyr::arrange(biomarker_label)

exam4_interaction_models <- purrr::map_dfr(
  exam4_all_markers,
  fit_interaction_cox,
  data = exam4_immune_data,
  base_covars = exam4_base_covars,
  time_var = "hf_time_days_e4",
  event_var = "hf_event_e4",
  race_var = "race_eth"
) %>%
  dplyr::mutate(
    family = "exam4_biomarker_immune",
    p_adj_holm = p.adjust(p_interaction, method = "holm"),
    q_adj_bh = p.adjust(p_interaction, method = "BH"),
    biomarker_label = dplyr::recode(biomarker, !!!exam4_marker_labels, .default = biomarker),
    analysis = "Exam4-forward separate analysis"
  ) %>%
  dplyr::arrange(biomarker_label)

exam4_additive_compact <- exam4_additive_models %>%
  dplyr::transmute(
    biomarker = biomarker_label,
    n,
    events,
    `HR (95% CI)` = glue::glue(
      "{formatC(hr, format = 'f', digits = 2)} ({formatC(lcl, format = 'f', digits = 2)}, {formatC(ucl, format = 'f', digits = 2)})"
    ),
    `Wald P` = formatC(p_value, format = "g", digits = 3)
  )

exam4_interaction_compact <- exam4_interaction_models %>%
  dplyr::transmute(
    biomarker = biomarker_label,
    n,
    events,
    `Ref HR (95% CI)` = glue::glue(
      "{formatC(hr_ref, format = 'f', digits = 2)} ({formatC(lcl_ref, format = 'f', digits = 2)}, {formatC(ucl_ref, format = 'f', digits = 2)})"
    ),
    `Ref Wald P` = formatC(p_ref, format = "g", digits = 3),
    `Interaction P` = formatC(p_interaction, format = "g", digits = 3),
    `Holm P` = formatC(p_adj_holm, format = "g", digits = 3),
    `BH q` = formatC(q_adj_bh, format = "g", digits = 3)
  )

readr::write_csv(
  exam4_additive_models,
  fs::path(OUT_DIR, "mesa_hf_exam4_forward_additive_models.csv")
)
readr::write_csv(
  exam4_interaction_models,
  fs::path(OUT_DIR, "mesa_hf_exam4_forward_interaction_models.csv")
)
readr::write_csv(
  exam4_additive_compact,
  fs::path(OUT_DIR, "mesa_hf_exam4_forward_additive_models_compact.csv")
)
readr::write_csv(
  exam4_interaction_compact,
  fs::path(OUT_DIR, "mesa_hf_exam4_forward_interaction_models_compact.csv")
)

exam4_forward_additive_models <- exam4_additive_models
exam4_forward_interaction_models <- exam4_interaction_models

print(exam4_additive_compact)
print(exam4_interaction_compact)

# 22) Black-White HF gap attenuation by biomarker elevation ------------------
estimate_black_white_gap <- function(data, time_var, event_var, covars,
                                     biomarkers = character(0), model_label = "") {
  vars_needed <- unique(c(time_var, event_var, "race_bw", covars, biomarkers))

  d <- data %>%
    dplyr::select(dplyr::any_of(vars_needed)) %>%
    dplyr::filter(stats::complete.cases(.), .data[[time_var]] > 0)

  if (nrow(d) == 0L || sum(d[[event_var]], na.rm = TRUE) == 0L) {
    return(tibble::tibble(
      model = model_label,
      n = nrow(d),
      events = sum(d[[event_var]], na.rm = TRUE),
      log_hr_black_vs_white = NA_real_,
      hr_black_vs_white = NA_real_,
      lcl = NA_real_,
      ucl = NA_real_,
      p_value = NA_real_
    ))
  }

  rhs <- paste(c("race_bw", covars, biomarkers), collapse = " + ")
  fml <- stats::as.formula(
    paste0("survival::Surv(", time_var, ", ", event_var, ") ~ ", rhs)
  )

  fit <- tryCatch(
    survival::coxph(fml, data = d, ties = "efron"),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(tibble::tibble(
      model = model_label,
      n = nrow(d),
      events = sum(d[[event_var]], na.rm = TRUE),
      log_hr_black_vs_white = NA_real_,
      hr_black_vs_white = NA_real_,
      lcl = NA_real_,
      ucl = NA_real_,
      p_value = NA_real_
    ))
  }

  td_log <- broom::tidy(fit, exponentiate = FALSE, conf.int = TRUE)
  td_hr <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)

  tr_log <- td_log %>% dplyr::filter(.data$term == "race_bwBlack")
  tr_hr <- td_hr %>% dplyr::filter(.data$term == "race_bwBlack")

  tibble::tibble(
    model = model_label,
    n = nrow(d),
    events = sum(d[[event_var]], na.rm = TRUE),
    log_hr_black_vs_white = if (nrow(tr_log) == 1) tr_log$estimate else NA_real_,
    hr_black_vs_white = if (nrow(tr_hr) == 1) tr_hr$estimate else NA_real_,
    lcl = if (nrow(tr_hr) == 1) tr_hr$conf.low else NA_real_,
    ucl = if (nrow(tr_hr) == 1) tr_hr$conf.high else NA_real_,
    p_value = if (nrow(tr_hr) == 1) tr_hr$p.value else NA_real_
  )
}

summarize_gap_attenuation <- function(tab) {
  base_log_hr <- tab$log_hr_black_vs_white[match("Clinical only", tab$model)]

  tab %>%
    dplyr::mutate(
      pct_excess_risk_explained = 100 * (base_log_hr - log_hr_black_vs_white) / base_log_hr,
      hr_ci = glue::glue(
        "{formatC(hr_black_vs_white, format = 'f', digits = 2)} ({formatC(lcl, format = 'f', digits = 2)}, {formatC(ucl, format = 'f', digits = 2)})"
      ),
      p_value_fmt = formatC(p_value, format = "g", digits = 3)
    )
}

# Exam 1: baseline inflammatory/hemostatic markers
exam1_gap_data <- mesa_main %>%
  dplyr::mutate(
    race_label = dplyr::recode(as.character(race_eth), !!!RACE_LABELS, .default = as.character(race_eth))
  ) %>%
  dplyr::filter(race_label %in% c("White", "Black")) %>%
  dplyr::mutate(race_bw = factor(race_label, levels = c("White", "Black")))

exam1_covars_gap <- base_covars[base_covars %in% names(exam1_gap_data)]
exam1_markers_gap <- c("z_log_crp1", "z_log_il61", "z_fib1", "z_log_ddimer1")
exam1_markers_gap <- exam1_markers_gap[exam1_markers_gap %in% names(exam1_gap_data)]

exam1_gap_models <- dplyr::bind_rows(
  estimate_black_white_gap(
    data = exam1_gap_data,
    time_var = "hf_time_days",
    event_var = "hf_event",
    covars = exam1_covars_gap,
    biomarkers = character(0),
    model_label = "Clinical only"
  ),
  purrr::map_dfr(exam1_markers_gap, function(bm) {
    estimate_black_white_gap(
      data = exam1_gap_data,
      time_var = "hf_time_days",
      event_var = "hf_event",
      covars = exam1_covars_gap,
      biomarkers = bm,
      model_label = paste0("Clinical + ", dplyr::recode(bm, !!!BIOMARKER_LABELS, .default = bm))
    )
  }),
  estimate_black_white_gap(
    data = exam1_gap_data,
    time_var = "hf_time_days",
    event_var = "hf_event",
    covars = exam1_covars_gap,
    biomarkers = exam1_markers_gap,
    model_label = "Clinical + all Exam 1 markers"
  )
) %>%
  summarize_gap_attenuation() %>%
  dplyr::mutate(analysis = "Exam 1")

# Exam 4-forward: core + immune phenotypes as separate analysis
exam4_gap_data <- exam4_immune_data %>%
  dplyr::mutate(
    race_label = dplyr::recode(as.character(race_eth), !!!RACE_LABELS, .default = as.character(race_eth))
  ) %>%
  dplyr::filter(race_label %in% c("White", "Black")) %>%
  dplyr::mutate(race_bw = factor(race_label, levels = c("White", "Black")))

exam4_covars_gap <- exam4_base_covars[exam4_base_covars %in% names(exam4_gap_data) & exam4_base_covars != "race_eth"]
exam4_core_gap <- c("z_log_crp4", "z_fib4", "z_log_ddimer4", "z_log_icam4")
exam4_core_gap <- exam4_core_gap[exam4_core_gap %in% names(exam4_gap_data)]
exam4_immune_gap <- immune_marker_z[immune_marker_z %in% names(exam4_gap_data)]

exam4_gap_models <- dplyr::bind_rows(
  estimate_black_white_gap(
    data = exam4_gap_data,
    time_var = "hf_time_days_e4",
    event_var = "hf_event_e4",
    covars = exam4_covars_gap,
    biomarkers = character(0),
    model_label = "Clinical only"
  ),
  purrr::map_dfr(exam4_core_gap, function(bm) {
    estimate_black_white_gap(
      data = exam4_gap_data,
      time_var = "hf_time_days_e4",
      event_var = "hf_event_e4",
      covars = exam4_covars_gap,
      biomarkers = bm,
      model_label = paste0("Clinical + ", dplyr::recode(bm, !!!exam4_marker_labels, .default = bm))
    )
  }),
  estimate_black_white_gap(
    data = exam4_gap_data,
    time_var = "hf_time_days_e4",
    event_var = "hf_event_e4",
    covars = exam4_covars_gap,
    biomarkers = exam4_core_gap,
    model_label = "Clinical + all Exam 4 core markers"
  ),
  estimate_black_white_gap(
    data = exam4_gap_data,
    time_var = "hf_time_days_e4",
    event_var = "hf_event_e4",
    covars = exam4_covars_gap,
    biomarkers = exam4_immune_gap,
    model_label = "Clinical + all immune phenotypes"
  ),
  estimate_black_white_gap(
    data = exam4_gap_data,
    time_var = "hf_time_days_e4",
    event_var = "hf_event_e4",
    covars = exam4_covars_gap,
    biomarkers = unique(c(exam4_core_gap, exam4_immune_gap)),
    model_label = "Clinical + core + immune phenotypes"
  )
) %>%
  summarize_gap_attenuation() %>%
  dplyr::mutate(analysis = "Exam 4 forward")

black_white_gap_mediation <- dplyr::bind_rows(exam1_gap_models, exam4_gap_models)

black_white_gap_mediation_compact <- black_white_gap_mediation %>%
  dplyr::transmute(
    analysis,
    model,
    n,
    events,
    `Black vs White HR (95% CI)` = hr_ci,
    `Wald P` = p_value_fmt,
    `% excess risk explained vs clinical-only` = formatC(pct_excess_risk_explained, format = "f", digits = 1)
  )

readr::write_csv(
  black_white_gap_mediation,
  fs::path(OUT_DIR, "mesa_hf_black_white_gap_mediation.csv")
)
readr::write_csv(
  black_white_gap_mediation_compact,
  fs::path(OUT_DIR, "mesa_hf_black_white_gap_mediation_compact.csv")
)

print(black_white_gap_mediation_compact)

# 23) Decompose significant Exam 4 race × biomarker interactions -------------
exam4_priority_markers <- c("z_log_ddimer4", "z_fib4")

exam4_sig_interactions <- exam4_interaction_models %>%
  dplyr::filter(
    .data$biomarker %in% exam4_priority_markers |
      (!is.na(.data$p_adj_holm) & .data$p_adj_holm < 0.05) |
      (!is.na(.data$q_adj_bh) & .data$q_adj_bh < 0.05)
  ) %>%
  dplyr::arrange(.data$p_interaction)

exam4_sig_markers <- unique(exam4_sig_interactions$biomarker)

exam4_group_specific_effects <- purrr::map_dfr(
  exam4_sig_markers,
  extract_group_slopes,
  data = exam4_immune_data,
  base_covars = exam4_base_covars,
  time_var = "hf_time_days_e4",
  event_var = "hf_event_e4",
  race_var = "race_eth"
) %>%
  dplyr::left_join(
    exam4_sig_interactions %>%
      dplyr::select(
        biomarker,
        p_interaction,
        p_adj_holm,
        q_adj_bh
      ),
    by = "biomarker"
  ) %>%
  dplyr::mutate(
    biomarker_label = dplyr::recode(biomarker, !!!exam4_marker_labels, .default = biomarker),
    race_label = dplyr::recode(as.character(race_level), !!!RACE_LABELS, .default = as.character(race_level)),
    hr_ci = glue::glue(
      "{formatC(hr, format = 'f', digits = 2)} ({formatC(lcl, format = 'f', digits = 2)}, {formatC(ucl, format = 'f', digits = 2)})"
    ),
    wald_p_fmt = formatC(wald_p, format = "g", digits = 3),
    contrast_p_vs_ref_fmt = dplyr::if_else(
      is.na(contrast_p_vs_ref),
      NA_character_,
      formatC(contrast_p_vs_ref, format = "g", digits = 3)
    ),
    interaction_p_fmt = formatC(p_interaction, format = "g", digits = 3),
    holm_p_fmt = formatC(p_adj_holm, format = "g", digits = 3),
    bh_q_fmt = formatC(q_adj_bh, format = "g", digits = 3)
  ) %>%
  dplyr::arrange(biomarker_label, race_level)

exam4_group_specific_effects_compact <- exam4_group_specific_effects %>%
  dplyr::transmute(
    biomarker = biomarker_label,
    race = race_label,
    n,
    events,
    `Group-specific HR (95% CI)` = hr_ci,
    `Group Wald P` = wald_p_fmt,
    `Contrast P vs White` = contrast_p_vs_ref_fmt,
    `Global interaction P` = interaction_p_fmt,
    `Holm P` = holm_p_fmt,
    `BH q` = bh_q_fmt
  )

exam4_ddimer4_fib4_group_effects <- exam4_group_specific_effects %>%
  dplyr::filter(.data$biomarker %in% c("z_log_ddimer4", "z_fib4"))

readr::write_csv(
  exam4_group_specific_effects,
  fs::path(OUT_DIR, "mesa_hf_exam4_significant_interactions_group_specific.csv")
)
readr::write_csv(
  exam4_group_specific_effects_compact,
  fs::path(OUT_DIR, "mesa_hf_exam4_significant_interactions_group_specific_compact.csv")
)
readr::write_csv(
  exam4_ddimer4_fib4_group_effects,
  fs::path(OUT_DIR, "mesa_hf_exam4_ddimer4_fib4_group_specific.csv")
)

print(exam4_group_specific_effects_compact)

# 24) Decompose significant Exam 1 race × biomarker interactions -------------
exam1_priority_markers <- c("z_log_il61", "z_log_ddimer1")

exam1_interaction_pool <- interaction_results %>%
  dplyr::filter(.data$biomarker %in% primary_biomarkers)

exam1_sig_interactions <- exam1_interaction_pool %>%
  dplyr::filter(
    .data$biomarker %in% exam1_priority_markers |
      (!is.na(.data$p_adj_holm) & .data$p_adj_holm < 0.05) |
      (!is.na(.data$q_adj_bh) & .data$q_adj_bh < 0.05)
  ) %>%
  dplyr::arrange(.data$p_interaction)

exam1_sig_markers <- unique(exam1_sig_interactions$biomarker)

exam1_group_specific_effects <- purrr::map_dfr(
  exam1_sig_markers,
  extract_group_slopes,
  data = mesa_main,
  base_covars = base_covars,
  time_var = "hf_time_days",
  event_var = "hf_event",
  race_var = "race_eth"
) %>%
  dplyr::left_join(
    exam1_sig_interactions %>%
      dplyr::select(
        biomarker,
        p_interaction,
        p_adj_holm,
        q_adj_bh
      ),
    by = "biomarker"
  ) %>%
  dplyr::mutate(
    biomarker_label = dplyr::recode(biomarker, !!!BIOMARKER_LABELS, .default = biomarker),
    race_label = dplyr::recode(as.character(race_level), !!!RACE_LABELS, .default = as.character(race_level)),
    hr_ci = glue::glue(
      "{formatC(hr, format = 'f', digits = 2)} ({formatC(lcl, format = 'f', digits = 2)}, {formatC(ucl, format = 'f', digits = 2)})"
    ),
    wald_p_fmt = formatC(wald_p, format = "g", digits = 3),
    contrast_p_vs_ref_fmt = dplyr::if_else(
      is.na(contrast_p_vs_ref),
      NA_character_,
      formatC(contrast_p_vs_ref, format = "g", digits = 3)
    ),
    interaction_p_fmt = formatC(p_interaction, format = "g", digits = 3),
    holm_p_fmt = formatC(p_adj_holm, format = "g", digits = 3),
    bh_q_fmt = formatC(q_adj_bh, format = "g", digits = 3)
  ) %>%
  dplyr::arrange(biomarker_label, race_level)

exam1_group_specific_effects_compact <- exam1_group_specific_effects %>%
  dplyr::transmute(
    biomarker = biomarker_label,
    race = race_label,
    n,
    events,
    `Group-specific HR (95% CI)` = hr_ci,
    `Group Wald P` = wald_p_fmt,
    `Contrast P vs White` = contrast_p_vs_ref_fmt,
    `Global interaction P` = interaction_p_fmt,
    `Holm P` = holm_p_fmt,
    `BH q` = bh_q_fmt
  )

exam1_il6_ddimer_group_effects <- exam1_group_specific_effects %>%
  dplyr::filter(.data$biomarker %in% c("z_log_il61", "z_log_ddimer1"))

readr::write_csv(
  exam1_group_specific_effects,
  fs::path(OUT_DIR, "mesa_hf_exam1_significant_interactions_group_specific.csv")
)
readr::write_csv(
  exam1_group_specific_effects_compact,
  fs::path(OUT_DIR, "mesa_hf_exam1_significant_interactions_group_specific_compact.csv")
)
readr::write_csv(
  exam1_il6_ddimer_group_effects,
  fs::path(OUT_DIR, "mesa_hf_exam1_il6_ddimer_group_specific.csv")
)

print(exam1_group_specific_effects_compact)

# 25) Fibrinogen × metabolic comorbidity interaction (Exam 1 & Exam 4) ------
prep_metabolic_groups <- function(data, dm_var = "dm031c", bmi_var = "bmi1c") {
  data %>%
    dplyr::mutate(
      dm_num = suppressWarnings(as.numeric(as.character(.data[[dm_var]]))),
      diabetes_status = dplyr::case_when(
        !is.na(dm_num) & dm_num == 1 ~ "Diabetes",
        !is.na(dm_num) & dm_num != 1 ~ "No diabetes",
        TRUE ~ NA_character_
      ),
      obesity_status = dplyr::case_when(
        is.na(.data[[bmi_var]]) ~ NA_character_,
        .data[[bmi_var]] >= 30 ~ "Obesity",
        TRUE ~ "No obesity"
      ),
      metabolic_status = dplyr::case_when(
        diabetes_status == "No diabetes" & obesity_status == "No obesity" ~ "No diabetes / No obesity",
        diabetes_status == "No diabetes" & obesity_status == "Obesity" ~ "No diabetes / Obesity",
        diabetes_status == "Diabetes" & obesity_status == "No obesity" ~ "Diabetes / No obesity",
        diabetes_status == "Diabetes" & obesity_status == "Obesity" ~ "Diabetes / Obesity",
        TRUE ~ NA_character_
      ),
      diabetes_status = factor(diabetes_status, levels = c("No diabetes", "Diabetes")),
      obesity_status = factor(obesity_status, levels = c("No obesity", "Obesity")),
      metabolic_status = factor(
        metabolic_status,
        levels = c(
          "No diabetes / No obesity",
          "No diabetes / Obesity",
          "Diabetes / No obesity",
          "Diabetes / Obesity"
        )
      )
    ) %>%
    dplyr::select(-dm_num)
}

run_fibrinogen_metabolic <- function(data, fibrinogen_var, covars, time_var, event_var, analysis_label,
                                     dm_var = "dm031c", bmi_var = "bmi1c") {
  d <- prep_metabolic_groups(data, dm_var = dm_var, bmi_var = bmi_var)

  covars_adj <- setdiff(covars, c(dm_var, bmi_var))

  dm_test <- fit_interaction_cox(
    data = d,
    biomarker = fibrinogen_var,
    base_covars = covars_adj,
    time_var = time_var,
    event_var = event_var,
    race_var = "diabetes_status"
  ) %>%
    dplyr::mutate(interaction_type = "Fibrinogen × diabetes")

  obesity_test <- fit_interaction_cox(
    data = d,
    biomarker = fibrinogen_var,
    base_covars = covars_adj,
    time_var = time_var,
    event_var = event_var,
    race_var = "obesity_status"
  ) %>%
    dplyr::mutate(interaction_type = "Fibrinogen × obesity")

  metabolic_test <- fit_interaction_cox(
    data = d,
    biomarker = fibrinogen_var,
    base_covars = covars_adj,
    time_var = time_var,
    event_var = event_var,
    race_var = "metabolic_status"
  ) %>%
    dplyr::mutate(interaction_type = "Fibrinogen × metabolic phenotype (4-level)")

  tests <- dplyr::bind_rows(dm_test, obesity_test, metabolic_test) %>%
    dplyr::mutate(
      analysis = analysis_label,
      p_adj_holm = p.adjust(p_interaction, method = "holm"),
      q_adj_bh = p.adjust(p_interaction, method = "BH"),
      p_interaction_fmt = formatC(p_interaction, format = "g", digits = 3),
      holm_p_fmt = formatC(p_adj_holm, format = "g", digits = 3),
      bh_q_fmt = formatC(q_adj_bh, format = "g", digits = 3)
    )

  metabolic_slopes <- extract_group_slopes(
    data = d,
    biomarker = fibrinogen_var,
    base_covars = covars_adj,
    time_var = time_var,
    event_var = event_var,
    race_var = "metabolic_status"
  ) %>%
    dplyr::left_join(
      tests %>%
        dplyr::filter(interaction_type == "Fibrinogen × metabolic phenotype (4-level)") %>%
        dplyr::select(p_interaction, p_adj_holm, q_adj_bh),
      by = character()
    ) %>%
    dplyr::mutate(
      analysis = analysis_label,
      fibrinogen = fibrinogen_var,
      subgroup = as.character(race_level),
      hr_ci = glue::glue(
        "{formatC(hr, format = 'f', digits = 2)} ({formatC(lcl, format = 'f', digits = 2)}, {formatC(ucl, format = 'f', digits = 2)})"
      ),
      wald_p_fmt = formatC(wald_p, format = "g", digits = 3),
      contrast_p_vs_ref_fmt = dplyr::if_else(
        is.na(contrast_p_vs_ref),
        NA_character_,
        formatC(contrast_p_vs_ref, format = "g", digits = 3)
      ),
      interaction_p_fmt = formatC(p_interaction, format = "g", digits = 3),
      holm_p_fmt = formatC(p_adj_holm, format = "g", digits = 3),
      bh_q_fmt = formatC(q_adj_bh, format = "g", digits = 3)
    )

  list(tests = tests, slopes = metabolic_slopes)
}

exam1_fib_metabolic <- run_fibrinogen_metabolic(
  data = mesa_main,
  fibrinogen_var = "z_fib1",
  covars = base_covars,
  time_var = "hf_time_days",
  event_var = "hf_event",
  analysis_label = "Exam 1"
)

exam4_fib_metabolic <- run_fibrinogen_metabolic(
  data = exam4_immune_data,
  fibrinogen_var = "z_fib4",
  covars = exam4_base_covars,
  time_var = "hf_time_days_e4",
  event_var = "hf_event_e4",
  analysis_label = "Exam 4 forward"
)

fibrinogen_metabolic_interaction_tests <- dplyr::bind_rows(
  exam1_fib_metabolic$tests,
  exam4_fib_metabolic$tests
) %>%
  dplyr::select(
    analysis,
    interaction_type,
    n,
    events,
    p_interaction,
    p_adj_holm,
    q_adj_bh,
    p_interaction_fmt,
    holm_p_fmt,
    bh_q_fmt
  )

fibrinogen_metabolic_group_slopes <- dplyr::bind_rows(
  exam1_fib_metabolic$slopes,
  exam4_fib_metabolic$slopes
) %>%
  dplyr::select(
    analysis,
    fibrinogen,
    subgroup,
    n,
    events,
    hr,
    lcl,
    ucl,
    wald_p,
    contrast_p_vs_ref,
    p_interaction,
    p_adj_holm,
    q_adj_bh,
    hr_ci,
    wald_p_fmt,
    contrast_p_vs_ref_fmt,
    interaction_p_fmt,
    holm_p_fmt,
    bh_q_fmt
  )

fibrinogen_metabolic_group_slopes_compact <- fibrinogen_metabolic_group_slopes %>%
  dplyr::transmute(
    analysis,
    fibrinogen,
    subgroup,
    n,
    events,
    `Group-specific HR (95% CI)` = hr_ci,
    `Group Wald P` = wald_p_fmt,
    `Contrast P vs reference` = contrast_p_vs_ref_fmt,
    `Global interaction P` = interaction_p_fmt,
    `Holm P` = holm_p_fmt,
    `BH q` = bh_q_fmt
  )

readr::write_csv(
  fibrinogen_metabolic_interaction_tests,
  fs::path(OUT_DIR, "mesa_hf_fibrinogen_metabolic_interaction_tests.csv")
)
readr::write_csv(
  fibrinogen_metabolic_group_slopes,
  fs::path(OUT_DIR, "mesa_hf_fibrinogen_metabolic_group_slopes.csv")
)
readr::write_csv(
  fibrinogen_metabolic_group_slopes_compact,
  fs::path(OUT_DIR, "mesa_hf_fibrinogen_metabolic_group_slopes_compact.csv")
)

print(fibrinogen_metabolic_interaction_tests)
print(fibrinogen_metabolic_group_slopes_compact)

# 26) Composite inflammatory burden score ------------------------------------
build_pca_burden_score <- function(data, markers, analysis_label, score_name) {
  markers <- markers[markers %in% names(data)]

  if (length(markers) < 2L) {
    return(list(
      score = rep(NA_real_, nrow(data)),
      loadings = tibble::tibble(
        analysis = analysis_label,
        score_name = score_name,
        marker = markers,
        loading_pc1 = NA_real_,
        abs_loading = NA_real_,
        variance_explained_pc1 = NA_real_,
        n_complete = 0L
      )
    ))
  }

  cc <- stats::complete.cases(data[, markers, drop = FALSE])
  n_complete <- sum(cc)

  if (n_complete < max(25L, length(markers) + 5L)) {
    return(list(
      score = rep(NA_real_, nrow(data)),
      loadings = tibble::tibble(
        analysis = analysis_label,
        score_name = score_name,
        marker = markers,
        loading_pc1 = NA_real_,
        abs_loading = NA_real_,
        variance_explained_pc1 = NA_real_,
        n_complete = n_complete
      )
    ))
  }

  pca_fit <- stats::prcomp(data[cc, markers, drop = FALSE], center = TRUE, scale. = TRUE)
  score_raw <- rep(NA_real_, nrow(data))
  score_raw[cc] <- pca_fit$x[, 1]

  mean_marker <- rowMeans(data[cc, markers, drop = FALSE], na.rm = TRUE)
  orient_sign <- if (stats::cor(score_raw[cc], mean_marker, use = "complete.obs") < 0) -1 else 1

  score_final <- make_z(score_raw * orient_sign)
  loadings <- tibble::tibble(
    analysis = analysis_label,
    score_name = score_name,
    marker = rownames(pca_fit$rotation),
    loading_pc1 = orient_sign * pca_fit$rotation[, 1],
    abs_loading = abs(orient_sign * pca_fit$rotation[, 1]),
    variance_explained_pc1 = summary(pca_fit)$importance[2, 1],
    n_complete = n_complete
  ) %>%
    dplyr::arrange(dplyr::desc(abs_loading))

  list(score = score_final, loadings = loadings)
}

exam1_burden_markers <- c("z_log_crp1", "z_log_il61", "z_fib1", "z_log_ddimer1")
exam1_burden <- build_pca_burden_score(
  data = mesa_main,
  markers = exam1_burden_markers,
  analysis_label = "Exam 1",
  score_name = "z_exam1_inflammatory_burden"
)

mesa_main <- mesa_main %>%
  dplyr::mutate(z_exam1_inflammatory_burden = exam1_burden$score)

exam4_burden_markers <- c(
  "z_log_crp4", "z_fib4", "z_log_ddimer4", "z_log_icam4",
  "z_epccd8la4", "z_epcwba4"
)
exam4_burden <- build_pca_burden_score(
  data = exam4_immune_data,
  markers = exam4_burden_markers,
  analysis_label = "Exam 4 forward",
  score_name = "z_exam4_multidomain_inflammatory_burden"
)

exam4_immune_data <- exam4_immune_data %>%
  dplyr::mutate(z_exam4_multidomain_inflammatory_burden = exam4_burden$score)

composite_burden_loadings <- dplyr::bind_rows(
  exam1_burden$loadings,
  exam4_burden$loadings
) %>%
  dplyr::mutate(
    marker_label = dplyr::case_when(
      marker %in% names(BIOMARKER_LABELS) ~ dplyr::recode(marker, !!!BIOMARKER_LABELS),
      marker %in% names(exam4_marker_labels) ~ dplyr::recode(marker, !!!exam4_marker_labels),
      TRUE ~ marker
    )
  )

composite_burden_additive <- dplyr::bind_rows(
  fit_additive_biomarker_cox(
    data = mesa_main,
    biomarker = "z_exam1_inflammatory_burden",
    base_covars = base_covars,
    time_var = "hf_time_days",
    event_var = "hf_event",
    race_var = "race_eth"
  ) %>% dplyr::mutate(analysis = "Exam 1", score_name = "Exam 1 inflammatory burden"),
  fit_additive_biomarker_cox(
    data = exam4_immune_data,
    biomarker = "z_exam4_multidomain_inflammatory_burden",
    base_covars = exam4_base_covars,
    time_var = "hf_time_days_e4",
    event_var = "hf_event_e4",
    race_var = "race_eth"
  ) %>% dplyr::mutate(analysis = "Exam 4 forward", score_name = "Exam 4 multidomain inflammatory burden")
) %>%
  dplyr::mutate(
    hr_ci = glue::glue(
      "{formatC(hr, format = 'f', digits = 2)} ({formatC(lcl, format = 'f', digits = 2)}, {formatC(ucl, format = 'f', digits = 2)})"
    ),
    p_value_fmt = formatC(p_value, format = "g", digits = 3)
  )

composite_burden_interaction <- dplyr::bind_rows(
  fit_interaction_cox(
    data = mesa_main,
    biomarker = "z_exam1_inflammatory_burden",
    base_covars = base_covars,
    time_var = "hf_time_days",
    event_var = "hf_event",
    race_var = "race_eth"
  ) %>% dplyr::mutate(analysis = "Exam 1", score_name = "Exam 1 inflammatory burden"),
  fit_interaction_cox(
    data = exam4_immune_data,
    biomarker = "z_exam4_multidomain_inflammatory_burden",
    base_covars = exam4_base_covars,
    time_var = "hf_time_days_e4",
    event_var = "hf_event_e4",
    race_var = "race_eth"
  ) %>% dplyr::mutate(analysis = "Exam 4 forward", score_name = "Exam 4 multidomain inflammatory burden")
) %>%
  dplyr::mutate(
    p_adj_holm = p.adjust(p_interaction, method = "holm"),
    q_adj_bh = p.adjust(p_interaction, method = "BH"),
    hr_ref_ci = glue::glue(
      "{formatC(hr_ref, format = 'f', digits = 2)} ({formatC(lcl_ref, format = 'f', digits = 2)}, {formatC(ucl_ref, format = 'f', digits = 2)})"
    ),
    p_ref_fmt = formatC(p_ref, format = "g", digits = 3),
    p_interaction_fmt = formatC(p_interaction, format = "g", digits = 3),
    holm_p_fmt = formatC(p_adj_holm, format = "g", digits = 3),
    bh_q_fmt = formatC(q_adj_bh, format = "g", digits = 3)
  )

composite_burden_additive_compact <- composite_burden_additive %>%
  dplyr::transmute(
    analysis,
    score = score_name,
    n,
    events,
    `HR (95% CI)` = hr_ci,
    `Wald P` = p_value_fmt
  )

composite_burden_interaction_compact <- composite_burden_interaction %>%
  dplyr::transmute(
    analysis,
    score = score_name,
    n,
    events,
    `Ref HR (95% CI)` = hr_ref_ci,
    `Ref Wald P` = p_ref_fmt,
    `Interaction P` = p_interaction_fmt,
    `Holm P` = holm_p_fmt,
    `BH q` = bh_q_fmt
  )

readr::write_csv(
  composite_burden_loadings,
  fs::path(OUT_DIR, "mesa_hf_composite_inflammatory_burden_loadings.csv")
)
readr::write_csv(
  composite_burden_additive,
  fs::path(OUT_DIR, "mesa_hf_composite_inflammatory_burden_additive.csv")
)
readr::write_csv(
  composite_burden_interaction,
  fs::path(OUT_DIR, "mesa_hf_composite_inflammatory_burden_interaction.csv")
)

print(composite_burden_additive_compact)
print(composite_burden_interaction_compact)
