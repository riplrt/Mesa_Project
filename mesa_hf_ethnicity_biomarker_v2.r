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
#   mesaaire1_drepos_20240603.csv and applies them before coalescing.
# - Added crp1_source/fib1_source provenance columns to mesa_main.
#------------------------------------------------------------------------------

# 0) Packages -----------------------------------------------------------------
packages <- c("tidyverse", "janitor", "fs", "glue", "survival", "broom")

to_install <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

invisible(lapply(packages, library, character.only = TRUE))
options(stringsAsFactors = FALSE)

# 1) Paths --------------------------------------------------------------------
BASE_DIR <- "data"
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

evt_secondary <- read_mesa_csv(paths$event_secondary, c("fuptt", "prebase", "exall", "chf", "chftt"))

exam1 <- read_mesa_csv(paths$exam1, c(
  "age1c", "gender1", "race1c", "bmi1c", "sbp1c", "htn1c", "dm031c",
  "chol1", "hdl1", "ldl1", "trig1", "cepgfr1c", "educ1", "income1",
  "cig1c", "pkyrs1c", "crp1", "crp1m", "il61", "fib1", "fib1m", "ddimer1",
  "olvef1"
))

site <- read_mesa_csv(paths$site, c("site1c"))

air_all <- read_mesa_csv(paths$air, c(
  "exam", "pm25_bl", "pm25_fu", "no2_bl", "no2_fu",
  "pm25_ugm3_1_yr_exam", "no2_ppb_1_yr_exam"
))

geocode_all <- read_mesa_csv(paths$geocode, c("exam", "accuracy"))

cardiac <- load_cardiac(paths$cardiac079, paths$cardiac244)

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
if (!all(is.na(aire1$crp1m_air))) aire1 <- aire1 %>% mutate(crp1_air = if_else(crp1m_air == 1, NA_real_, crp1_air))
if (!all(is.na(aire1$fib1m_air))) aire1 <- aire1 %>% mutate(fib1_air = if_else(fib1m_air == 1, NA_real_, fib1_air))

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
if ("crp1m" %in% names(exam1)) exam1 <- exam1 %>% mutate(crp1 = if_else(crp1m == 1, NA_real_, crp1))
if ("fib1m" %in% names(exam1)) exam1 <- exam1 %>% mutate(fib1 = if_else(fib1m == 1, NA_real_, fib1))

if ("crp4m" %in% names(exam4_biom)) exam4_biom <- exam4_biom %>% mutate(crp4 = if_else(crp4m == 1, NA_real_, crp4))
if ("fib4m" %in% names(exam4_biom)) exam4_biom <- exam4_biom %>% mutate(fib4 = if_else(fib4m == 1, NA_real_, fib4))

# 5) Primary Exam 1 analytic dataset ------------------------------------------
mesa_main <- exam1 %>%
  left_join(aire1, by = "mesaid") %>%
  left_join(site, by = "mesaid") %>%
  left_join(evt_primary, by = "mesaid") %>%
  left_join(evt_secondary, by = "mesaid") %>%
  left_join(cardiac, by = "mesaid") %>%
  mutate(
    crp1_source = dplyr::case_when(
      !is.na(crp1)     ~ "exam1",
      !is.na(crp1_air) ~ "aire1",
      TRUE             ~ NA_character_
    ),
    fib1_source = dplyr::case_when(
      !is.na(fib1)     ~ "exam1",
      !is.na(fib1_air) ~ "aire1",
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
# Re-extract main-effect HRs from the full interaction models for each
# available biomarker so every result is on a single plot.

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