# =============================================================================
# Study 1 — GLMM Analysis: Head Impact Acceleration in Judo Athletes
#
# Run from the PROJECT ROOT:
#   Rscript code/01_glmm_analysis.R
#
# Outputs (saved to output/models/):
#   fits_A1.rds, fits_A2.rds
#   fits_A1_mech.rds
#   int_results_A1_pav.rds, int_results_A1_dpav.rds
#   int_results_A2_pav.rds, int_results_A2_dpav.rds
#   df_A1_int.rds, df_A2_int.rds
#   summary_A1.csv, summary_A2.csv, summary_A1_mech.csv
#   aic_A1.csv, aic_A2.csv
#
# Model formula:
#   outcome ~ skill_level + (1 | team / device_sn) + (1 | team:date)
# Family selected by AIC: Gaussian (identity) / Gaussian (log) / Gamma (log)
# =============================================================================

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(tidyverse)
  library(glmmTMB)
  library(performance)
})

# =============================================================================
# Paths  (relative to project root)
# =============================================================================

DATA_DIR    <- "data"
MODEL_DIR   <- "output/models"
dir.create(MODEL_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Configuration
# =============================================================================

CONFIG <- list(
  A1 = list(file = "child_A1_high_ang_1.csv", label = "A1"),
  A2 = list(file = "child_A2_high_lin.csv",   label = "A2"),
  skill_file            = "skill_level.csv",
  mechanics_all_file    = "master_mechanics_all.csv",
  mechanics_all_A2_file = "master_mechanics_all_A2.csv",
  metrics  = c("PLA_g", "PAA_rad_s2", "PLV_m_s", "PAV_rad_s"),
  labels   = c("PLA (g)", "PAA (rad/s2)", "PLV (m/s)", "PAV (rad/s)"),
  directions      = c("FRONT", "REAR", "LEFT", "RIGHT", "TOP", "BOTTOM"),
  side_directions = c("FRONT", "REAR", "SIDE", "TOP", "BOTTOM")
)

# =============================================================================
# Helper functions
# =============================================================================

read_impact <- function(filename) {
  read_csv(file.path(DATA_DIR, filename), show_col_types = FALSE) |>
    (\(df) if ("situation1" %in% names(df)) select(df, -situation1) else df)() |>
    rename_with(tolower) |>
    rename(
      impact_id  = `impact#`,
      device_sn  = `device sn`,
      PLA_g      = `peak linear acceleration (pla, g)`,
      PAA_rad_s2 = `peak angular acceleration (paa, rad/sec^2)`,
      PLV_m_s    = `peak linear velocity (plv, m/sec)`,
      PAV_rad_s  = `peak angular velocity (pav, rad/sec)`,
      Work_J     = `work (j)`,
      situation  = situation1
    ) |>
    mutate(
      situation = as.character(situation),
      direction = str_extract(location, "^[A-Z]+")
    )
}

prep_data <- function(df_raw, skill_tbl) {
  df_raw |>
    left_join(skill_tbl, by = "device_sn") |>
    filter(!situation %in% c("NG", "N") | is.na(situation)) |>
    mutate(
      skill_level    = factor(skill_level, levels = c("beginner", "experienced")),
      gender         = factor(gender,      levels = c("male", "female")),
      team           = factor(team),
      device_sn      = factor(device_sn),
      date           = as.factor(date),
      direction      = factor(direction, levels = CONFIG$directions),
      direction_side = factor(
        if_else(direction %in% c("LEFT", "RIGHT"), "SIDE", as.character(direction)),
        levels = CONFIG$side_directions
      )
    )
}

# Compare Gaussian(identity), Gaussian(log), Gamma(log) by AIC
compare_families <- function(response, data) {
  d   <- data |> filter(!is.na(.data[[response]]), !is.na(skill_level))
  rhs <- "skill_level + (1 | team / device_sn) + (1 | team:date)"
  fml <- as.formula(paste0(response, " ~ ", rhs))
  results <- list()
  tryCatch({ m <- lmer(fml, data = d, REML = FALSE)
    results[["Gaussian (identity)"]] <- AIC(m) }, error = function(e) NULL)
  tryCatch({ m <- glmer(fml, data = d, family = gaussian(link = "log"))
    results[["Gaussian (log)"]]      <- AIC(m) }, error = function(e) NULL)
  tryCatch({ m <- glmer(fml, data = d, family = Gamma(link = "log"))
    results[["Gamma (log)"]]         <- AIC(m) }, error = function(e) NULL)
  tibble(Model = names(results), AIC = round(unlist(results), 1)) |>
    arrange(AIC) |>
    mutate(Best = if_else(row_number() == 1L, "yes", ""))
}

# Fit the AIC-best model and return list(model, family, aic, response)
fit_best_model <- function(response, data) {
  d          <- data |> filter(!is.na(.data[[response]]), !is.na(skill_level))
  rhs        <- "skill_level + (1 | team / device_sn) + (1 | team:date)"
  fml        <- as.formula(paste0(response, " ~ ", rhs))
  best_aic   <- Inf; best_model <- NULL; best_family <- NULL
  tryCatch({ m <- lmer(fml, data = d, REML = FALSE)
    if (AIC(m) < best_aic) { best_aic <- AIC(m); best_model <- m
      best_family <- "Gaussian (identity)" } }, error = function(e) NULL)
  tryCatch({ m <- glmer(fml, data = d, family = gaussian(link = "log"))
    if (AIC(m) < best_aic) { best_aic <- AIC(m); best_model <- m
      best_family <- "Gaussian (log)" } }, error = function(e) NULL)
  tryCatch({ m <- glmer(fml, data = d, family = Gamma(link = "log"))
    if (AIC(m) < best_aic) { best_aic <- AIC(m); best_model <- m
      best_family <- "Gamma (log)" } }, error = function(e) NULL)
  list(model = best_model, family = best_family, aic = best_aic, response = response)
}

# Extract skill_level fixed effect with 95% CI and exp(estimate)
extract_skill_effect <- function(fit, label) {
  m <- fit$model
  if (is.null(m)) return(NULL)
  coef_tbl <- tryCatch(
    coef(summary(m)) |> as.data.frame() |> rownames_to_column("term"),
    error = function(e) NULL)
  if (is.null(coef_tbl)) return(NULL)
  row <- coef_tbl |> filter(grepl("skill_level", term))
  if (nrow(row) == 0) return(NULL)
  p_col <- intersect(c("Pr(>|t|)", "Pr(>|z|)"), colnames(coef_tbl))
  p_val <- if (length(p_col) > 0) row[[p_col[1]]] else NA_real_
  est   <- row$Estimate; se <- row[["Std. Error"]]
  tibble(
    Metric       = label,
    Model        = fit$family,
    Estimate     = round(est, 4),
    SE           = round(se,  4),
    CI_low       = round(est - 1.96 * se, 4),
    CI_high      = round(est + 1.96 * se, 4),
    Exp_Est      = round(exp(est), 4),
    Exp_CI_low   = round(exp(est - 1.96 * se), 4),
    Exp_CI_high  = round(exp(est + 1.96 * se), 4),
    Direction    = if_else(est > 0, "Exp > Beg", "Exp < Beg"),
    p.value      = round(p_val, 4),
    Sig          = case_when(
      p_val < 0.001 ~ "***", p_val < 0.01 ~ "**",
      p_val < 0.05  ~ "*",   p_val < 0.10 ~ ".",
      TRUE ~ "n.s.")
  )
}

# Interaction GLMM: outcome ~ skill_level * situation_type + nested RE (Gamma, log)
run_interaction_glmm <- function(metric, label, data) {
  d   <- data |> filter(!is.na(.data[[metric]]), .data[[metric]] > 0)
  fml <- as.formula(paste0(
    metric, " ~ skill_level * situation_type",
    " + (1 | team / device_sn) + (1 | team:date)"))
  fit <- tryCatch(
    glmmTMB(fml, data = d, family = Gamma(link = "log")),
    error = function(e) { message(sprintf("[%s] convergence failed: %s", label, e$message)); NULL })
  if (is.null(fit)) return(NULL)
  coef(summary(fit))$cond |> as.data.frame() |> rownames_to_column("Term") |>
    mutate(
      Metric       = label,
      beta         = round(Estimate, 4),
      SE           = round(`Std. Error`, 4),
      CI_low       = round(Estimate - 1.96 * `Std. Error`, 4),
      CI_high      = round(Estimate + 1.96 * `Std. Error`, 4),
      RR           = round(exp(Estimate), 4),
      RR_CI_low    = round(exp(Estimate - 1.96 * `Std. Error`), 4),
      RR_CI_high   = round(exp(Estimate + 1.96 * `Std. Error`), 4),
      p_value      = `Pr(>|z|)`,
      p_display    = if_else(`Pr(>|z|)` < 0.001, "<0.001", sprintf("%.3f", `Pr(>|z|)`)),
      Term = case_when(
        Term == "(Intercept)"                                                     ~ "Intercept",
        grepl("skill_level.*:.*situation_type|situation_type.*:.*skill_level", Term) ~
          "skill x situation (non-BRE)",
        grepl("skill_level",    Term) ~ "Skill (experienced)",
        grepl("situation_type", Term) ~ "Situation (non-BRE)",
        TRUE ~ Term),
      .before = 1
    ) |>
    select(Metric, Term, beta, SE, CI_low, CI_high, RR, RR_CI_low, RR_CI_high, p_value, p_display)
}

# Build BRE / non-BRE combined dataset for interaction models
make_int_df <- function(df) {
  bind_rows(
    df |> filter(situation == "TRE")             |> mutate(situation_type = "BRE"),
    df |> filter(situation %in% c("CIE", "GWE")) |> mutate(situation_type = "non-BRE")
  ) |>
    filter(!is.na(skill_level)) |>
    mutate(
      situation_type = factor(situation_type, levels = c("BRE", "non-BRE")),
      device_sn = factor(device_sn),
      team      = factor(team),
      date      = factor(date)
    )
}

# =============================================================================
# Data loading & preprocessing
# =============================================================================

skill_tbl <- read_csv(file.path(DATA_DIR, CONFIG$skill_file), show_col_types = FALSE)

df_A1 <- read_impact(CONFIG$A1$file) |> prep_data(skill_tbl)
df_A2 <- read_impact(CONFIG$A2$file) |> prep_data(skill_tbl)

# Join biomechanical metrics (DeltaPAV, PeakAlpha, PeakOmega)
join_mech <- function(df, mech_file) {
  mech <- read_csv(file.path(DATA_DIR, mech_file), show_col_types = FALSE) |>
    transmute(
      device_sn_chr = as.character(DeviceSN),
      date_chr      = as.character(Date),
      `local time`  = LocalTime,
      DeltaPAV, PeakAlpha, PeakOmega
    ) |>
    distinct(device_sn_chr, date_chr, `local time`, .keep_all = TRUE)
  df |>
    mutate(device_sn_chr = as.character(device_sn), date_chr = as.character(date)) |>
    left_join(mech, by = c("device_sn_chr", "date_chr", "local time")) |>
    select(-device_sn_chr, -date_chr)
}

df_A1 <- join_mech(df_A1, CONFIG$mechanics_all_file)
df_A2 <- join_mech(df_A2, CONFIG$mechanics_all_A2_file)

cat(sprintf("A1: %d events (DeltaPAV non-NA: %d)\n", nrow(df_A1), sum(!is.na(df_A1$DeltaPAV))))
cat(sprintf("A2: %d events (DeltaPAV non-NA: %d)\n", nrow(df_A2), sum(!is.na(df_A2$DeltaPAV))))

labels_map   <- setNames(CONFIG$labels, CONFIG$metrics)
mech_metrics <- c("DeltaPAV", "PeakAlpha", "PeakOmega")
mech_labels  <- c("DeltaPAV (rad/s)", "PeakAlpha (rad/s2)", "PeakOmega (rad/s)")

# =============================================================================
# AIC model selection
# =============================================================================

message("--- AIC comparison ---")
aic_A1 <- map_dfr(CONFIG$metrics, function(m) {
  compare_families(m, df_A1) |> mutate(Metric = labels_map[m], .before = 1)
})
aic_A2 <- map_dfr(CONFIG$metrics, function(m) {
  compare_families(m, df_A2) |> mutate(Metric = labels_map[m], .before = 1)
})
write_csv(aic_A1, file.path(MODEL_DIR, "aic_A1.csv"))
write_csv(aic_A2, file.path(MODEL_DIR, "aic_A2.csv"))
cat("AIC results saved.\n")

# =============================================================================
# Main GLMM — skill level main effect
# =============================================================================

message("--- Fitting main effect models ---")
fits_A1 <- setNames(map(CONFIG$metrics, fit_best_model, data = df_A1), CONFIG$metrics)
fits_A2 <- setNames(map(CONFIG$metrics, fit_best_model, data = df_A2), CONFIG$metrics)

fits_A1_mech <- setNames(map(mech_metrics, fit_best_model, data = df_A1), mech_metrics)

summary_A1      <- map2_dfr(fits_A1,      CONFIG$labels, extract_skill_effect)
summary_A2      <- map2_dfr(fits_A2,      CONFIG$labels, extract_skill_effect)
summary_A1_mech <- map2_dfr(fits_A1_mech, mech_labels,   extract_skill_effect)

write_csv(summary_A1,      file.path(MODEL_DIR, "summary_A1.csv"))
write_csv(summary_A2,      file.path(MODEL_DIR, "summary_A2.csv"))
write_csv(summary_A1_mech, file.path(MODEL_DIR, "summary_A1_mech.csv"))

saveRDS(fits_A1,      file.path(MODEL_DIR, "fits_A1.rds"))
saveRDS(fits_A2,      file.path(MODEL_DIR, "fits_A2.rds"))
saveRDS(fits_A1_mech, file.path(MODEL_DIR, "fits_A1_mech.rds"))

print(summary_A1)
print(summary_A2)
print(summary_A1_mech)

# =============================================================================
# Interaction GLMM — skill_level x event type (BRE / non-BRE)
# =============================================================================

int_metrics_pav  <- c("PLA_g", "PAA_rad_s2", "PAV_rad_s")
int_labels_pav   <- c("PLA (g)", "PAA (rad/s2)", "PAV (rad/s)")
int_metrics_dpav <- c("PLA_g", "PAA_rad_s2", "DeltaPAV")
int_labels_dpav  <- c("PLA (g)", "PAA (rad/s2)", "DeltaPAV (rad/s)")

df_A1_int <- make_int_df(df_A1)
df_A2_int <- make_int_df(df_A2)

message("--- Fitting interaction models ---")
int_results_A1_pav  <- map2_dfr(int_metrics_pav,  int_labels_pav,  ~ run_interaction_glmm(.x, .y, df_A1_int))
int_results_A1_dpav <- map2_dfr(int_metrics_dpav, int_labels_dpav, ~ run_interaction_glmm(.x, .y, df_A1_int))
int_results_A2_pav  <- map2_dfr(int_metrics_pav,  int_labels_pav,  ~ run_interaction_glmm(.x, .y, df_A2_int))
int_results_A2_dpav <- map2_dfr(int_metrics_dpav, int_labels_dpav, ~ run_interaction_glmm(.x, .y, df_A2_int))

saveRDS(int_results_A1_pav,  file.path(MODEL_DIR, "int_results_A1_pav.rds"))
saveRDS(int_results_A1_dpav, file.path(MODEL_DIR, "int_results_A1_dpav.rds"))
saveRDS(int_results_A2_pav,  file.path(MODEL_DIR, "int_results_A2_pav.rds"))
saveRDS(int_results_A2_dpav, file.path(MODEL_DIR, "int_results_A2_dpav.rds"))
saveRDS(df_A1_int, file.path(MODEL_DIR, "df_A1_int.rds"))
saveRDS(df_A2_int, file.path(MODEL_DIR, "df_A2_int.rds"))

write_csv(int_results_A1_pav,  file.path(MODEL_DIR, "int_results_A1_pav.csv"))
write_csv(int_results_A1_dpav, file.path(MODEL_DIR, "int_results_A1_dpav.csv"))
write_csv(int_results_A2_pav,  file.path(MODEL_DIR, "int_results_A2_pav.csv"))
write_csv(int_results_A2_dpav, file.path(MODEL_DIR, "int_results_A2_dpav.csv"))

print(int_results_A1_pav)
print(int_results_A1_dpav)
print(int_results_A2_pav)
print(int_results_A2_dpav)

message("Done. All models and results saved to: ", MODEL_DIR)
