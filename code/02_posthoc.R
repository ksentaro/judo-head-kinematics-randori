# =============================================================================
# Study 1 — Post-hoc Analyses: emmeans & Pairwise Comparisons
#
# Run AFTER 01_glmm_analysis.R (requires output/models/*.rds).
# Run from the PROJECT ROOT:
#   Rscript code/02_posthoc.R
#
# Contents:
#   1. Main effect models — EMM & pairwise (skill_level)
#   2. Interaction models — interaction contrasts
#   3. Interaction models — EMM & pairwise (skill_level | situation_type)
#   4. ICC by random effect group
# =============================================================================

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(tidyverse)
  library(glmmTMB)
  library(performance)
})

MODEL_DIR  <- "output/models"
OUTPUT_DIR <- "output/tables"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Load fitted models
# =============================================================================

fits_A1      <- readRDS(file.path(MODEL_DIR, "fits_A1.rds"))
fits_A2      <- readRDS(file.path(MODEL_DIR, "fits_A2.rds"))
fits_A1_mech <- readRDS(file.path(MODEL_DIR, "fits_A1_mech.rds"))

int_results_A1_pav  <- readRDS(file.path(MODEL_DIR, "int_results_A1_pav.rds"))
int_results_A1_dpav <- readRDS(file.path(MODEL_DIR, "int_results_A1_dpav.rds"))
int_results_A2_pav  <- readRDS(file.path(MODEL_DIR, "int_results_A2_pav.rds"))
int_results_A2_dpav <- readRDS(file.path(MODEL_DIR, "int_results_A2_dpav.rds"))

df_A1_int <- readRDS(file.path(MODEL_DIR, "df_A1_int.rds"))
df_A2_int <- readRDS(file.path(MODEL_DIR, "df_A2_int.rds"))

# Metric labels for readability
metrics      <- c("PLA_g", "PAA_rad_s2", "PLV_m_s", "PAV_rad_s")
labels       <- c("PLA (g)", "PAA (rad/s2)", "PLV (m/s)", "PAV (rad/s)")
labels_map   <- setNames(labels, metrics)
mech_metrics <- c("DeltaPAV", "PeakAlpha", "PeakOmega")
mech_labels  <- c("DeltaPAV (rad/s)", "PeakAlpha (rad/s2)", "PeakOmega (rad/s)")

int_metrics_pav  <- c("PLA_g", "PAA_rad_s2", "PAV_rad_s")
int_labels_pav   <- c("PLA (g)", "PAA (rad/s2)", "PAV (rad/s)")
int_metrics_dpav <- c("PLA_g", "PAA_rad_s2", "DeltaPAV")
int_labels_dpav  <- c("PLA (g)", "PAA (rad/s2)", "DeltaPAV (rad/s)")

# =============================================================================
# Helper: refit interaction model (needed for emmeans on glmmTMB)
# =============================================================================

refit_interaction <- function(metric, data) {
  d   <- data |> filter(!is.na(.data[[metric]]), .data[[metric]] > 0)
  fml <- as.formula(paste0(
    metric, " ~ skill_level * situation_type",
    " + (1 | team / device_sn) + (1 | team:date)"))
  tryCatch(glmmTMB(fml, data = d, family = Gamma(link = "log")), error = function(e) NULL)
}

# =============================================================================
# 1. Main effect: EMM & pairwise (skill_level)
# =============================================================================

run_main_emm <- function(fits, labels, section) {
  cat(sprintf("\n%s\n%s\n", section, strrep("=", nchar(section))))
  emm_rows <- map2(fits, labels, function(fit, lbl) {
    if (is.null(fit$model)) { cat(sprintf("%s: model not available\n", lbl)); return(NULL) }
    emm <- emmeans(fit$model, ~ skill_level, type = "response")
    cat(sprintf("\n--- %s: Estimated marginal means ---\n", lbl))
    print(summary(emm))
    pw <- pairs(emm, adjust = "holm")
    cat(sprintf("\n--- %s: Pairwise comparison (Holm) ---\n", lbl))
    print(pw)
    print(confint(pw))
    as.data.frame(summary(emm)) |> mutate(Metric = lbl, .before = 1)
  })
  bind_rows(emm_rows)
}

emm_A1      <- run_main_emm(fits_A1,      labels,      "A1: EMM & pairwise — impact metrics")
emm_A2      <- run_main_emm(fits_A2,      labels,      "A2: EMM & pairwise — impact metrics")
emm_A1_mech <- run_main_emm(fits_A1_mech, mech_labels, "A1: EMM & pairwise — biomechanical metrics")

write_csv(emm_A1,      file.path(OUTPUT_DIR, "emm_A1.csv"))
write_csv(emm_A2,      file.path(OUTPUT_DIR, "emm_A2.csv"))
write_csv(emm_A1_mech, file.path(OUTPUT_DIR, "emm_A1_mech.csv"))

# =============================================================================
# 2. ICC by random effect group (main effect models)
# =============================================================================

get_icc <- function(fits, labels_map_or_vec, section) {
  cat(sprintf("\n%s\n%s\n", section, strrep("=", nchar(section))))
  map_dfr(names(fits), function(m) {
    fit <- fits[[m]]
    if (is.null(fit$model)) return(NULL)
    icc_res <- tryCatch(icc(fit$model, by_group = TRUE), error = function(e) NULL)
    if (is.null(icc_res)) return(NULL)
    lbl <- if (is.null(names(labels_map_or_vec))) labels_map_or_vec[m] else labels_map_or_vec[[m]]
    as.data.frame(icc_res) |> mutate(Metric = lbl, .before = 1)
  })
}

icc_A1      <- get_icc(fits_A1,      labels_map,                   "A1: ICC — impact metrics")
icc_A2      <- get_icc(fits_A2,      labels_map,                   "A2: ICC — impact metrics")
icc_A1_mech <- get_icc(fits_A1_mech, setNames(mech_labels, mech_metrics), "A1: ICC — biomechanical metrics")

print(icc_A1); print(icc_A2); print(icc_A1_mech)
write_csv(icc_A1,      file.path(OUTPUT_DIR, "icc_A1.csv"))
write_csv(icc_A2,      file.path(OUTPUT_DIR, "icc_A2.csv"))
write_csv(icc_A1_mech, file.path(OUTPUT_DIR, "icc_A1_mech.csv"))

# =============================================================================
# 3. Interaction models: interaction contrasts & EMM by situation
# =============================================================================

run_interaction_posthoc <- function(metric, label, df_int, section_prefix) {
  fit <- refit_interaction(metric, df_int)
  if (is.null(fit)) { cat(sprintf("%s: model not available\n", label)); return(invisible(NULL)) }

  cat(sprintf("\n%s — %s\n%s\n", section_prefix, label, strrep("-", 60)))

  # Step 1: interaction contrast
  cat("Interaction contrast (skill_level x situation_type, Holm):\n")
  interact <- emmeans(fit, ~ skill_level * situation_type, type = "response")
  int_con  <- contrast(interact, interaction = c("pairwise", "pairwise"), adjust = "holm")
  print(int_con)
  int_p <- as.data.frame(summary(int_con))$p.value
  if (!is.null(int_p) && !is.na(int_p[1])) {
    int_p <- int_p[1]
    if      (int_p < 0.05) cat(sprintf("  >> Significant interaction (p = %.3f)\n", int_p))
    else if (int_p < 0.10) cat(sprintf("  >> Marginal interaction (p = %.3f)\n", int_p))
    else                   cat(sprintf("  >> Non-significant interaction (p = %.3f)\n", int_p))
  }

  # Step 2: EMM by situation
  cat("\nEstimated marginal means (skill_level | situation_type):\n")
  emm <- emmeans(fit, ~ skill_level | situation_type, type = "response")
  print(summary(emm))

  cat("\nPairwise comparisons (Holm):\n")
  pw <- pairs(emm, adjust = "holm")
  print(pw)
  print(confint(pw))
}

# A1 — PAV set
message("\n--- A1: Interaction post-hoc (PAV set) ---")
walk2(int_metrics_pav, int_labels_pav,
  ~ run_interaction_posthoc(.x, .y, df_A1_int, "A1"))

# A1 — DeltaPAV set
message("\n--- A1: Interaction post-hoc (DeltaPAV set) ---")
walk2(int_metrics_dpav, int_labels_dpav,
  ~ run_interaction_posthoc(.x, .y, df_A1_int, "A1"))

# A2 — PAV set
message("\n--- A2: Interaction post-hoc (PAV set) ---")
if (!is.null(int_results_A2_pav) && nrow(int_results_A2_pav) > 0)
  walk2(int_metrics_pav, int_labels_pav,
    ~ run_interaction_posthoc(.x, .y, df_A2_int, "A2"))

# A2 — DeltaPAV set
message("\n--- A2: Interaction post-hoc (DeltaPAV set) ---")
if (!is.null(int_results_A2_dpav) && nrow(int_results_A2_dpav) > 0)
  walk2(int_metrics_dpav, int_labels_dpav,
    ~ run_interaction_posthoc(.x, .y, df_A2_int, "A2"))

message("\nDone. EMM tables saved to: ", OUTPUT_DIR)
