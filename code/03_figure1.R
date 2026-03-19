# =============================================================================
# Study 1 — Figure 1: Skill Level x Event Type Interaction Plots
#
# Run AFTER 01_glmm_analysis.R (requires output/models/*.rds).
# Run from the PROJECT ROOT:
#   Rscript code/03_figure1.R
#
# Outputs (saved to output/figures/):
#   fig1a_interaction_A1_pav.png  / .pdf   — PLA, PAA, PAV  (A1)
#   fig1b_interaction_A1_dpav.png / .pdf   — PLA, PAA, DeltaPAV (A1)
#   fig1c_interaction_A2_pav.png  / .pdf   — PLA, PAA, PAV  (A2)
#   fig1d_interaction_A2_dpav.png / .pdf   — PLA, PAA, DeltaPAV (A2)
#
# Plot design:
#   x-axis  : Event type (BRE / non-BRE)
#   y-axis  : Estimated marginal mean (GLMM, Gamma log-link, response scale)
#   colour  : Skill level (Novice / Experienced)
#   error bar: 95% CI
#   faint dots: 50% random sample of raw data (jitter-dodge)
#   p-label : interaction p-value per panel
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(emmeans)
  library(glmmTMB)
  library(patchwork)
})

MODEL_DIR  <- "output/models"
FIG_DIR    <- "output/figures"
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Figure dimensions
FIG_W <- 9; FIG_H <- 4.5; FIG_DPI <- 300

# =============================================================================
# Load pre-fitted interaction results and data
# =============================================================================

int_results_A1_pav  <- readRDS(file.path(MODEL_DIR, "int_results_A1_pav.rds"))
int_results_A1_dpav <- readRDS(file.path(MODEL_DIR, "int_results_A1_dpav.rds"))
int_results_A2_pav  <- readRDS(file.path(MODEL_DIR, "int_results_A2_pav.rds"))
int_results_A2_dpav <- readRDS(file.path(MODEL_DIR, "int_results_A2_dpav.rds"))
df_A1_int <- readRDS(file.path(MODEL_DIR, "df_A1_int.rds"))
df_A2_int <- readRDS(file.path(MODEL_DIR, "df_A2_int.rds"))

# =============================================================================
# Colour / shape palette
# =============================================================================

PAL  <- c("Novice" = "#2980B9", "Experienced" = "#C0392B")
SHP  <- c("Novice" = 21L,       "Experienced" = 24L)
DODGE <- position_dodge(width = 0.2)

recode_skill <- function(x)
  factor(x, levels = c("beginner", "experienced"), labels = c("Novice", "Experienced"))

recode_sit <- function(x) factor(x, levels = c("BRE", "non-BRE"))

# =============================================================================
# Helper: refit interaction model (glmmTMB) for emmeans
# =============================================================================

refit_interaction <- function(metric, data) {
  d   <- data |> filter(!is.na(.data[[metric]]), .data[[metric]] > 0)
  fml <- as.formula(paste0(
    metric, " ~ skill_level * situation_type",
    " + (1 | team / device_sn) + (1 | team:date)"))
  tryCatch(glmmTMB(fml, data = d, family = Gamma(link = "log")), error = function(e) NULL)
}

# =============================================================================
# Core plot function
# =============================================================================

make_interaction_plot <- function(int_results, df_int,
                                   plot_metrics, plot_labels,
                                   title_suffix = "") {
  # --- Compute EMM -----------------------------------------------------------
  emm_list <- map2(plot_metrics, plot_labels, function(metric, label) {
    fit <- refit_interaction(metric, df_int)
    if (is.null(fit)) return(NULL)
    emmeans(fit, ~ skill_level * situation_type, type = "response")
  })

  plot_df <- map2_dfr(plot_labels, emm_list, function(lbl, emm) {
    if (is.null(emm)) return(NULL)
    df <- as.data.frame(summary(emm, infer = c(TRUE, FALSE), level = 0.95))
    # Normalise CI column names across emmeans versions
    if (!"lower.CL" %in% names(df)) {
      lo <- grep("lower|lcl", names(df), value = TRUE, ignore.case = TRUE)[1]
      hi <- grep("upper|ucl", names(df), value = TRUE, ignore.case = TRUE)[1]
      df <- rename(df, lower.CL = all_of(lo), upper.CL = all_of(hi))
    }
    df |> mutate(
      Metric         = factor(lbl, levels = plot_labels),
      skill_level    = recode_skill(skill_level),
      situation_type = recode_sit(situation_type)
    )
  })

  if (nrow(plot_df) == 0) {
    message("No EMM data available for: ", title_suffix)
    return(invisible(NULL))
  }

  # --- Annotation: interaction p-value per panel ----------------------------
  annot_df <- int_results |>
    filter(grepl("skill x situation", Term)) |>
    transmute(
      Metric  = factor(Metric, levels = plot_labels),
      p_label = paste0("p = ", p_display)
    )

  # --- Raw data (50% jittered subsample) ------------------------------------
  raw_df <- df_int |>
    filter(!is.na(skill_level), !is.na(situation_type)) |>
    slice_sample(prop = 0.5) |>
    select(situation_type, skill_level, all_of(plot_metrics)) |>
    pivot_longer(all_of(plot_metrics), names_to = "mkey", values_to = "value") |>
    filter(!is.na(value), value > 0) |>
    mutate(
      Metric         = factor(mkey, levels = plot_metrics, labels = plot_labels),
      skill_level    = recode_skill(skill_level),
      situation_type = recode_sit(situation_type)
    )

  # --- Build plot ------------------------------------------------------------
  ggplot(plot_df,
      aes(x = situation_type, y = response,
          colour = skill_level, fill = skill_level,
          group  = skill_level, shape = skill_level)) +
    # Raw data jittered in background
    geom_point(
      data = raw_df,
      aes(x = situation_type, y = value,
          colour = skill_level, group = skill_level),
      position    = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.2),
      alpha       = 0.07, size = 0.6, shape = 16,
      inherit.aes = FALSE, show.legend = FALSE
    ) +
    # EMM: connecting lines, points, error bars
    geom_line(linewidth = 1.1, position = DODGE) +
    geom_point(size = 3.2, stroke = 1.6, position = DODGE) +
    geom_errorbar(
      aes(ymin = lower.CL, ymax = upper.CL),
      width = 0.12, linewidth = 0.9, position = DODGE
    ) +
    # Interaction p-value annotation (top-right of each panel)
    geom_text(
      data = annot_df,
      aes(x = Inf, y = Inf, label = p_label),
      inherit.aes = FALSE,
      hjust = 1.05, vjust = 1.6, size = 3.3,
      colour = "gray25", fontface = "italic"
    ) +
    facet_wrap(~ Metric, scales = "free_y", ncol = length(plot_metrics)) +
    scale_colour_manual(values = PAL, name = NULL) +
    scale_fill_manual(values   = PAL, name = NULL) +
    scale_shape_manual(values  = SHP, name = NULL) +
    labs(
      title   = paste0("Skill level \u00d7 event type interaction", title_suffix),
      x       = "Event type",
      y       = "Estimated marginal mean (response scale)",
      caption = paste(
        "Points & lines: GLMM estimated marginal means (Gamma, log-link); bars: 95% CI.",
        "Faint dots: 50% random sample of raw observations.",
        "p-values: Holm-adjusted interaction contrast."
      )
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position    = "top",
      legend.direction   = "horizontal",
      strip.background   = element_rect(fill = "#2C3E50", colour = NA),
      strip.text         = element_text(colour = "white", face = "bold", size = 11),
      panel.grid.major.y = element_line(colour = "gray92", linewidth = 0.4),
      panel.grid.minor   = element_blank(),
      plot.title         = element_text(face = "bold", size = 13),
      plot.caption       = element_text(colour = "gray50", size = 8)
    )
}

# =============================================================================
# Generate and save figures
# =============================================================================

save_fig <- function(p, stem) {
  if (is.null(p)) return(invisible(NULL))
  ggsave(file.path(FIG_DIR, paste0(stem, ".png")), p, width = FIG_W, height = FIG_H, dpi = FIG_DPI)
  ggsave(file.path(FIG_DIR, paste0(stem, ".pdf")), p, width = FIG_W, height = FIG_H)
  message("Saved: ", stem)
}

int_metrics_pav  <- c("PLA_g", "PAA_rad_s2", "PAV_rad_s")
int_labels_pav   <- c("PLA (g)", "PAA (rad/s\u00b2)", "PAV (rad/s)")
int_metrics_dpav <- c("PLA_g", "PAA_rad_s2", "DeltaPAV")
int_labels_dpav  <- c("PLA (g)", "PAA (rad/s\u00b2)", "\u0394PAV (rad/s)")

message("--- Generating Figure 1 panels ---")

# Fig 1a: A1 — PLA, PAA, PAV
p1a <- make_interaction_plot(
  int_results_A1_pav, df_A1_int,
  int_metrics_pav, int_labels_pav,
  title_suffix = " \u2014 Dataset A1 (PLA, PAA, PAV)"
)
save_fig(p1a, "fig1a_interaction_A1_pav")

# Fig 1b: A1 — PLA, PAA, DeltaPAV
p1b <- make_interaction_plot(
  int_results_A1_dpav, df_A1_int,
  int_metrics_dpav, int_labels_dpav,
  title_suffix = " \u2014 Dataset A1 (PLA, PAA, \u0394PAV)"
)
save_fig(p1b, "fig1b_interaction_A1_dpav")

# Fig 1c: A2 — PLA, PAA, PAV
p1c <- make_interaction_plot(
  int_results_A2_pav, df_A2_int,
  int_metrics_pav, int_labels_pav,
  title_suffix = " \u2014 Dataset A2 (PLA, PAA, PAV)"
)
save_fig(p1c, "fig1c_interaction_A2_pav")

# Fig 1d: A2 — PLA, PAA, DeltaPAV
p1d <- make_interaction_plot(
  int_results_A2_dpav, df_A2_int,
  int_metrics_dpav, int_labels_dpav,
  title_suffix = " \u2014 Dataset A2 (PLA, PAA, \u0394PAV)"
)
save_fig(p1d, "fig1d_interaction_A2_dpav")

message("\nDone. Figures saved to: ", FIG_DIR)
