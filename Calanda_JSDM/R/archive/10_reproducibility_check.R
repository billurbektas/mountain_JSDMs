# ==============================================================================
# Script: 10_reproducibility_check.R
# Purpose: Assess reproducibility of sjSDM variance partitioning by comparing
#          multiple fits with identical config (LIN_XY_XY, alpha=1, lambda=0.01)
#
#          Sources compared:
#            - Repro Run1 & Run2: same Colab session, same seed
#            - Exp1: LIN_XY_XY_a1_l0.01 (earlier session)
#            - Exp2: 3 decoupled runs at lambda_vary=0.01, alpha=1 (another session)
#
# Inputs:
#   - output/results/repro_check_run1.rds
#   - output/results/repro_check_run2.rds
#   - output/results/runs/LIN_XY_XY_a1_l0.01.rds
#   - output/results/decoupled/decoupled_LIN_XY_XY_lambda_{env,sp,bio}_0.01_a1.rds
#
# Outputs:
#   - plot/repro_check.pdf
#
# Requires: R/00_setup/functions_calanda.R
# ==============================================================================

library(tidyverse)
library(conflicted)
library(here)
library(patchwork)

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::select)

source(here("Calanda_JSDM", "R", "00_setup", "functions_calanda.R"))

# Colors
color_env     = "#81caf3"
color_spa     = "#d00000"
color_codist  = "#00bd89"
color_overall = "black"

comp_colors = c(env = color_env, spa = color_spa, codist = color_codist, r2 = color_overall)
comp_labels = c(env = "Environment", spa = "Spatial", codist = "Biotic", r2 = "R\u00B2")

# ==============================================================================
# LOAD ALL RUNS
# ==============================================================================
cat("\n=== Loading all runs with config: LIN_XY_XY, alpha=1, lambda=0.01 ===\n")

runs = list()

# Repro runs (same Colab session)
runs[["Repro 1"]] = readRDS(here("Calanda_JSDM", "output", "results", "repro_check_run1.rds"))
runs[["Repro 2"]] = readRDS(here("Calanda_JSDM", "output", "results", "repro_check_run2.rds"))

# Exp1 run
runs[["Exp1"]] = readRDS(here("Calanda_JSDM", "output", "results", "runs",
                               "LIN_XY_XY_a1_l0.01.rds"))

# Exp2 decoupled runs (all 3 at anchor = identical config)
runs[["Exp2 env"]] = readRDS(here("Calanda_JSDM", "output", "results", "decoupled",
                                   "decoupled_LIN_XY_XY_lambda_env_0.01_a1.rds"))
runs[["Exp2 spa"]] = readRDS(here("Calanda_JSDM", "output", "results", "decoupled",
                                   "decoupled_LIN_XY_XY_lambda_sp_0.01_a1.rds"))
runs[["Exp2 bio"]] = readRDS(here("Calanda_JSDM", "output", "results", "decoupled",
                                   "decoupled_LIN_XY_XY_lambda_bio_0.01_a1.rds"))

cat(sprintf("Loaded %d runs\n", length(runs)))

# ==============================================================================
# PAIRWISE COMPARISONS
# ==============================================================================
cat("\n=== Pairwise comparisons ===\n")

run_names = names(runs)
pairs = combn(run_names, 2, simplify = FALSE)

comparison_rows = list()

for (pair in pairs) {
  a_name = pair[1]
  b_name = pair[2]
  a = runs[[a_name]]
  b = runs[[b_name]]

  r2_diff = b$partition$R2 - a$partition$R2

  cat(sprintf("\n  %s vs %s\n", a_name, b_name))
  cat(sprintf("    Overall R2: %.6f vs %.6f  (diff = %.6f)\n",
              a$partition$R2, b$partition$R2, r2_diff))

  for (col in c("env", "spa", "codist", "r2")) {
    sp_r = cor(a$partition$species[, col], b$partition$species[, col], use = "complete.obs")
    si_r = cor(a$partition$sites[, col],   b$partition$sites[, col],   use = "complete.obs")

    comparison_rows = c(comparison_rows, list(tibble(
      run_a     = a_name,
      run_b     = b_name,
      pair      = paste0(a_name, " vs ", b_name),
      component = col,
      species_r = sp_r,
      sites_r   = si_r,
      r2_diff   = r2_diff
    )))

    cat(sprintf("    %-7s species r = %.4f   sites r = %.4f\n", col, sp_r, si_r))
  }
}

df_comp = bind_rows(comparison_rows) %>%
  mutate(component = factor(component, levels = c("r2", "env", "spa", "codist")))

# Classify pairs by type
df_comp = df_comp %>%
  mutate(
    pair_type = case_when(
      grepl("Repro.*Repro", pair)    ~ "Same session",
      grepl("Exp2.*Exp2", pair)      ~ "Same session (Exp2)",
      TRUE                           ~ "Cross-session"
    ),
    pair_type = factor(pair_type, levels = c("Same session", "Same session (Exp2)", "Cross-session"))
  )

# ==============================================================================
# FIGURE
# ==============================================================================
cat("\n=== Creating reproducibility figure ===\n")

# --- Panel A: Overall R2 differences per pair ---
df_r2 = df_comp %>%
  distinct(pair, pair_type, r2_diff)

p_r2 = ggplot(df_r2, aes(x = pair, y = r2_diff, fill = pair_type)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  geom_col(width = 0.6, alpha = 0.7) +
  geom_text(aes(label = sprintf("%.4f", r2_diff)),
            vjust = ifelse(df_r2$r2_diff >= 0, -0.5, 1.5), size = 3) +
  scale_fill_manual(values = c("Same session" = "#2ecc71",
                                "Same session (Exp2)" = "#3498db",
                                "Cross-session" = "#e74c3c")) +
  labs(
    title = "A) Overall R\u00B2 difference per pair",
    x = NULL, y = expression(Delta * R^2),
    fill = "Comparison type"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

# --- Panel B: Species-level correlations ---
p_sp = ggplot(df_comp, aes(x = pair, y = species_r, fill = component)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  geom_hline(yintercept = 1, linewidth = 0.3, linetype = "dashed", color = "grey60") +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.7) +
  scale_fill_manual(values = comp_colors, labels = comp_labels) +
  labs(
    title = "B) Species-level correlation per pair",
    x = NULL, y = "Pearson r",
    fill = "Component"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

# --- Panel C: Site-level correlations ---
p_si = ggplot(df_comp, aes(x = pair, y = sites_r, fill = component)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  geom_hline(yintercept = 1, linewidth = 0.3, linetype = "dashed", color = "grey60") +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.7) +
  scale_fill_manual(values = comp_colors, labels = comp_labels) +
  scale_y_continuous(limits = c(0, 1.08)) +
  labs(
    title = "C) Site-level correlation per pair",
    x = NULL, y = "Pearson r",
    fill = "Component"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

# --- Combine ---
p_repro = p_r2 / p_sp / p_si +
  plot_annotation(
    title = "Reproducibility of variance partitioning across fits (LIN_XY_XY, alpha = 1, lambda = 0.01)",
    subtitle = "Same session = same seed in one Colab run. Cross-session = different Colab runs.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10)
    )
  )

pdf(here("Calanda_JSDM", "plot", "repro_check.pdf"),
    width = 14, height = 14)
print(p_repro)
dev.off()
cat("  Saved repro_check.pdf\n")

cat("\n=== Reproducibility check complete ===\n")
