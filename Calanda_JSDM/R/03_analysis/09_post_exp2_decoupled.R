# ==============================================================================
# Script: 09_post_exp2_decoupled.R
# Purpose: Post-analysis of sjSDM Experiment 2 - decoupled lambda sensitivity
#          LIN_XY_XY only, 3 compartments x 3 lambdas x 3 alphas = 27 runs
#          One compartment's lambda varies (0.001, 0.01, 0.1), others anchored
#          at 0.01. Alpha (0, 0.5, 1) shared across compartments.
#
# Inputs:
#   - output/results/decoupled/*.rds (27 cached runs)
#
# Outputs:
#   - plot/exp2_convergence.pdf
#   - plot/exp2_species_violins.pdf
#   - plot/exp2_sites_violins.pdf
#   - plot/exp2_anova_R2.pdf
#   - output/results/exp2_model_metrics.csv
#
# Requires: R/00_setup/functions_calanda.R
# ==============================================================================

library(tidyverse)
library(conflicted)
library(here)
library(patchwork)
library(ggrepel)

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::select)

source(here("Calanda_JSDM", "R", "00_setup", "functions_calanda.R"))

# Colors (consistent with project conventions)
color_env     = "#81caf3"
color_spa     = "#d00000"
color_codist  = "#00bd89"
color_overall = "black"

# Ensure output directories exist
dir.create(here("Calanda_JSDM", "plot"),
           showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# STEP 1: LOAD DATA
# ==============================================================================
cat("\n=== Step 1: Loading data ===\n")

run_files = list.files(
  here("Calanda_JSDM", "output", "results", "decoupled"),
  pattern = "\\.rds$", full.names = TRUE
)
cat(sprintf("Found %d run files\n", length(run_files)))

# Extract data from each run
history_list  = list()
metrics_list  = list()
species_list  = list()
sites_list    = list()

for (f in run_files) {
  run = readRDS(f)
  cfg = run$config
  rid = cfg$run_id
  vary = cfg$vary
  lv   = cfg$lambda_vary
  a    = cfg$alpha_common

  # Training history
  history_list[[rid]] = tibble(
    run_id       = rid,
    vary         = vary,
    lambda_vary  = lv,
    alpha        = a,
    iteration    = seq_along(run$model$history),
    loss         = as.numeric(run$model$history)
  )

  # Model-level McFadden R² from anova results
  anova_res = run$partition$anova$results
  r2_mcf    = setNames(anova_res$`R2 McFadden`, anova_res$models)

  metrics_list[[rid]] = tibble(
    run_id       = rid,
    vary         = vary,
    lambda_vary  = lv,
    alpha        = a,
    metric_type  = c("overall_R2", "env_R2", "spa_R2", "bio_R2",
                     "env_bio_R2", "env_spa_R2", "bio_spa_R2", "shared_R2"),
    value        = c(r2_mcf["Full"], r2_mcf["F_A"], r2_mcf["F_S"], r2_mcf["F_B"],
                     r2_mcf["F_AB"], r2_mcf["F_AS"], r2_mcf["F_BS"], r2_mcf["F_ABS"])
  )

  # Species-level components (352 x 4)
  species_list[[rid]] = as_tibble(run$partition$species) %>%
    mutate(run_id = rid, vary = vary, lambda_vary = lv, alpha = a,
           species_idx = row_number())

  # Site-level components (576 x 4)
  sites_list[[rid]] = as_tibble(run$partition$sites) %>%
    mutate(run_id = rid, vary = vary, lambda_vary = lv, alpha = a,
           site_idx = row_number())
}

df_history = bind_rows(history_list)
df_metrics = bind_rows(metrics_list)
df_species = bind_rows(species_list)
df_sites   = bind_rows(sites_list)

cat(sprintf("  History: %d rows\n", nrow(df_history)))
cat(sprintf("  Metrics: %d rows\n", nrow(df_metrics)))
cat(sprintf("  Species: %d rows (%d species x %d runs)\n",
            nrow(df_species), n_distinct(df_species$species_idx),
            n_distinct(df_species$run_id)))
cat(sprintf("  Sites:   %d rows (%d sites x %d runs)\n",
            nrow(df_sites), n_distinct(df_sites$site_idx),
            n_distinct(df_sites$run_id)))

# Nice labels for the "vary" facet
vary_labels = c(
  lambda_env = "Vary lambda[env]",
  lambda_sp  = "Vary lambda[spa]",
  lambda_bio = "Vary lambda[bio]"
)

# ==============================================================================
# SANITY CHECK: Exp2 lambda_vary=0.01 vs Exp1 LIN_XY_XY l=0.01
# ==============================================================================
cat("\n=== Sanity check: Exp2 (lambda_vary=0.01) vs Exp1 (LIN_XY_XY, l=0.01) ===\n")

sanity_rows = list()

for (alpha_val in c(0, 0.5, 1)) {

  # --- Exp1 reference ---
  exp1_file = here("Calanda_JSDM", "output", "results", "runs",
                   paste0("LIN_XY_XY_a", alpha_val, "_l0.01.rds"))
  if (!file.exists(exp1_file)) {
    cat("  Exp1 file not found for alpha =", alpha_val, "- skipping\n")
    next
  }
  exp1 = readRDS(exp1_file)
  exp1_r2      = exp1$partition$R2
  exp1_species = exp1$partition$species
  exp1_sites   = exp1$partition$sites

  cat(sprintf("\n  alpha = %s\n", alpha_val))
  cat("  ", strrep("-", 50), "\n")

  # --- Exp2: pick first vary run (all three are identical at anchor) ---
  exp2_rid = df_species %>%
    filter(lambda_vary == 0.01, alpha == alpha_val) %>%
    pull(run_id) %>% unique() %>% .[1]

  exp2 = readRDS(here("Calanda_JSDM", "output", "results", "decoupled",
                      paste0(exp2_rid, ".rds")))
  exp2_r2      = exp2$partition$R2
  exp2_species = exp2$partition$species
  exp2_sites   = exp2$partition$sites

  r2_diff = exp2_r2 - exp1_r2
  cat(sprintf("  Overall R2: Exp1=%.6f  Exp2=%.6f  diff=%.6f\n",
              exp1_r2, exp2_r2, r2_diff))

  for (col in c("env", "spa", "codist", "r2")) {
    sp_cor = cor(exp1_species[, col], exp2_species[, col], use = "complete.obs")
    si_cor = cor(exp1_sites[, col],   exp2_sites[, col],   use = "complete.obs")

    sanity_rows = c(sanity_rows, list(tibble(
      alpha     = alpha_val,
      component = col,
      level     = c("Species", "Sites"),
      correlation = c(sp_cor, si_cor)
    )))
  }

  # Overall R2 as a separate entry
  sanity_rows = c(sanity_rows, list(tibble(
    alpha       = alpha_val,
    component   = "overall_R2",
    level       = "Model",
    correlation = NA_real_,
    r2_diff     = r2_diff
  )))

  sp_cors = sapply(c("env", "spa", "codist", "r2"), function(col)
    cor(exp1_species[, col], exp2_species[, col], use = "complete.obs"))
  si_cors = sapply(c("env", "spa", "codist", "r2"), function(col)
    cor(exp1_sites[, col], exp2_sites[, col], use = "complete.obs"))

  cat(sprintf("  Species cor:  env=%.4f  spa=%.4f  codist=%.4f  r2=%.4f\n",
              sp_cors["env"], sp_cors["spa"], sp_cors["codist"], sp_cors["r2"]))
  cat(sprintf("  Sites cor:    env=%.4f  spa=%.4f  codist=%.4f  r2=%.4f\n",
              si_cors["env"], si_cors["spa"], si_cors["codist"], si_cors["r2"]))
}

df_sanity = bind_rows(sanity_rows)

cat("\n  (Correlations should be ~1.0 and R2 diffs ~0 if runs are reproducible)\n")

# --- Sanity check figure ---
cat("\n=== Creating sanity check figure ===\n")

comp_labels = c(env = "Environment", spa = "Spatial", codist = "Biotic", r2 = "R\u00B2")
comp_colors_sanity = c(env = color_env, spa = color_spa, codist = color_codist, r2 = color_overall)

# Panel A: Overall R2 difference
df_r2_diff = df_sanity %>%
  filter(level == "Model") %>%
  mutate(alpha_lab = paste0("alpha = ", alpha))

p_r2 = ggplot(df_r2_diff, aes(x = factor(alpha), y = r2_diff)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  geom_col(fill = "grey50", width = 0.5) +
  geom_text(aes(label = sprintf("%.4f", r2_diff)),
            vjust = ifelse(df_r2_diff$r2_diff >= 0, -0.5, 1.5), size = 3.5) +
  labs(
    title = "A) Overall R\u00B2 difference (Exp2 - Exp1)",
    x = expression(alpha), y = expression(Delta * R^2)
  ) +
  theme_bw(base_size = 11)

# Panel B: Species-level correlations
df_sp_cor = df_sanity %>%
  filter(level == "Species") %>%
  mutate(
    component_lab = recode(component, !!!comp_labels),
    component_lab = factor(component_lab, levels = comp_labels)
  )

p_sp = ggplot(df_sp_cor, aes(x = factor(alpha), y = correlation, fill = component)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  geom_hline(yintercept = 1, linewidth = 0.3, linetype = "dashed", color = "grey60") +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.7) +
  geom_text(aes(label = sprintf("%.2f", correlation)),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 2.8) +
  scale_fill_manual(values = comp_colors_sanity, labels = comp_labels) +
  labs(
    title = "B) Species-level correlation (Exp1 vs Exp2)",
    x = expression(alpha), y = "Pearson r",
    fill = "Component"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# Panel C: Site-level correlations
df_si_cor = df_sanity %>%
  filter(level == "Sites") %>%
  mutate(
    component_lab = recode(component, !!!comp_labels),
    component_lab = factor(component_lab, levels = comp_labels)
  )

p_si = ggplot(df_si_cor, aes(x = factor(alpha), y = correlation, fill = component)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  geom_hline(yintercept = 1, linewidth = 0.3, linetype = "dashed", color = "grey60") +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.7) +
  geom_text(aes(label = sprintf("%.2f", correlation)),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 2.8) +
  scale_fill_manual(values = comp_colors_sanity, labels = comp_labels) +
  scale_y_continuous(limits = c(0, 1.08)) +
  labs(
    title = "C) Site-level correlation (Exp1 vs Exp2)",
    x = expression(alpha), y = "Pearson r",
    fill = "Component"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# Combine with patchwork
p_sanity = p_r2 / (p_sp | p_si) +
  plot_annotation(
    title = "Reproducibility check: Exp1 vs Exp2 at identical config (LIN_XY_XY, lambda = 0.01)",
    subtitle = "Same hyperparameters, independent model fits",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 11)
    )
  )

pdf(here("Calanda_JSDM", "plot", "exp2_sanity_check.pdf"),
    width = 12, height = 10)
print(p_sanity)
dev.off()
cat("  Saved exp2_sanity_check.pdf\n")

# ==============================================================================
# STEP 2: CONVERGENCE CHECK
# ==============================================================================
cat("\n=== Step 2: Convergence plot ===\n")

p_conv = ggplot(df_history, aes(x = iteration, y = loss, color = run_id)) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "solid", color = "grey40") +
  geom_line(alpha = 0.6, linewidth = 0.4) +
  labs(title = "Experiment 2: Training convergence - decoupled lambda",
       x = "Iteration", y = "Loss") +
  theme_bw() +
  theme(legend.position = "none")

pdf(here("Calanda_JSDM", "plot", "exp2_convergence.pdf"),
    width = 10, height = 6)
print(p_conv)
dev.off()
cat("  Saved exp2_convergence.pdf\n")

# ==============================================================================
# STEP 3: MODEL-LEVEL METRICS TABLE
# ==============================================================================
cat("\n=== Step 3: Model-level metrics ===\n")

write_csv(
  df_metrics,
  here("Calanda_JSDM", "output", "results", "exp2_model_metrics.csv")
)
cat("  Saved exp2_model_metrics.csv\n")

# ==============================================================================
# STEP 4: SPECIES VIOLIN PLOT
# ==============================================================================
cat("\n=== Step 4: Species violin plot ===\n")

# Violin color mappings (unit-level: env, spa, codist, r2)
violin_colors = c(
  "Overall"     = color_overall,
  "Environment" = color_env,
  "Spatial"     = color_spa,
  "Biotic"      = color_codist
)

# Helper: pivot unit-level data to long format
pivot_unit_long = function(df) {
  df %>%
    pivot_longer(
      cols      = c(env, spa, codist, r2),
      names_to  = "component_raw",
      values_to = "value"
    ) %>%
    mutate(
      component = recode(component_raw,
        r2 = "Overall", env = "Environment",
        spa = "Spatial", codist = "Biotic"
      ),
      component = factor(component, levels = names(violin_colors)),
      lambda_f  = factor(lambda_vary),
      alpha_lab = paste0("alpha = ", alpha)
    )
}

# Helper: build violin-only plot for exp2
make_violin_plot = function(unit_long, title_label) {
  ggplot(unit_long, aes(x = lambda_f, y = value, fill = component)) +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = "solid", color = "grey40") +
    geom_violin(
      alpha = 0.3, scale = "width",
      position = position_dodge(width = 0.7),
      draw_quantiles = 0.5, linewidth = 0.2
    ) +
    scale_fill_manual(values = violin_colors) +
    facet_grid(alpha_lab ~ vary, labeller = labeller(
      vary = vary_labels
    )) +
    labs(
      title = sprintf("Experiment 2: %s-level variance partitioning (decoupled lambda)", title_label),
      subtitle = "Violins with median line. Others anchored at lambda = 0.01",
      x = expression(lambda[varied]), y = "Value",
      fill = "Component"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold")
    )
}

df_species_long = pivot_unit_long(df_species)

p_species = make_violin_plot(df_species_long, "Species")

pdf(here("Calanda_JSDM", "plot", "exp2_species_violins.pdf"),
    width = 14, height = 10)
print(p_species)
dev.off()
cat("  Saved exp2_species_violins.pdf\n")

# ==============================================================================
# STEP 5: SITE VIOLIN PLOT
# ==============================================================================
cat("\n=== Step 5: Site violin plot ===\n")

df_sites_long = pivot_unit_long(df_sites)

p_sites = make_violin_plot(df_sites_long, "Site")

pdf(here("Calanda_JSDM", "plot", "exp2_sites_violins.pdf"),
    width = 14, height = 10)
print(p_sites)
dev.off()
cat("  Saved exp2_sites_violins.pdf\n")

# ==============================================================================
# STEP 6: MODEL-LEVEL ANOVA R2 (connected dots)
# ==============================================================================
cat("\n=== Step 6: Model-level anova R2 plot ===\n")

# Color mappings for all 8 anova components
color_env_bio = "#41a4ae"
color_env_spa = "#b0656a"
color_bio_spa = "#6b5f45"
color_shared  = "#9b59b6"

anova_colors = c(
  "Overall"     = color_overall,
  "Environment" = color_env,
  "Spatial"     = color_spa,
  "Biotic"      = color_codist,
  "Env x Bio"   = color_env_bio,
  "Env x Spa"   = color_env_spa,
  "Bio x Spa"   = color_bio_spa,
  "Shared"      = color_shared
)

df_model_lines = df_metrics %>%
  mutate(
    component = recode(metric_type,
      overall_R2 = "Overall", env_R2 = "Environment",
      spa_R2 = "Spatial", bio_R2 = "Biotic",
      env_bio_R2 = "Env x Bio", env_spa_R2 = "Env x Spa",
      bio_spa_R2 = "Bio x Spa", shared_R2 = "Shared"
    ),
    component = factor(component, levels = names(anova_colors)),
    lambda_f  = factor(lambda_vary),
    alpha_lab = paste0("alpha = ", alpha)
  )

# Labels at the end (rightmost lambda) of each line
df_model_labels = df_model_lines %>%
  filter(lambda_vary == max(lambda_vary))

p_anova = ggplot(df_model_lines,
                 aes(x = lambda_f, y = value, color = component, group = component)) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "solid", color = "grey40") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_text_repel(
    data = df_model_labels,
    aes(label = component),
    nudge_x = 0.3, direction = "y", hjust = 0,
    size = 3, segment.size = 0.3, show.legend = FALSE
  ) +
  scale_color_manual(values = anova_colors) +
  facet_grid(alpha_lab ~ vary, labeller = labeller(
    vary = vary_labels
  )) +
  labs(
    title = "Experiment 2: Model-level McFadden R\u00B2 (decoupled lambda)",
    subtitle = "Connected dots per anova Venn fraction. Others anchored at lambda = 0.01",
    x = expression(lambda[varied]), y = "McFadden R\u00B2",
    color = "Component"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  ) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5.5, 40, 5.5, 5.5))

pdf(here("Calanda_JSDM", "plot", "exp2_anova_R2.pdf"),
    width = 14, height = 10)
print(p_anova)
dev.off()
cat("  Saved exp2_anova_R2.pdf\n")

cat("\n=== Post Experiment 2 analysis complete ===\n")
