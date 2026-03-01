# ==============================================================================
# Script: 08_post_exp1_grid.R
# Purpose: Post-analysis of sjSDM Experiment 1 grid search results
#          (2 spatial forms x 3 alphas x 5 lambdas = 30 runs)
#
# Inputs:
#   - output/results/summary_exp1_grid.csv
#   - output/results/runs/*.rds (30 cached model runs)
#
# Outputs:
#   - plot/exp1_convergence.pdf
#   - plot/exp1_species_violins.pdf
#   - plot/exp1_sites_violins.pdf
#   - plot/exp1_anova_R2.pdf
#   - output/results/exp1_model_metrics.csv
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

summary_df = read_csv(
  here("Calanda_JSDM", "output", "results", "summary_exp1_grid.csv"),
  show_col_types = FALSE
)

run_files = list.files(
  here("Calanda_JSDM", "output", "results", "runs"),
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
  sf  = cfg$spatial_form
  a   = cfg$alpha_common
  l   = cfg$lambda_common

  # Training history
  history_list[[rid]] = tibble(
    run_id       = rid,
    spatial_form = sf,
    alpha        = a,
    lambda       = l,
    iteration    = seq_along(run$model$history),
    loss         = as.numeric(run$model$history)
  )

  # Model-level McFadden R² from anova results
  anova_res = run$partition$anova$results
  r2_mcf    = setNames(anova_res$`R2 McFadden`, anova_res$models)

  metrics_list[[rid]] = tibble(
    run_id       = rid,
    spatial_form = sf,
    alpha        = a,
    lambda       = l,
    metric_type  = c("overall_R2", "env_R2", "spa_R2", "bio_R2",
                     "env_bio_R2", "env_spa_R2", "bio_spa_R2", "shared_R2"),
    value        = c(r2_mcf["Full"], r2_mcf["F_A"], r2_mcf["F_S"], r2_mcf["F_B"],
                     r2_mcf["F_AB"], r2_mcf["F_AS"], r2_mcf["F_BS"], r2_mcf["F_ABS"])
  )

  # Species-level components (352 x 4 matrix: env, spa, codist, r2)
  species_list[[rid]] = as_tibble(run$partition$species) %>%
    mutate(run_id = rid, spatial_form = sf, alpha = a, lambda = l,
           species_idx = row_number())

  # Site-level components (576 x 4 matrix: env, spa, codist, r2)
  sites_list[[rid]] = as_tibble(run$partition$sites) %>%
    mutate(run_id = rid, spatial_form = sf, alpha = a, lambda = l,
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

# ==============================================================================
# STEP 2: CONVERGENCE CHECK
# ==============================================================================
cat("\n=== Step 2: Convergence plot ===\n")

p_conv = ggplot(df_history, aes(x = iteration, y = loss, color = run_id)) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "solid", color = "grey40") +
  geom_line(alpha = 0.6, linewidth = 0.4) +
  labs(title = "Training convergence - all 30 runs",
       x = "Iteration", y = "Loss") +
  theme_bw() +
  theme(legend.position = "none")

pdf(here("Calanda_JSDM", "plot", "exp1_convergence.pdf"),
    width = 10, height = 6)
print(p_conv)
dev.off()
cat("  Saved exp1_convergence.pdf\n")

# ==============================================================================
# STEP 3: MODEL-LEVEL METRICS TABLE
# ==============================================================================
cat("\n=== Step 3: Model-level metrics ===\n")

write_csv(
  df_metrics,
  here("Calanda_JSDM", "output", "results", "exp1_model_metrics.csv")
)
cat("  Saved exp1_model_metrics.csv\n")

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
      lambda_f  = factor(lambda),
      alpha_lab = paste0("alpha = ", alpha)
    )
}

# Helper: build violin-only plot
make_violin_plot = function(unit_long, title_label) {
  ggplot(unit_long, aes(x = lambda_f, y = value, fill = component)) +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = "solid", color = "grey40") +
    geom_violin(
      alpha = 0.3, scale = "width",
      position = position_dodge(width = 0.7),
      draw_quantiles = 0.5, linewidth = 0.2
    ) +
    scale_fill_manual(values = violin_colors) +
    facet_grid(spatial_form ~ alpha_lab) +
    labs(
      title = sprintf("Experiment 1: %s-level variance partitioning", title_label),
      subtitle = "Violins with median line",
      x = expression(lambda), y = "Value",
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

pdf(here("Calanda_JSDM", "plot", "exp1_species_violins.pdf"),
    width = 14, height = 8)
print(p_species)
dev.off()
cat("  Saved exp1_species_violins.pdf\n")

# ==============================================================================
# STEP 5: SITE VIOLIN PLOT
# ==============================================================================
cat("\n=== Step 5: Site violin plot ===\n")

df_sites_long = pivot_unit_long(df_sites)

p_sites = make_violin_plot(df_sites_long, "Site")

pdf(here("Calanda_JSDM", "plot", "exp1_sites_violins.pdf"),
    width = 14, height = 8)
print(p_sites)
dev.off()
cat("  Saved exp1_sites_violins.pdf\n")

# ==============================================================================
# STEP 6: MODEL-LEVEL ANOVA R2 (connected dots)
# ==============================================================================
cat("\n=== Step 6: Model-level anova R2 plot ===\n")

# Color mappings for all 8 anova components
# Pairwise blends derived from the two parent colors
color_env_bio = "#41a4ae"  # Env x Bio (blue-green blend)
color_env_spa = "#b0656a"  # Env x Spa (blue-red blend)
color_bio_spa = "#6b5f45"  # Bio x Spa (green-red blend)
color_shared  = "#9b59b6"  # All three (purple)

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
    lambda_f  = factor(lambda),
    alpha_lab = paste0("alpha = ", alpha)
  )

# Labels at the end (rightmost lambda) of each line
df_model_labels = df_model_lines %>%
  filter(lambda == max(lambda))

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
  facet_grid(spatial_form ~ alpha_lab) +
  labs(
    title = "Experiment 1: Model-level McFadden R\u00B2 across grid",
    subtitle = "Connected dots per anova Venn fraction",
    x = expression(lambda), y = "McFadden R\u00B2",
    color = "Component"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  ) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5.5, 40, 5.5, 5.5))

pdf(here("Calanda_JSDM", "plot", "exp1_anova_R2.pdf"),
    width = 14, height = 8)
print(p_anova)
dev.off()
cat("  Saved exp1_anova_R2.pdf\n")

cat("\n=== Post Experiment 1 analysis complete ===\n")
