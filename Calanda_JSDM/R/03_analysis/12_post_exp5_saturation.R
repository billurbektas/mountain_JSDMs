# ==============================================================================
# Script: 12_post_exp5_saturation.R
# Purpose: Post-analysis of sjSDM Experiment 5 - sampling saturation
#          Same sampling value used for both model fitting (sampling_fit) and
#          variance partitioning (anova_samples). Largest run (10000) is the
#          reference. For each smaller run, compute:
#            - Model-level: anova R2 (8 Venn fractions) + convergence
#            - Species-level: correlations + MAD vs reference, violins
#            - Site-level: correlations + MAD vs reference, violins
#
# Inputs:
#   - output/results/saturation/saturation_*.rds (5 runs: 100, 500, 1000, 5000, 10000)
#
# Outputs:
#   - plot/exp5_convergence.pdf
#   - plot/exp5_anova_R2.pdf
#   - plot/exp5_species_violins.pdf
#   - plot/exp5_sites_violins.pdf
#   - plot/exp5_species_cor_mad.pdf
#   - plot/exp5_sites_cor_mad.pdf
#   - output/results/exp5_saturation_summary.csv
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

unit_colors = c(
  "Overall"     = color_overall,
  "Environment" = color_env,
  "Spatial"     = color_spa,
  "Biotic"      = color_codist
)

dir.create(here("Calanda_JSDM", "plot"), showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# STEP 1: LOAD DATA
# ==============================================================================
cat("\n=== Step 1: Loading saturation results ===\n")

sat_files = list.files(
  here("Calanda_JSDM", "output", "results", "saturation"),
  pattern = "^saturation_.*\\.rds$", full.names = TRUE
)
cat(sprintf("  Found %d saturation files\n", length(sat_files)))

# Parse sampling size from filename and load
sat_data = list()
history_list = list()
metrics_list = list()
species_list = list()
sites_list   = list()

for (f in sat_files) {
  run = readRDS(f)
  ns  = run$sampling
  ns_chr = as.character(ns)

  sat_data[[ns_chr]] = run

  # Training history
  history_list[[ns_chr]] = tibble(
    sampling  = ns,
    iteration = seq_along(run$model$history),
    loss      = as.numeric(run$model$history)
  )

  # Model-level McFadden R2 from anova results
  anova_res = run$partition$anova$results
  r2_mcf    = setNames(anova_res$`R2 McFadden`, anova_res$models)

  metrics_list[[ns_chr]] = tibble(
    sampling    = ns,
    metric_type = c("overall_R2", "env_R2", "spa_R2", "bio_R2",
                    "env_bio_R2", "env_spa_R2", "bio_spa_R2", "shared_R2"),
    value       = c(r2_mcf["Full"], r2_mcf["F_A"], r2_mcf["F_S"], r2_mcf["F_B"],
                    r2_mcf["F_AB"], r2_mcf["F_AS"], r2_mcf["F_BS"], r2_mcf["F_ABS"])
  )

  # Species-level components (352 x 4)
  species_list[[ns_chr]] = as_tibble(run$partition$species) %>%
    mutate(sampling = ns, species_idx = row_number())

  # Site-level components (576 x 4)
  sites_list[[ns_chr]] = as_tibble(run$partition$sites) %>%
    mutate(sampling = ns, site_idx = row_number())
}

df_history = bind_rows(history_list)
df_metrics = bind_rows(metrics_list)
df_species = bind_rows(species_list)
df_sites   = bind_rows(sites_list)

sample_sizes = sort(unique(df_species$sampling))
ref_size     = max(sample_sizes)

cat(sprintf("  Sample sizes: %s\n", paste(sample_sizes, collapse = ", ")))
cat(sprintf("  Reference size: %d\n", ref_size))
cat(sprintf("  History: %d rows\n", nrow(df_history)))
cat(sprintf("  Metrics: %d rows\n", nrow(df_metrics)))
cat(sprintf("  Species: %d rows (%d species x %d runs)\n",
            nrow(df_species), n_distinct(df_species$species_idx),
            n_distinct(df_species$sampling)))
cat(sprintf("  Sites:   %d rows (%d sites x %d runs)\n",
            nrow(df_sites), n_distinct(df_sites$site_idx),
            n_distinct(df_sites$sampling)))

# ==============================================================================
# STEP 2: CONVERGENCE CHECK
# ==============================================================================
cat("\n=== Step 2: Convergence plot ===\n")

p_conv = ggplot(df_history,
                aes(x = iteration, y = loss, color = factor(sampling))) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "solid", color = "grey40") +
  geom_line(alpha = 0.7, linewidth = 0.5) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Experiment 5: Training convergence across sampling sizes",
       x = "Iteration", y = "Loss",
       color = "Sampling") +
  theme_bw() +
  theme(legend.position = "bottom")

pdf(here("Calanda_JSDM", "plot", "exp5_convergence.pdf"),
    width = 10, height = 6)
print(p_conv)
dev.off()
cat("  Saved exp5_convergence.pdf\n")

# ==============================================================================
# STEP 3: MODEL-LEVEL ANOVA R2 (connected dots across sampling sizes)
# ==============================================================================
cat("\n=== Step 3: Model-level anova R2 plot ===\n")

df_model_lines = df_metrics %>%
  mutate(
    component = recode(metric_type,
      overall_R2 = "Overall", env_R2 = "Environment",
      spa_R2 = "Spatial", bio_R2 = "Biotic",
      env_bio_R2 = "Env x Bio", env_spa_R2 = "Env x Spa",
      bio_spa_R2 = "Bio x Spa", shared_R2 = "Shared"
    ),
    component = factor(component, levels = names(anova_colors))
  )

df_model_labels = df_model_lines %>%
  filter(sampling == ref_size)

p_anova = ggplot(df_model_lines,
                 aes(x = sampling, y = value, color = component, group = component)) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "solid", color = "grey40") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_text_repel(
    data = df_model_labels,
    aes(label = component),
    nudge_x = 0.15, direction = "y", hjust = 0,
    size = 3, segment.size = 0.3, show.legend = FALSE
  ) +
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(values = anova_colors) +
  labs(
    title = "Experiment 5: Model-level McFadden R\u00B2 across sampling sizes",
    subtitle = "Same sampling for model fitting and ANOVA variance partitioning",
    x = "Sampling size (log scale)", y = "McFadden R\u00B2",
    color = "Component"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  ) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5.5, 40, 5.5, 5.5))

pdf(here("Calanda_JSDM", "plot", "exp5_anova_R2.pdf"),
    width = 12, height = 7)
print(p_anova)
dev.off()
cat("  Saved exp5_anova_R2.pdf\n")

# ==============================================================================
# STEP 4: SPECIES VIOLIN PLOT
# ==============================================================================
cat("\n=== Step 4: Species violin plot ===\n")

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
      component  = factor(component, levels = names(unit_colors)),
      sampling_f = factor(sampling)
    )
}

make_violin_plot = function(unit_long, title_label) {
  ggplot(unit_long, aes(x = sampling_f, y = value, fill = component)) +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = "solid", color = "grey40") +
    geom_violin(
      alpha = 0.3, scale = "width",
      position = position_dodge(width = 0.7),
      draw_quantiles = 0.5, linewidth = 0.2
    ) +
    scale_fill_manual(values = unit_colors) +
    labs(
      title = sprintf("Experiment 5: %s-level variance partitioning across sampling sizes", title_label),
      subtitle = "Violins with median line",
      x = "Sampling size", y = "Value",
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

pdf(here("Calanda_JSDM", "plot", "exp5_species_violins.pdf"),
    width = 12, height = 7)
print(p_species)
dev.off()
cat("  Saved exp5_species_violins.pdf\n")

# ==============================================================================
# STEP 5: SITE VIOLIN PLOT
# ==============================================================================
cat("\n=== Step 5: Site violin plot ===\n")

df_sites_long = pivot_unit_long(df_sites)
p_sites = make_violin_plot(df_sites_long, "Site")

pdf(here("Calanda_JSDM", "plot", "exp5_sites_violins.pdf"),
    width = 12, height = 7)
print(p_sites)
dev.off()
cat("  Saved exp5_sites_violins.pdf\n")

# ==============================================================================
# STEP 6: CORRELATIONS + MAD vs REFERENCE
# ==============================================================================
cat("\n=== Step 6: Correlations and MAD vs reference ===\n")

ref_species = sat_data[[as.character(ref_size)]]$partition$species
ref_sites   = sat_data[[as.character(ref_size)]]$partition$sites

corr_mad_list = list()

for (ns in sample_sizes) {
  p_sp = sat_data[[as.character(ns)]]$partition$species
  p_si = sat_data[[as.character(ns)]]$partition$sites

  for (col in c("r2", "env", "spa", "codist")) {
    comp_label = recode(col, r2 = "Overall", env = "Environment",
                        spa = "Spatial", codist = "Biotic")

    corr_mad_list = c(corr_mad_list, list(tibble(
      sampling    = ns,
      component   = comp_label,
      species_r   = cor(p_sp[, col], ref_species[, col], use = "complete.obs"),
      species_mad = mean(abs(p_sp[, col] - ref_species[, col]), na.rm = TRUE),
      sites_r     = cor(p_si[, col], ref_sites[, col], use = "complete.obs"),
      sites_mad   = mean(abs(p_si[, col] - ref_sites[, col]), na.rm = TRUE)
    )))
  }
}

df_corr_mad = bind_rows(corr_mad_list) %>%
  mutate(component = factor(component, levels = names(unit_colors)))

# --- Species: correlation + MAD ---
df_sp_long = df_corr_mad %>%
  select(sampling, component, species_r, species_mad) %>%
  pivot_longer(cols = c(species_r, species_mad),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
    species_r   = "Correlation (Pearson r)",
    species_mad = "Mean absolute difference"
  ))

p_sp_cor = df_sp_long %>%
  filter(metric == "Correlation (Pearson r)") %>%
  ggplot(aes(x = sampling, y = value, color = component)) +
  geom_hline(yintercept = 1, linewidth = 0.3, linetype = "dashed", color = "grey60") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(values = unit_colors) +
  labs(
    title = "A) Species-level correlation with reference",
    x = "Sampling size (log scale)", y = "Pearson r",
    color = "Component"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

p_sp_mad = df_sp_long %>%
  filter(metric == "Mean absolute difference") %>%
  ggplot(aes(x = sampling, y = value, color = component)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(values = unit_colors) +
  labs(
    title = "B) Species-level MAD from reference",
    x = "Sampling size (log scale)", y = "Mean absolute difference",
    color = "Component"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

p_species_corr_mad = p_sp_cor / p_sp_mad +
  plot_annotation(
    title = sprintf("Experiment 5: Species-level saturation (reference = %s)", format(ref_size, big.mark = ",")),
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

pdf(here("Calanda_JSDM", "plot", "exp5_species_cor_mad.pdf"),
    width = 10, height = 10)
print(p_species_corr_mad)
dev.off()
cat("  Saved exp5_species_cor_mad.pdf\n")

# --- Sites: correlation + MAD ---
df_si_long = df_corr_mad %>%
  select(sampling, component, sites_r, sites_mad) %>%
  pivot_longer(cols = c(sites_r, sites_mad),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
    sites_r   = "Correlation (Pearson r)",
    sites_mad = "Mean absolute difference"
  ))

p_si_cor = df_si_long %>%
  filter(metric == "Correlation (Pearson r)") %>%
  ggplot(aes(x = sampling, y = value, color = component)) +
  geom_hline(yintercept = 1, linewidth = 0.3, linetype = "dashed", color = "grey60") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(values = unit_colors) +
  labs(
    title = "A) Site-level correlation with reference",
    x = "Sampling size (log scale)", y = "Pearson r",
    color = "Component"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

p_si_mad = df_si_long %>%
  filter(metric == "Mean absolute difference") %>%
  ggplot(aes(x = sampling, y = value, color = component)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(values = unit_colors) +
  labs(
    title = "B) Site-level MAD from reference",
    x = "Sampling size (log scale)", y = "Mean absolute difference",
    color = "Component"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

p_sites_corr_mad = p_si_cor / p_si_mad +
  plot_annotation(
    title = sprintf("Experiment 5: Site-level saturation (reference = %s)", format(ref_size, big.mark = ",")),
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

pdf(here("Calanda_JSDM", "plot", "exp5_sites_cor_mad.pdf"),
    width = 10, height = 10)
print(p_sites_corr_mad)
dev.off()
cat("  Saved exp5_sites_cor_mad.pdf\n")

# ==============================================================================
# STEP 7: SUMMARY TABLE
# ==============================================================================
cat("\n=== Step 7: Summary table ===\n")

# Anova R2 per sampling size (wide)
df_anova_wide = df_metrics %>%
  mutate(
    metric_type = recode(metric_type,
      overall_R2 = "anova_Overall", env_R2 = "anova_Environment",
      spa_R2 = "anova_Spatial", bio_R2 = "anova_Biotic",
      env_bio_R2 = "anova_Env_x_Bio", env_spa_R2 = "anova_Env_x_Spa",
      bio_spa_R2 = "anova_Bio_x_Spa", shared_R2 = "anova_Shared"
    )
  ) %>%
  pivot_wider(names_from = metric_type, values_from = value)

# Correlations + MAD (wide)
df_corr_mad_wide = df_corr_mad %>%
  pivot_wider(
    names_from  = component,
    values_from = c(species_r, species_mad, sites_r, sites_mad),
    names_glue  = "{.value}_{component}"
  )

df_summary = df_anova_wide %>%
  left_join(df_corr_mad_wide, by = "sampling")

write_csv(df_summary,
          here("Calanda_JSDM", "output", "results", "exp5_saturation_summary.csv"))
cat("  Saved exp5_saturation_summary.csv\n")
print(df_summary)

cat("\n=== Post Experiment 5 analysis complete ===\n")
