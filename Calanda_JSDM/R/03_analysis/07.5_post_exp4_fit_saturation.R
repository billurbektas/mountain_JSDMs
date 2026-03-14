# ==============================================================================
# Script: 15_post_exp4_fit_saturation.R
# Purpose: Post-analysis of sjSDM Experiment 4 - Model fit sampling saturation
#          Fixed ANOVA samples (30k), varying model fit sampling
#          (5k, 10k, 20k, 30k). Reference = 30k.
#          Isolates the effect of model fitting MC sampling from anova sampling.
#
# Inputs:
#   - output/results/fit_saturation/exp4_anova_fit*_av30000_<param_tag>.rds
#   - param_tag format: a<alpha>_l<lambda>_le<lambda_env>
#
# Outputs:
#   - plot/exp4_anova_R2_<param_tag>.pdf
#   - plot/exp4_species_violins_<param_tag>.pdf
#   - plot/exp4_sites_violins_<param_tag>.pdf
#   - plot/exp4_species_cor_mad_<param_tag>.pdf
#   - plot/exp4_sites_cor_mad_<param_tag>.pdf
#   - plot/exp4_timing_<param_tag>.pdf
#   - output/results/exp4_fit_saturation_summary_<param_tag>.csv
#
# Requires: R/00_setup/functions_calanda.R
# ==============================================================================

library(sjSDM)
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

dir.create(here("Calanda_JSDM", "plot", "detailed_modeling_experiment_figures"), showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# STEP 1: LOAD DATA
# ==============================================================================
cat("\n=== Step 1: Loading exp4 fit saturation results ===\n")

# Detect available param tags from VP files
all_vp_files = list.files(
  here("Calanda_JSDM", "output", "results", "fit_saturation"),
  pattern = "^exp4_anova_fit.*\\.rds$", full.names = TRUE
)

# Extract param_tag from first file (format: exp4_anova_fit<size>_av<anova>_<param_tag>.rds)
tag_match = regmatches(basename(all_vp_files[1]),
                       regexpr("a[0-9.]+_l[0-9.]+_le[0-9.]+(?=\\.)", basename(all_vp_files[1]), perl = TRUE))
if (length(tag_match) == 0 || tag_match == "") {
  stop("Could not detect param_tag from filenames. Expected format: exp4_anova_fit<size>_av<anova>_a<alpha>_l<lambda>_le<lambda_env>.rds")
}
param_tag = tag_match
cat(sprintf("  Detected param_tag: %s\n", param_tag))

# Filter to only files matching this param_tag
vp_files = all_vp_files[grepl(param_tag, basename(all_vp_files), fixed = TRUE)]
cat(sprintf("  Found %d VP files for %s\n", length(vp_files), param_tag))

exp4_data    = list()
metrics_list = list()
species_list = list()
sites_list   = list()

# Extract hyperparameters from first file
first_run = readRDS(vp_files[1])
hp_alpha      = first_run$alpha
hp_lambda     = first_run$lambda
hp_lambda_env = first_run$lambda_env
fixed_anova   = first_run$anova_samples
hp_subtitle   = sprintf("alpha = %s, lambda = %s, lambda_env = %s | anova samples = %s",
                         hp_alpha, hp_lambda, hp_lambda_env, format(fixed_anova, big.mark = ","))
cat(sprintf("  Hyperparameters: %s\n", hp_subtitle))

for (f in vp_files) {
  run = readRDS(f)
  fs  = run$fit_sampling
  fs_chr = as.character(fs)

  exp4_data[[fs_chr]] = run

  # Model-level McFadden R2 from anova results
  anova_res = run$partition$anova$results
  r2_mcf    = setNames(anova_res$`R2 McFadden`, anova_res$models)

  metrics_list[[fs_chr]] = tibble(
    fit_sampling  = fs,
    vp_time_min   = run$vp_time_min,
    metric_type   = c("overall_R2", "env_R2", "spa_R2", "bio_R2",
                      "env_bio_R2", "env_spa_R2", "bio_spa_R2", "shared_R2"),
    value         = c(r2_mcf["Full"], r2_mcf["F_A"], r2_mcf["F_S"], r2_mcf["F_B"],
                      r2_mcf["F_AB"], r2_mcf["F_AS"], r2_mcf["F_BS"], r2_mcf["F_ABS"])
  )

  # Species-level components
  species_list[[fs_chr]] = as_tibble(run$partition$species) %>%
    mutate(fit_sampling = fs, species_idx = row_number())

  # Site-level components
  sites_list[[fs_chr]] = as_tibble(run$partition$sites) %>%
    mutate(fit_sampling = fs, site_idx = row_number())
}

df_metrics = bind_rows(metrics_list)
df_species = bind_rows(species_list)
df_sites   = bind_rows(sites_list)

fit_sizes = sort(unique(df_species$fit_sampling))
ref_size  = max(fit_sizes)

cat(sprintf("  Fit sampling sizes: %s\n", paste(fit_sizes, collapse = ", ")))
cat(sprintf("  Reference size: %d\n", ref_size))
cat(sprintf("  Metrics: %d rows\n", nrow(df_metrics)))
cat(sprintf("  Species: %d rows (%d species x %d runs)\n",
            nrow(df_species), n_distinct(df_species$species_idx),
            n_distinct(df_species$fit_sampling)))
cat(sprintf("  Sites:   %d rows (%d sites x %d runs)\n",
            nrow(df_sites), n_distinct(df_sites$site_idx),
            n_distinct(df_sites$fit_sampling)))

# ==============================================================================
# STEP 2: MODEL-LEVEL ANOVA R2 (connected dots)
# ==============================================================================
cat("\n=== Step 2: Model-level anova R2 plot ===\n")

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
  filter(fit_sampling == ref_size)

p_anova = ggplot(df_model_lines,
                 aes(x = fit_sampling, y = value, color = component, group = component)) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "solid", color = "grey40") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_text_repel(
    data = df_model_labels,
    aes(label = component),
    nudge_x = 0.1, direction = "y", hjust = 0,
    size = 3, segment.size = 0.3, show.legend = FALSE
  ) +
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(values = anova_colors) +
  labs(
    title = "Experiment 4: Model-level McFadden R\u00B2 across fit sampling sizes",
    subtitle = hp_subtitle,
    x = "Model fit sampling (log scale)", y = "McFadden R\u00B2",
    color = "Component"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  ) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5.5, 40, 5.5, 5.5))

pdf(here("Calanda_JSDM", "plot", "detailed_modeling_experiment_figures", paste0("exp4_anova_R2_", param_tag, ".pdf")),
    width = 12, height = 7)
print(p_anova)
dev.off()
cat(sprintf("  Saved exp4_anova_R2_%s.pdf\n", param_tag))

# ==============================================================================
# STEP 3: SPECIES VIOLIN PLOT
# ==============================================================================
cat("\n=== Step 3: Species violin plot ===\n")

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
      component       = factor(component, levels = names(unit_colors)),
      fit_sampling_f  = factor(fit_sampling, levels = fit_sizes,
                               labels = scales::comma(fit_sizes))
    )
}

make_violin_plot = function(unit_long, title_label) {
  ggplot(unit_long, aes(x = fit_sampling_f, y = value, fill = component)) +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = "solid", color = "grey40") +
    geom_violin(
      alpha = 0.3, scale = "width",
      position = position_dodge(width = 0.7),
      draw_quantiles = 0.5, linewidth = 0.2
    ) +
    scale_fill_manual(values = unit_colors) +
    labs(
      title = sprintf("Experiment 4: %s-level variance partitioning across fit sampling sizes", title_label),
      subtitle = hp_subtitle,
      x = "Model fit sampling", y = "Value",
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

pdf(here("Calanda_JSDM", "plot", "detailed_modeling_experiment_figures", paste0("exp4_species_violins_", param_tag, ".pdf")),
    width = 12, height = 7)
print(p_species)
dev.off()
cat(sprintf("  Saved exp4_species_violins_%s.pdf\n", param_tag))

# ==============================================================================
# STEP 4: SITE VIOLIN PLOT
# ==============================================================================
cat("\n=== Step 4: Site violin plot ===\n")

df_sites_long = pivot_unit_long(df_sites)
p_sites = make_violin_plot(df_sites_long, "Site")

pdf(here("Calanda_JSDM", "plot", "detailed_modeling_experiment_figures", paste0("exp4_sites_violins_", param_tag, ".pdf")),
    width = 12, height = 7)
print(p_sites)
dev.off()
cat(sprintf("  Saved exp4_sites_violins_%s.pdf\n", param_tag))

# ==============================================================================
# STEP 5: CORRELATIONS + MAD vs REFERENCE
# ==============================================================================
cat("\n=== Step 5: Correlations and MAD vs reference ===\n")

ref_species = exp4_data[[as.character(ref_size)]]$partition$species
ref_sites   = exp4_data[[as.character(ref_size)]]$partition$sites

corr_mad_list = list()

for (fs in fit_sizes) {
  p_sp = exp4_data[[as.character(fs)]]$partition$species
  p_si = exp4_data[[as.character(fs)]]$partition$sites

  for (col in c("r2", "env", "spa", "codist")) {
    comp_label = recode(col, r2 = "Overall", env = "Environment",
                        spa = "Spatial", codist = "Biotic")

    corr_mad_list = c(corr_mad_list, list(tibble(
      fit_sampling  = fs,
      component     = comp_label,
      species_r     = cor(p_sp[, col], ref_species[, col], use = "complete.obs"),
      species_mad   = mean(abs(p_sp[, col] - ref_species[, col]), na.rm = TRUE),
      sites_r       = cor(p_si[, col], ref_sites[, col], use = "complete.obs"),
      sites_mad     = mean(abs(p_si[, col] - ref_sites[, col]), na.rm = TRUE)
    )))
  }
}

df_corr_mad = bind_rows(corr_mad_list) %>%
  mutate(component = factor(component, levels = names(unit_colors)))

# --- Species: correlation + MAD ---
p_sp_cor = df_corr_mad %>%
  ggplot(aes(x = fit_sampling, y = species_r, color = component)) +
  geom_hline(yintercept = 1, linewidth = 0.3, linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = 0.99, linewidth = 0.3, linetype = "dotted", color = "grey50") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(values = unit_colors) +
  labs(
    title = "A) Species-level correlation with reference",
    x = "Model fit sampling (log scale)", y = "Pearson r",
    color = "Component"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

p_sp_mad = df_corr_mad %>%
  ggplot(aes(x = fit_sampling, y = species_mad, color = component)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(values = unit_colors) +
  labs(
    title = "B) Species-level MAD from reference",
    x = "Model fit sampling (log scale)", y = "Mean absolute difference",
    color = "Component"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

p_species_corr_mad = p_sp_cor / p_sp_mad +
  plot_annotation(
    title = sprintf("Experiment 4: Species-level fit saturation (reference = %s fit sampling)",
                    format(ref_size, big.mark = ",")),
    subtitle = hp_subtitle,
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10)
    )
  )

pdf(here("Calanda_JSDM", "plot", "detailed_modeling_experiment_figures", paste0("exp4_species_cor_mad_", param_tag, ".pdf")),
    width = 10, height = 10)
print(p_species_corr_mad)
dev.off()
cat(sprintf("  Saved exp4_species_cor_mad_%s.pdf\n", param_tag))

# --- Sites: correlation + MAD ---
p_si_cor = df_corr_mad %>%
  ggplot(aes(x = fit_sampling, y = sites_r, color = component)) +
  geom_hline(yintercept = 1, linewidth = 0.3, linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = 0.99, linewidth = 0.3, linetype = "dotted", color = "grey50") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(values = unit_colors) +
  labs(
    title = "A) Site-level correlation with reference",
    x = "Model fit sampling (log scale)", y = "Pearson r",
    color = "Component"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

p_si_mad = df_corr_mad %>%
  ggplot(aes(x = fit_sampling, y = sites_mad, color = component)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(values = unit_colors) +
  labs(
    title = "B) Site-level MAD from reference",
    x = "Model fit sampling (log scale)", y = "Mean absolute difference",
    color = "Component"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

p_sites_corr_mad = p_si_cor / p_si_mad +
  plot_annotation(
    title = sprintf("Experiment 4: Site-level fit saturation (reference = %s fit sampling)",
                    format(ref_size, big.mark = ",")),
    subtitle = hp_subtitle,
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10)
    )
  )

pdf(here("Calanda_JSDM", "plot", "detailed_modeling_experiment_figures", paste0("exp4_sites_cor_mad_", param_tag, ".pdf")),
    width = 10, height = 10)
print(p_sites_corr_mad)
dev.off()
cat(sprintf("  Saved exp4_sites_cor_mad_%s.pdf\n", param_tag))

# ==============================================================================
# STEP 6: TIMING PLOT
# ==============================================================================
cat("\n=== Step 6: Timing plot ===\n")

df_timing = df_metrics %>%
  distinct(fit_sampling, vp_time_min)

p_timing = ggplot(df_timing, aes(x = fit_sampling, y = vp_time_min)) +
  geom_line(linewidth = 0.8, color = "grey30") +
  geom_point(size = 3, color = "grey30") +
  geom_text(aes(label = paste0(round(vp_time_min, 1), " min")),
            vjust = -1, size = 3.5) +
  scale_x_log10(labels = scales::comma) +
  labs(
    title = "Experiment 4: Variance partitioning wall time",
    subtitle = paste0(hp_subtitle, "\nIncludes refitting 8 reduced models + MC evaluation"),
    x = "Model fit sampling (log scale)", y = "VP time (minutes)"
  ) +
  theme_bw(base_size = 11) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(10, 5.5, 5.5, 5.5))

pdf(here("Calanda_JSDM", "plot", "detailed_modeling_experiment_figures", paste0("exp4_timing_", param_tag, ".pdf")),
    width = 8, height = 5)
print(p_timing)
dev.off()
cat(sprintf("  Saved exp4_timing_%s.pdf\n", param_tag))

# ==============================================================================
# STEP 7: SUMMARY TABLE
# ==============================================================================
cat("\n=== Step 7: Summary table ===\n")

# Anova R2 per fit_sampling (wide)
df_anova_wide = df_metrics %>%
  mutate(
    metric_type = recode(metric_type,
      overall_R2 = "anova_Overall", env_R2 = "anova_Environment",
      spa_R2 = "anova_Spatial", bio_R2 = "anova_Biotic",
      env_bio_R2 = "anova_Env_x_Bio", env_spa_R2 = "anova_Env_x_Spa",
      bio_spa_R2 = "anova_Bio_x_Spa", shared_R2 = "anova_Shared"
    )
  ) %>%
  select(-vp_time_min) %>%
  pivot_wider(names_from = metric_type, values_from = value)

# Timing
df_timing_wide = df_metrics %>%
  distinct(fit_sampling, vp_time_min)

# Correlations + MAD (wide)
df_corr_mad_wide = df_corr_mad %>%
  pivot_wider(
    names_from  = component,
    values_from = c(species_r, species_mad, sites_r, sites_mad),
    names_glue  = "{.value}_{component}"
  )

df_summary = df_anova_wide %>%
  left_join(df_timing_wide, by = "fit_sampling") %>%
  left_join(df_corr_mad_wide, by = "fit_sampling") %>%
  mutate(
    anova_samples = fixed_anova,
    alpha = hp_alpha, lambda = hp_lambda, lambda_env = hp_lambda_env,
    .after = fit_sampling
  )

summary_file = paste0("exp4_fit_saturation_summary_", param_tag, ".csv")
write_csv(df_summary,
          here("Calanda_JSDM", "output", "results", summary_file))
cat(sprintf("  Saved %s\n", summary_file))
print(df_summary)

cat("\n=== Post Experiment 4 analysis complete ===\n")
