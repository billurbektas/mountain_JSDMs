# ==============================================================================
# Script: 13_post_exp3_anova_saturation.R
# Purpose: Post-analysis of sjSDM Experiment 3 - ANOVA sampling saturation
#          Same fitted model (sampling = 5000), varying anova samples
#          (10k, 20k, 30k, 50k). Reference = 50k.
#          Isolates the effect of anova MC sampling from model fitting.
#
# Inputs:
#   - output/results/anova_saturation/exp3_anova_*_<param_tag>.rds
#   - param_tag format: a<alpha>_l<lambda>_le<lambda_env>
#
# Outputs:
#   - plot/exp3_anova_R2_<param_tag>.pdf
#   - plot/exp3_species_violins_<param_tag>.pdf
#   - plot/exp3_sites_violins_<param_tag>.pdf
#   - plot/exp3_species_cor_mad_<param_tag>.pdf
#   - plot/exp3_sites_cor_mad_<param_tag>.pdf
#   - plot/exp3_timing_<param_tag>.pdf
#   - output/results/exp3_anova_saturation_summary_<param_tag>.csv
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

dir.create(here("Calanda_JSDM", "plot", "detailed_modeling_experiment_figures"), showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# STEP 1: LOAD DATA
# ==============================================================================
cat("\n=== Step 1: Loading exp3 anova saturation results ===\n")

# Detect available param tags from anova files
all_anova_files = list.files(
  here("Calanda_JSDM", "output", "results", "anova_saturation"),
  pattern = "^exp3_anova_.*\\.rds$", full.names = TRUE
)

# Extract param_tag from first file (format: exp3_anova_<size>_<param_tag>.rds)
tag_match = regmatches(basename(all_anova_files[1]),
                       regexpr("a[0-9.]+_l[0-9.]+_le[0-9.]+(?=\\.)", basename(all_anova_files[1]), perl = TRUE))
if (length(tag_match) == 0 || tag_match == "") {
  stop("Could not detect param_tag from filenames. Expected format: exp3_anova_<size>_a<alpha>_l<lambda>_le<lambda_env>.rds")
}
param_tag = tag_match
cat(sprintf("  Detected param_tag: %s\n", param_tag))

# Filter to only files matching this param_tag
vp_files = all_anova_files[grepl(param_tag, basename(all_anova_files), fixed = TRUE)]
cat(sprintf("  Found %d VP files for %s\n", length(vp_files), param_tag))

exp3_data    = list()
metrics_list = list()
species_list = list()
sites_list   = list()

# Extract hyperparameters from first file
first_run = readRDS(vp_files[1])
hp_alpha      = first_run$alpha
hp_lambda     = first_run$lambda
hp_lambda_env = first_run$lambda_env
hp_subtitle   = sprintf("alpha = %s, lambda = %s, lambda_env = %s | fit sampling = 5000",
                         hp_alpha, hp_lambda, hp_lambda_env)
cat(sprintf("  Hyperparameters: %s\n", hp_subtitle))

for (f in vp_files) {
  run = readRDS(f)
  ns  = run$anova_samples
  ns_chr = as.character(ns)

  exp3_data[[ns_chr]] = run

  # Model-level McFadden R2 from anova results
  anova_res = run$partition$anova$results
  r2_mcf    = setNames(anova_res$`R2 McFadden`, anova_res$models)

  metrics_list[[ns_chr]] = tibble(
    anova_samples = ns,
    vp_time_min   = run$vp_time_min,
    metric_type   = c("overall_R2", "env_R2", "spa_R2", "bio_R2",
                      "env_bio_R2", "env_spa_R2", "bio_spa_R2", "shared_R2"),
    value         = c(r2_mcf["Full"], r2_mcf["F_A"], r2_mcf["F_S"], r2_mcf["F_B"],
                      r2_mcf["F_AB"], r2_mcf["F_AS"], r2_mcf["F_BS"], r2_mcf["F_ABS"])
  )

  # Species-level components (352 x 4)
  species_list[[ns_chr]] = as_tibble(run$partition$species) %>%
    mutate(anova_samples = ns, species_idx = row_number())

  # Site-level components (576 x 4)
  sites_list[[ns_chr]] = as_tibble(run$partition$sites) %>%
    mutate(anova_samples = ns, site_idx = row_number())
}

df_metrics = bind_rows(metrics_list)
df_species = bind_rows(species_list)
df_sites   = bind_rows(sites_list)

anova_sizes = sort(unique(df_species$anova_samples))
ref_size    = max(anova_sizes)

cat(sprintf("  ANOVA sample sizes: %s\n", paste(anova_sizes, collapse = ", ")))
cat(sprintf("  Reference size: %d\n", ref_size))
cat(sprintf("  Metrics: %d rows\n", nrow(df_metrics)))
cat(sprintf("  Species: %d rows (%d species x %d runs)\n",
            nrow(df_species), n_distinct(df_species$species_idx),
            n_distinct(df_species$anova_samples)))
cat(sprintf("  Sites:   %d rows (%d sites x %d runs)\n",
            nrow(df_sites), n_distinct(df_sites$site_idx),
            n_distinct(df_sites$anova_samples)))

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
  filter(anova_samples == ref_size)

p_anova = ggplot(df_model_lines,
                 aes(x = anova_samples, y = value, color = component, group = component)) +
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
    title = "Experiment 3: Model-level McFadden R\u00B2 across ANOVA sample sizes",
    subtitle = hp_subtitle,
    x = "ANOVA samples (log scale)", y = "McFadden R\u00B2",
    color = "Component"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  ) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5.5, 40, 5.5, 5.5))

pdf(here("Calanda_JSDM", "plot", "detailed_modeling_experiment_figures", paste0("exp3_anova_R2_", param_tag, ".pdf")),
    width = 12, height = 7)
print(p_anova)
dev.off()
cat(sprintf("  Saved exp3_anova_R2_%s.pdf\n", param_tag))

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
      component      = factor(component, levels = names(unit_colors)),
      anova_samples_f = factor(anova_samples, levels = anova_sizes,
                               labels = scales::comma(anova_sizes))
    )
}

make_violin_plot = function(unit_long, title_label) {
  ggplot(unit_long, aes(x = anova_samples_f, y = value, fill = component)) +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = "solid", color = "grey40") +
    geom_violin(
      alpha = 0.3, scale = "width",
      position = position_dodge(width = 0.7),
      draw_quantiles = 0.5, linewidth = 0.2
    ) +
    scale_fill_manual(values = unit_colors) +
    labs(
      title = sprintf("Experiment 3: %s-level variance partitioning across ANOVA sample sizes", title_label),
      subtitle = hp_subtitle,
      x = "ANOVA samples", y = "Value",
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

pdf(here("Calanda_JSDM", "plot", "detailed_modeling_experiment_figures", paste0("exp3_species_violins_", param_tag, ".pdf")),
    width = 12, height = 7)
print(p_species)
dev.off()
cat(sprintf("  Saved exp3_species_violins_%s.pdf\n", param_tag))

# ==============================================================================
# STEP 4: SITE VIOLIN PLOT
# ==============================================================================
cat("\n=== Step 4: Site violin plot ===\n")

df_sites_long = pivot_unit_long(df_sites)
p_sites = make_violin_plot(df_sites_long, "Site")

pdf(here("Calanda_JSDM", "plot", "detailed_modeling_experiment_figures", paste0("exp3_sites_violins_", param_tag, ".pdf")),
    width = 12, height = 7)
print(p_sites)
dev.off()
cat(sprintf("  Saved exp3_sites_violins_%s.pdf\n", param_tag))

# ==============================================================================
# STEP 5: CORRELATIONS + MAD vs REFERENCE
# ==============================================================================
cat("\n=== Step 5: Correlations and MAD vs reference ===\n")

ref_species = exp3_data[[as.character(ref_size)]]$partition$species
ref_sites   = exp3_data[[as.character(ref_size)]]$partition$sites

corr_mad_list = list()

for (ns in anova_sizes) {
  p_sp = exp3_data[[as.character(ns)]]$partition$species
  p_si = exp3_data[[as.character(ns)]]$partition$sites

  for (col in c("r2", "env", "spa", "codist")) {
    comp_label = recode(col, r2 = "Overall", env = "Environment",
                        spa = "Spatial", codist = "Biotic")

    corr_mad_list = c(corr_mad_list, list(tibble(
      anova_samples = ns,
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
  ggplot(aes(x = anova_samples, y = species_r, color = component)) +
  geom_hline(yintercept = 1, linewidth = 0.3, linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = 0.99, linewidth = 0.3, linetype = "dotted", color = "grey50") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(values = unit_colors) +
  labs(
    title = "A) Species-level correlation with reference",
    x = "ANOVA samples (log scale)", y = "Pearson r",
    color = "Component"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

p_sp_mad = df_corr_mad %>%
  ggplot(aes(x = anova_samples, y = species_mad, color = component)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(values = unit_colors) +
  labs(
    title = "B) Species-level MAD from reference",
    x = "ANOVA samples (log scale)", y = "Mean absolute difference",
    color = "Component"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

p_species_corr_mad = p_sp_cor / p_sp_mad +
  plot_annotation(
    title = sprintf("Experiment 3: Species-level ANOVA saturation (reference = %s samples)",
                    format(ref_size, big.mark = ",")),
    subtitle = hp_subtitle,
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10)
    )
  )

pdf(here("Calanda_JSDM", "plot", "detailed_modeling_experiment_figures", paste0("exp3_species_cor_mad_", param_tag, ".pdf")),
    width = 10, height = 10)
print(p_species_corr_mad)
dev.off()
cat(sprintf("  Saved exp3_species_cor_mad_%s.pdf\n", param_tag))

# --- Sites: correlation + MAD ---
p_si_cor = df_corr_mad %>%
  ggplot(aes(x = anova_samples, y = sites_r, color = component)) +
  geom_hline(yintercept = 1, linewidth = 0.3, linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = 0.99, linewidth = 0.3, linetype = "dotted", color = "grey50") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(values = unit_colors) +
  labs(
    title = "A) Site-level correlation with reference",
    x = "ANOVA samples (log scale)", y = "Pearson r",
    color = "Component"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

p_si_mad = df_corr_mad %>%
  ggplot(aes(x = anova_samples, y = sites_mad, color = component)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(values = unit_colors) +
  labs(
    title = "B) Site-level MAD from reference",
    x = "ANOVA samples (log scale)", y = "Mean absolute difference",
    color = "Component"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

p_sites_corr_mad = p_si_cor / p_si_mad +
  plot_annotation(
    title = sprintf("Experiment 3: Site-level ANOVA saturation (reference = %s samples)",
                    format(ref_size, big.mark = ",")),
    subtitle = hp_subtitle,
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10)
    )
  )

pdf(here("Calanda_JSDM", "plot", "detailed_modeling_experiment_figures", paste0("exp3_sites_cor_mad_", param_tag, ".pdf")),
    width = 10, height = 10)
print(p_sites_corr_mad)
dev.off()
cat(sprintf("  Saved exp3_sites_cor_mad_%s.pdf\n", param_tag))

# ==============================================================================
# STEP 5b: PER-SPECIES RAW VALUE TRAJECTORIES
# ==============================================================================
cat("\n=== Step 5b: Per-species raw value trajectories ===\n")

sp_raw_list = list()
for (ns in anova_sizes) {
  p_sp = exp3_data[[as.character(ns)]]$partition$species
  for (col in c("env", "spa", "codist", "r2")) {
    comp_label = recode(col, r2 = "Overall", env = "Environment",
                        spa = "Spatial", codist = "Biotic")
    sp_raw_list = c(sp_raw_list, list(tibble(
      anova_samples = ns, component = comp_label,
      idx = seq_len(nrow(p_sp)), value = p_sp[, col]
    )))
  }
}

df_sp_raw = bind_rows(sp_raw_list) %>%
  mutate(component = factor(component, levels = names(unit_colors)))

# Compute per-species abs difference between the two largest anova sizes
two_largest = sort(anova_sizes, decreasing = TRUE)[1:2]
val_at_max  = df_sp_raw %>% filter(anova_samples == two_largest[1]) %>% select(component, idx, v1 = value)
val_at_prev = df_sp_raw %>% filter(anova_samples == two_largest[2]) %>% select(component, idx, v2 = value)
df_sp_last_diff = val_at_max %>%
  left_join(val_at_prev, by = c("component", "idx")) %>%
  mutate(last_abs_diff = abs(v1 - v2)) %>%
  select(component, idx, last_abs_diff)

df_sp_raw = df_sp_raw %>%
  left_join(df_sp_last_diff, by = c("component", "idx"))

diff_label = paste0("|diff| last two\n(",
                    scales::comma(two_largest[2]), " vs ",
                    scales::comma(two_largest[1]), ")")

p_sp_raw = ggplot(df_sp_raw, aes(x = anova_samples, y = value,
                                  group = idx, color = last_abs_diff)) +
  geom_line(alpha = 0.4, linewidth = 0.3) +
  facet_wrap(~ component, ncol = 1) +
  scale_x_log10(labels = scales::comma) +
  scale_color_viridis_c(option = "cividis", direction = -1, name = diff_label) +
  labs(x = "ANOVA samples (log scale)", y = "VP value") +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(face = "bold"))

# Comparisons: last two runs + ref vs 130k, 100k, 80k
ref_sp   = exp3_data[[as.character(ref_size)]]$partition$species
sp_130k  = exp3_data[["130000"]]$partition$species
sp_100k  = exp3_data[["100000"]]$partition$species
sp_80k   = exp3_data[["80000"]]$partition$species

sp_hist_list = list()
for (col in c("env", "spa", "codist", "r2")) {
  comp_label = recode(col, r2 = "Overall", env = "Environment",
                      spa = "Spatial", codist = "Biotic")
  n_sp = nrow(ref_sp)
  sp_hist_list = c(sp_hist_list, list(
    tibble(component = comp_label, idx = seq_len(n_sp),
           comparison = paste0(scales::comma(two_largest[2]), " vs ", scales::comma(two_largest[1])),
           abs_diff = df_sp_last_diff %>% filter(component == comp_label) %>% pull(last_abs_diff)),
    tibble(component = comp_label, idx = seq_len(n_sp),
           comparison = paste0("130,000 vs ", scales::comma(ref_size)),
           abs_diff = abs(sp_130k[, col] - ref_sp[, col])),
    tibble(component = comp_label, idx = seq_len(n_sp),
           comparison = paste0("100,000 vs ", scales::comma(ref_size)),
           abs_diff = abs(sp_100k[, col] - ref_sp[, col])),
    tibble(component = comp_label, idx = seq_len(n_sp),
           comparison = paste0("80,000 vs ", scales::comma(ref_size)),
           abs_diff = abs(sp_80k[, col] - ref_sp[, col]))
  ))
}
df_sp_hist = bind_rows(sp_hist_list) %>%
  mutate(
    component = factor(component, levels = names(unit_colors)),
    comparison = factor(comparison, levels = c(
      paste0(scales::comma(two_largest[2]), " vs ", scales::comma(two_largest[1])),
      paste0("130,000 vs ", scales::comma(ref_size)),
      paste0("100,000 vs ", scales::comma(ref_size)),
      paste0("80,000 vs ", scales::comma(ref_size))
    ))
  )

# Pre-bin and compute mean color per bin
df_sp_hist = df_sp_hist %>%
  group_by(component, comparison) %>%
  mutate(bin = cut(abs_diff, breaks = 30, labels = FALSE, include.lowest = TRUE)) %>%
  ungroup()

df_sp_bin_summary = df_sp_hist %>%
  group_by(component, comparison, bin) %>%
  summarise(count = n(),
            mean_abs_diff = mean(abs_diff, na.rm = TRUE),
            .groups = "drop")

# Median + mean per facet
df_sp_hist_stats = df_sp_hist %>%
  group_by(component, comparison) %>%
  summarise(
    median_diff = median(abs_diff, na.rm = TRUE),
    mean_diff   = mean(abs_diff, na.rm = TRUE),
    .groups = "drop"
  )

p_sp_hist = ggplot(df_sp_bin_summary, aes(x = mean_abs_diff, y = count, fill = mean_abs_diff)) +
  geom_col(width = max(df_sp_hist$abs_diff, na.rm = TRUE) / 32, show.legend = FALSE) +
  geom_vline(data = df_sp_hist_stats, aes(xintercept = median_diff),
             linetype = "dashed", color = "grey30", linewidth = 0.5) +
  geom_vline(data = df_sp_hist_stats, aes(xintercept = mean_diff),
             linetype = "dotted", color = "grey30", linewidth = 0.5) +
  geom_text(data = df_sp_hist_stats, aes(x = median_diff, y = Inf,
            label = sprintf("med=%.4f", median_diff)),
            inherit.aes = FALSE, vjust = 1.5, hjust = -0.05, size = 2.2, color = "grey20") +
  geom_text(data = df_sp_hist_stats, aes(x = mean_diff, y = Inf,
            label = sprintf("mean=%.4f", mean_diff)),
            inherit.aes = FALSE, vjust = 3, hjust = -0.05, size = 2.2, color = "grey20") +
  facet_grid(component ~ comparison, scales = "free") +
  scale_fill_viridis_c(option = "cividis", direction = -1) +
  labs(x = "| diff |", y = "Count") +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(face = "bold"))

p_sp_combined = (p_sp_raw | p_sp_hist) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    title = "Experiment 3: Per-species VP values across ANOVA sample sizes",
    subtitle = hp_subtitle,
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10)
    )
  )

pdf(here("Calanda_JSDM", "plot", "detailed_modeling_experiment_figures", paste0("exp3_species_raw_trajectories_", param_tag, ".pdf")),
    width = 22, height = 12)
print(p_sp_combined)
dev.off()
cat(sprintf("  Saved exp3_species_raw_trajectories_%s.pdf\n", param_tag))

# ==============================================================================
# STEP 5c: PER-SITE RAW VALUE TRAJECTORIES
# ==============================================================================
cat("\n=== Step 5c: Per-site raw value trajectories ===\n")

si_raw_list = list()
for (ns in anova_sizes) {
  p_si = exp3_data[[as.character(ns)]]$partition$sites
  for (col in c("env", "spa", "codist", "r2")) {
    comp_label = recode(col, r2 = "Overall", env = "Environment",
                        spa = "Spatial", codist = "Biotic")
    si_raw_list = c(si_raw_list, list(tibble(
      anova_samples = ns, component = comp_label,
      idx = seq_len(nrow(p_si)), value = p_si[, col]
    )))
  }
}

df_si_raw = bind_rows(si_raw_list) %>%
  mutate(component = factor(component, levels = names(unit_colors)))

si_at_max  = df_si_raw %>% filter(anova_samples == two_largest[1]) %>% select(component, idx, v1 = value)
si_at_prev = df_si_raw %>% filter(anova_samples == two_largest[2]) %>% select(component, idx, v2 = value)
df_si_last_diff = si_at_max %>%
  left_join(si_at_prev, by = c("component", "idx")) %>%
  mutate(last_abs_diff = abs(v1 - v2)) %>%
  select(component, idx, last_abs_diff)

df_si_raw = df_si_raw %>%
  left_join(df_si_last_diff, by = c("component", "idx"))

p_si_raw = ggplot(df_si_raw, aes(x = anova_samples, y = value,
                                  group = idx, color = last_abs_diff)) +
  geom_line(alpha = 0.2, linewidth = 0.2) +
  facet_wrap(~ component, ncol = 1) +
  scale_x_log10(labels = scales::comma) +
  scale_color_viridis_c(option = "cividis", direction = -1, name = diff_label) +
  labs(x = "ANOVA samples (log scale)", y = "VP value") +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(face = "bold"))

# Comparisons: last two runs + ref vs 130k, 100k, 80k
ref_si   = exp3_data[[as.character(ref_size)]]$partition$sites
si_130k  = exp3_data[["130000"]]$partition$sites
si_100k  = exp3_data[["100000"]]$partition$sites
si_80k   = exp3_data[["80000"]]$partition$sites

si_hist_list = list()
for (col in c("env", "spa", "codist", "r2")) {
  comp_label = recode(col, r2 = "Overall", env = "Environment",
                      spa = "Spatial", codist = "Biotic")
  n_si = nrow(ref_si)
  si_hist_list = c(si_hist_list, list(
    tibble(component = comp_label, idx = seq_len(n_si),
           comparison = paste0(scales::comma(two_largest[2]), " vs ", scales::comma(two_largest[1])),
           abs_diff = df_si_last_diff %>% filter(component == comp_label) %>% pull(last_abs_diff)),
    tibble(component = comp_label, idx = seq_len(n_si),
           comparison = paste0("130,000 vs ", scales::comma(ref_size)),
           abs_diff = abs(si_130k[, col] - ref_si[, col])),
    tibble(component = comp_label, idx = seq_len(n_si),
           comparison = paste0("100,000 vs ", scales::comma(ref_size)),
           abs_diff = abs(si_100k[, col] - ref_si[, col])),
    tibble(component = comp_label, idx = seq_len(n_si),
           comparison = paste0("80,000 vs ", scales::comma(ref_size)),
           abs_diff = abs(si_80k[, col] - ref_si[, col]))
  ))
}
df_si_hist = bind_rows(si_hist_list) %>%
  mutate(
    component = factor(component, levels = names(unit_colors)),
    comparison = factor(comparison, levels = c(
      paste0(scales::comma(two_largest[2]), " vs ", scales::comma(two_largest[1])),
      paste0("130,000 vs ", scales::comma(ref_size)),
      paste0("100,000 vs ", scales::comma(ref_size)),
      paste0("80,000 vs ", scales::comma(ref_size))
    ))
  )

# Pre-bin and compute mean color per bin
df_si_hist = df_si_hist %>%
  group_by(component, comparison) %>%
  mutate(bin = cut(abs_diff, breaks = 30, labels = FALSE, include.lowest = TRUE)) %>%
  ungroup()

df_si_bin_summary = df_si_hist %>%
  group_by(component, comparison, bin) %>%
  summarise(count = n(),
            mean_abs_diff = mean(abs_diff, na.rm = TRUE),
            .groups = "drop")

# Median + mean per facet
df_si_hist_stats = df_si_hist %>%
  group_by(component, comparison) %>%
  summarise(
    median_diff = median(abs_diff, na.rm = TRUE),
    mean_diff   = mean(abs_diff, na.rm = TRUE),
    .groups = "drop"
  )

p_si_hist = ggplot(df_si_bin_summary, aes(x = mean_abs_diff, y = count, fill = mean_abs_diff)) +
  geom_col(width = max(df_si_hist$abs_diff, na.rm = TRUE) / 32, show.legend = FALSE) +
  geom_vline(data = df_si_hist_stats, aes(xintercept = median_diff),
             linetype = "dashed", color = "grey30", linewidth = 0.5) +
  geom_vline(data = df_si_hist_stats, aes(xintercept = mean_diff),
             linetype = "dotted", color = "grey30", linewidth = 0.5) +
  geom_text(data = df_si_hist_stats, aes(x = median_diff, y = Inf,
            label = sprintf("med=%.4f", median_diff)),
            inherit.aes = FALSE, vjust = 1.5, hjust = -0.05, size = 2.2, color = "grey20") +
  geom_text(data = df_si_hist_stats, aes(x = mean_diff, y = Inf,
            label = sprintf("mean=%.4f", mean_diff)),
            inherit.aes = FALSE, vjust = 3, hjust = -0.05, size = 2.2, color = "grey20") +
  facet_grid(component ~ comparison, scales = "free") +
  scale_fill_viridis_c(option = "cividis", direction = -1) +
  labs(x = "| diff |", y = "Count") +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(face = "bold"))

p_si_combined = (p_si_raw | p_si_hist) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    title = "Experiment 3: Per-site VP values across ANOVA sample sizes",
    subtitle = hp_subtitle,
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10)
    )
  )

pdf(here("Calanda_JSDM", "plot", "detailed_modeling_experiment_figures", paste0("exp3_sites_raw_trajectories_", param_tag, ".pdf")),
    width = 22, height = 12)
print(p_si_combined)
dev.off()
cat(sprintf("  Saved exp3_sites_raw_trajectories_%s.pdf\n", param_tag))

# ==============================================================================
# STEP 5d: PER-SPECIES ABSOLUTE DIFFERENCE FROM REFERENCE
# ==============================================================================
cat("\n=== Step 5d: Per-species absolute difference from reference ===\n")

ref_species_mat = exp3_data[[as.character(ref_size)]]$partition$species

sp_absdiff_list = list()
for (ns in anova_sizes) {
  p_sp = exp3_data[[as.character(ns)]]$partition$species
  for (col in c("env", "spa", "codist", "r2")) {
    comp_label = recode(col, r2 = "Overall", env = "Environment",
                        spa = "Spatial", codist = "Biotic")
    sp_absdiff_list = c(sp_absdiff_list, list(tibble(
      anova_samples = ns, component = comp_label,
      idx = seq_len(nrow(p_sp)),
      abs_diff = abs(p_sp[, col] - ref_species_mat[, col])
    )))
  }
}

df_sp_absdiff = bind_rows(sp_absdiff_list) %>%
  mutate(component = factor(component, levels = names(unit_colors)))

p_sp_absdiff = ggplot(df_sp_absdiff, aes(x = anova_samples, y = abs_diff, group = idx)) +
  geom_line(alpha = 0.3, linewidth = 0.25, color = "grey30") +
  facet_wrap(~ component, scales = "free_y") +
  scale_x_log10(labels = scales::comma) +
  labs(
    title = sprintf("Experiment 3: Per-species |difference| from reference (%s samples)",
                    format(ref_size, big.mark = ",")),
    subtitle = hp_subtitle,
    x = "ANOVA samples (log scale)", y = "| value - reference |"
  ) +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(face = "bold"))

pdf(here("Calanda_JSDM", "plot", "detailed_modeling_experiment_figures", paste0("exp3_species_absdiff_", param_tag, ".pdf")),
    width = 12, height = 9)
print(p_sp_absdiff)
dev.off()
cat(sprintf("  Saved exp3_species_absdiff_%s.pdf\n", param_tag))

# ==============================================================================
# STEP 5e: PER-SITE ABSOLUTE DIFFERENCE FROM REFERENCE
# ==============================================================================
cat("\n=== Step 5e: Per-site absolute difference from reference ===\n")

ref_sites_mat = exp3_data[[as.character(ref_size)]]$partition$sites

si_absdiff_list = list()
for (ns in anova_sizes) {
  p_si = exp3_data[[as.character(ns)]]$partition$sites
  for (col in c("env", "spa", "codist", "r2")) {
    comp_label = recode(col, r2 = "Overall", env = "Environment",
                        spa = "Spatial", codist = "Biotic")
    si_absdiff_list = c(si_absdiff_list, list(tibble(
      anova_samples = ns, component = comp_label,
      idx = seq_len(nrow(p_si)),
      abs_diff = abs(p_si[, col] - ref_sites_mat[, col])
    )))
  }
}

df_si_absdiff = bind_rows(si_absdiff_list) %>%
  mutate(component = factor(component, levels = names(unit_colors)))

p_si_absdiff = ggplot(df_si_absdiff, aes(x = anova_samples, y = abs_diff, group = idx)) +
  geom_line(alpha = 0.15, linewidth = 0.2, color = "grey30") +
  facet_wrap(~ component, scales = "free_y") +
  scale_x_log10(labels = scales::comma) +
  labs(
    title = sprintf("Experiment 3: Per-site |difference| from reference (%s samples)",
                    format(ref_size, big.mark = ",")),
    subtitle = hp_subtitle,
    x = "ANOVA samples (log scale)", y = "| value - reference |"
  ) +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(face = "bold"))

pdf(here("Calanda_JSDM", "plot", "detailed_modeling_experiment_figures", paste0("exp3_sites_absdiff_", param_tag, ".pdf")),
    width = 12, height = 9)
print(p_si_absdiff)
dev.off()
cat(sprintf("  Saved exp3_sites_absdiff_%s.pdf\n", param_tag))

# ==============================================================================
# STEP 6: TIMING PLOT
# ==============================================================================
cat("\n=== Step 6: Timing plot ===\n")

df_timing = df_metrics %>%
  distinct(anova_samples, vp_time_min)

p_timing = ggplot(df_timing, aes(x = anova_samples, y = vp_time_min)) +
  geom_line(linewidth = 0.8, color = "grey30") +
  geom_point(size = 3, color = "grey30") +
  geom_text(aes(label = paste0(round(vp_time_min, 1), " min")),
            vjust = -1, size = 3.5) +
  scale_x_log10(labels = scales::comma) +
  labs(
    title = "Experiment 3: Variance partitioning wall time",
    subtitle = paste0(hp_subtitle, "\nIncludes refitting 8 reduced models + MC evaluation"),
    x = "ANOVA samples (log scale)", y = "VP time (minutes)"
  ) +
  theme_bw(base_size = 11) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(10, 5.5, 5.5, 5.5))

pdf(here("Calanda_JSDM", "plot", "detailed_modeling_experiment_figures", paste0("exp3_timing_", param_tag, ".pdf")),
    width = 8, height = 5)
print(p_timing)
dev.off()
cat(sprintf("  Saved exp3_timing_%s.pdf\n", param_tag))

# ==============================================================================
# STEP 7: SUMMARY TABLE
# ==============================================================================
cat("\n=== Step 7: Summary table ===\n")

# Anova R2 per anova_samples (wide)
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
  distinct(anova_samples, vp_time_min)

# Correlations + MAD (wide)
df_corr_mad_wide = df_corr_mad %>%
  pivot_wider(
    names_from  = component,
    values_from = c(species_r, species_mad, sites_r, sites_mad),
    names_glue  = "{.value}_{component}"
  )

df_summary = df_anova_wide %>%
  left_join(df_timing_wide, by = "anova_samples") %>%
  left_join(df_corr_mad_wide, by = "anova_samples") %>%
  mutate(alpha = hp_alpha, lambda = hp_lambda, lambda_env = hp_lambda_env, .after = anova_samples)

summary_file = paste0("exp3_anova_saturation_summary_", param_tag, ".csv")
write_csv(df_summary,
          here("Calanda_JSDM", "output", "results", summary_file))
cat(sprintf("  Saved %s\n", summary_file))
print(df_summary)

cat("\n=== Post Experiment 3 analysis complete ===\n")
