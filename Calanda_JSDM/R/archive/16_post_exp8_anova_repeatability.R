# ==============================================================================
# Script: 16_post_exp8_anova_repeatability.R
# Purpose: Post-analysis of sjSDM Experiment 8 - ANOVA repeatability
#          Same fitted model (sampling = 5000), same anova samples (80k),
#          10 repetitions with different seeds.
#          Goal: species-level confidence intervals for variance partitioning.
#
# Inputs:
#   - output/results/anova_repeatability/exp8_anova_rep*_av80000_<param_tag>.rds
#
# Outputs:
#   - plot/exp8_model_R2_reps_<param_tag>.pdf
#   - plot/exp8_species_sd_hist_<param_tag>.pdf
#   - plot/exp8_species_ci_ranked_<param_tag>.pdf
#   - plot/exp8_sites_sd_hist_<param_tag>.pdf
#   - plot/exp8_sites_ci_ranked_<param_tag>.pdf
#   - plot/exp8_species_pairwise_cor_<param_tag>.pdf
#   - output/results/exp8_species_ci_summary_<param_tag>.csv
#   - output/results/exp8_sites_ci_summary_<param_tag>.csv
#
# Requires: R/00_setup/functions_calanda.R
# ==============================================================================

library(sjSDM)
library(tidyverse)
library(conflicted)
library(here)
library(patchwork)

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::select)

source(here("Calanda_JSDM", "R", "00_setup", "functions_calanda.R"))

# Colors (consistent with project conventions)
color_env     = "#81caf3"
color_spa     = "#d00000"
color_codist  = "#00bd89"
color_overall = "black"

unit_colors = c(
  "Overall"     = color_overall,
  "Environment" = color_env,
  "Spatial"     = color_spa,
  "Biotic"      = color_codist
)

comp_recode = c(r2 = "Overall", env = "Environment", spa = "Spatial", codist = "Biotic")

dir.create(here("Calanda_JSDM", "plot"), showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# STEP 1: LOAD DATA
# ==============================================================================
cat("\n=== Step 1: Loading exp8 anova repeatability results ===\n")

all_vp_files = list.files(
  here("Calanda_JSDM", "output", "results", "anova_repeatability"),
  pattern = "^exp8_anova_rep.*\\.rds$", full.names = TRUE
)

# Extract param_tag
tag_match = regmatches(basename(all_vp_files[1]),
                       regexpr("a[0-9.]+_l[0-9.]+_le[0-9.]+(?=\\.)", basename(all_vp_files[1]), perl = TRUE))
if (length(tag_match) == 0 || tag_match == "") {
  stop("Could not detect param_tag from filenames.")
}
param_tag = tag_match
cat(sprintf("  Detected param_tag: %s\n", param_tag))

vp_files = all_vp_files[grepl(param_tag, basename(all_vp_files), fixed = TRUE)]
cat(sprintf("  Found %d rep files for %s\n", length(vp_files), param_tag))

# Load all reps
exp8_data    = list()
metrics_list = list()
species_list = list()
sites_list   = list()

first_run = readRDS(vp_files[1])
hp_alpha      = first_run$alpha
hp_lambda     = first_run$lambda
hp_lambda_env = first_run$lambda_env
fixed_anova   = first_run$anova_samples
fit_sampling  = first_run$fit_sampling
hp_subtitle   = sprintf("alpha = %s, lambda = %s, lambda_env = %s | fit = %s, anova = %s",
                         hp_alpha, hp_lambda, hp_lambda_env,
                         format(fit_sampling, big.mark = ","),
                         format(fixed_anova, big.mark = ","))
cat(sprintf("  Hyperparameters: %s\n", hp_subtitle))

for (f in vp_files) {
  run    = readRDS(f)
  rep_i  = run$rep
  ri_chr = as.character(rep_i)

  exp8_data[[ri_chr]] = run

  # Model-level McFadden R2
  anova_res = run$partition$anova$results
  r2_mcf    = setNames(anova_res$`R2 McFadden`, anova_res$models)

  metrics_list[[ri_chr]] = tibble(
    rep        = rep_i,
    seed       = run$seed,
    vp_time_min = run$vp_time_min,
    metric_type = c("overall_R2", "env_R2", "spa_R2", "bio_R2",
                    "env_bio_R2", "env_spa_R2", "bio_spa_R2", "shared_R2"),
    value       = c(r2_mcf["Full"], r2_mcf["F_A"], r2_mcf["F_S"], r2_mcf["F_B"],
                    r2_mcf["F_AB"], r2_mcf["F_AS"], r2_mcf["F_BS"], r2_mcf["F_ABS"])
  )

  # Species-level
  species_list[[ri_chr]] = as_tibble(run$partition$species) %>%
    mutate(rep = rep_i, species_idx = row_number())

  # Site-level
  sites_list[[ri_chr]] = as_tibble(run$partition$sites) %>%
    mutate(rep = rep_i, site_idx = row_number())
}

df_metrics = bind_rows(metrics_list)
df_species = bind_rows(species_list)
df_sites   = bind_rows(sites_list)

n_reps    = n_distinct(df_species$rep)
n_species = n_distinct(df_species$species_idx)
n_sites   = n_distinct(df_sites$site_idx)

cat(sprintf("  Reps: %d | Species: %d | Sites: %d\n", n_reps, n_species, n_sites))

# ==============================================================================
# STEP 2: MODEL-LEVEL R2 ACROSS REPS (dot + range)
# ==============================================================================
cat("\n=== Step 2: Model-level R2 across reps ===\n")

df_model = df_metrics %>%
  mutate(
    component = recode(metric_type,
      overall_R2 = "Overall", env_R2 = "Environment",
      spa_R2 = "Spatial", bio_R2 = "Biotic",
      env_bio_R2 = "Env x Bio", env_spa_R2 = "Env x Spa",
      bio_spa_R2 = "Bio x Spa", shared_R2 = "Shared"
    )
  ) %>%
  filter(component %in% c("Overall", "Environment", "Spatial", "Biotic"))

df_model_summary = df_model %>%
  group_by(component) %>%
  summarise(
    mean = mean(value), sd = sd(value),
    min = min(value), max = max(value),
    .groups = "drop"
  ) %>%
  mutate(component = factor(component, levels = names(unit_colors)))

p_model_reps = ggplot(df_model %>%
                        mutate(component = factor(component, levels = names(unit_colors))),
                      aes(x = component, y = value, color = component)) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  geom_crossbar(
    data = df_model_summary,
    aes(x = component, y = mean, ymin = mean - sd, ymax = mean + sd),
    width = 0.4, linewidth = 0.5, fatten = 2, alpha = 0.3
  ) +
  scale_color_manual(values = unit_colors) +
  labs(
    title = sprintf("Experiment 8: Model-level McFadden R\u00B2 across %d reps", n_reps),
    subtitle = hp_subtitle,
    x = NULL, y = "McFadden R\u00B2",
    caption = "Crossbar = mean \u00B1 1 SD"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

pdf(here("Calanda_JSDM", "plot", paste0("exp8_model_R2_reps_", param_tag, ".pdf")),
    width = 8, height = 6)
print(p_model_reps)
dev.off()
cat(sprintf("  Saved exp8_model_R2_reps_%s.pdf\n", param_tag))

# ==============================================================================
# STEP 3: SPECIES-LEVEL SD HISTOGRAMS
# ==============================================================================
cat("\n=== Step 3: Species-level SD histograms ===\n")

df_species_stats = df_species %>%
  pivot_longer(cols = c(env, spa, codist, r2),
               names_to = "component_raw", values_to = "value") %>%
  mutate(component = recode(component_raw, !!!comp_recode),
         component = factor(component, levels = names(unit_colors))) %>%
  group_by(species_idx, component) %>%
  summarise(
    mean = mean(value), sd = sd(value),
    lo = quantile(value, 0.025), hi = quantile(value, 0.975),
    range = max(value) - min(value),
    .groups = "drop"
  )

p_sp_sd = ggplot(df_species_stats, aes(x = sd, fill = component)) +
  geom_histogram(bins = 40, alpha = 0.7, position = "identity") +
  facet_wrap(~ component, scales = "free_y") +
  scale_fill_manual(values = unit_colors) +
  labs(
    title = sprintf("Experiment 8: Distribution of species-level SD across %d reps", n_reps),
    subtitle = hp_subtitle,
    x = "Standard deviation (across reps)", y = "Number of species"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

pdf(here("Calanda_JSDM", "plot", paste0("exp8_species_sd_hist_", param_tag, ".pdf")),
    width = 10, height = 7)
print(p_sp_sd)
dev.off()
cat(sprintf("  Saved exp8_species_sd_hist_%s.pdf\n", param_tag))

# ==============================================================================
# STEP 4: SPECIES-LEVEL RANKED CI PLOT (mean +/- 95% interval)
# ==============================================================================
cat("\n=== Step 4: Species-level ranked CI plot ===\n")

make_ranked_ci_plot = function(stats_df, unit_label) {
  plots = list()
  for (comp in names(unit_colors)) {
    d = stats_df %>%
      filter(component == comp) %>%
      arrange(mean) %>%
      mutate(rank = row_number())

    plots[[comp]] = ggplot(d, aes(x = rank, y = mean)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), fill = unit_colors[comp], alpha = 0.25) +
      geom_line(color = unit_colors[comp], linewidth = 0.4) +
      geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
      labs(
        title = comp,
        x = paste0(unit_label, " (ranked by mean)"),
        y = "Value (mean \u00B1 95% interval)"
      ) +
      theme_bw(base_size = 10)
  }

  (plots[["Environment"]] + plots[["Spatial"]]) /
    (plots[["Biotic"]] + plots[["Overall"]]) +
    plot_annotation(
      title = sprintf("Experiment 8: %s-level variance partitioning (mean \u00B1 95%% CI, %d reps)",
                      unit_label, n_reps),
      subtitle = hp_subtitle,
      theme = theme(
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10)
      )
    )
}

p_sp_ci = make_ranked_ci_plot(df_species_stats, "Species")

pdf(here("Calanda_JSDM", "plot", paste0("exp8_species_ci_ranked_", param_tag, ".pdf")),
    width = 12, height = 9)
print(p_sp_ci)
dev.off()
cat(sprintf("  Saved exp8_species_ci_ranked_%s.pdf\n", param_tag))

# ==============================================================================
# STEP 5: SITE-LEVEL SD HISTOGRAMS
# ==============================================================================
cat("\n=== Step 5: Site-level SD histograms ===\n")

df_sites_stats = df_sites %>%
  pivot_longer(cols = c(env, spa, codist, r2),
               names_to = "component_raw", values_to = "value") %>%
  mutate(component = recode(component_raw, !!!comp_recode),
         component = factor(component, levels = names(unit_colors))) %>%
  group_by(site_idx, component) %>%
  summarise(
    mean = mean(value), sd = sd(value),
    lo = quantile(value, 0.025), hi = quantile(value, 0.975),
    range = max(value) - min(value),
    .groups = "drop"
  )

p_si_sd = ggplot(df_sites_stats, aes(x = sd, fill = component)) +
  geom_histogram(bins = 40, alpha = 0.7, position = "identity") +
  facet_wrap(~ component, scales = "free_y") +
  scale_fill_manual(values = unit_colors) +
  labs(
    title = sprintf("Experiment 8: Distribution of site-level SD across %d reps", n_reps),
    subtitle = hp_subtitle,
    x = "Standard deviation (across reps)", y = "Number of sites"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

pdf(here("Calanda_JSDM", "plot", paste0("exp8_sites_sd_hist_", param_tag, ".pdf")),
    width = 10, height = 7)
print(p_si_sd)
dev.off()
cat(sprintf("  Saved exp8_sites_sd_hist_%s.pdf\n", param_tag))

# ==============================================================================
# STEP 6: SITE-LEVEL RANKED CI PLOT
# ==============================================================================
cat("\n=== Step 6: Site-level ranked CI plot ===\n")

p_si_ci = make_ranked_ci_plot(df_sites_stats, "Site")

pdf(here("Calanda_JSDM", "plot", paste0("exp8_sites_ci_ranked_", param_tag, ".pdf")),
    width = 12, height = 9)
print(p_si_ci)
dev.off()
cat(sprintf("  Saved exp8_sites_ci_ranked_%s.pdf\n", param_tag))

# ==============================================================================
# STEP 7: PAIRWISE CORRELATION ACROSS REPS (species-level)
# ==============================================================================
cat("\n=== Step 7: Pairwise rep correlation heatmap (species-level) ===\n")

# For each component, compute pairwise correlation matrix across reps
rep_ids = sort(unique(df_species$rep))

cor_list = list()
for (comp_raw in c("env", "spa", "codist", "r2")) {
  comp_label = comp_recode[comp_raw]

  # Build matrix: species x reps
  mat = matrix(NA, nrow = n_species, ncol = n_reps)
  for (j in seq_along(rep_ids)) {
    vals = df_species %>%
      filter(rep == rep_ids[j]) %>%
      arrange(species_idx) %>%
      pull(!!sym(comp_raw))
    mat[, j] = vals
  }
  colnames(mat) = paste0("rep", rep_ids)

  cor_mat = cor(mat, use = "complete.obs")
  cor_long = as_tibble(as.data.frame(as.table(cor_mat))) %>%
    rename(rep_x = Var1, rep_y = Var2, r = Freq) %>%
    mutate(component = comp_label)

  cor_list[[comp_raw]] = cor_long
}

df_cor_pairs = bind_rows(cor_list) %>%
  mutate(component = factor(component, levels = names(unit_colors)))

p_cor_heat = ggplot(df_cor_pairs, aes(x = rep_x, y = rep_y, fill = r)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.3f", r)), size = 2.2) +
  facet_wrap(~ component, nrow = 2) +
  scale_fill_gradient2(low = "#d73027", mid = "white", high = "#1a9850",
                       midpoint = 1, limits = c(NA, 1)) +
  labs(
    title = sprintf("Experiment 8: Pairwise correlation of species-level VP across %d reps", n_reps),
    subtitle = hp_subtitle,
    x = NULL, y = NULL, fill = "Pearson r"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )

pdf(here("Calanda_JSDM", "plot", paste0("exp8_species_pairwise_cor_", param_tag, ".pdf")),
    width = 10, height = 9)
print(p_cor_heat)
dev.off()
cat(sprintf("  Saved exp8_species_pairwise_cor_%s.pdf\n", param_tag))

# ==============================================================================
# STEP 8: SUMMARY TABLES
# ==============================================================================
cat("\n=== Step 8: Summary tables ===\n")

# Species summary
write_csv(df_species_stats %>%
            mutate(alpha = hp_alpha, lambda = hp_lambda, lambda_env = hp_lambda_env,
                   anova_samples = fixed_anova, fit_sampling = fit_sampling,
                   n_reps = n_reps, .before = 1),
          here("Calanda_JSDM", "output", "results",
               paste0("exp8_species_ci_summary_", param_tag, ".csv")))

# Sites summary
write_csv(df_sites_stats %>%
            mutate(alpha = hp_alpha, lambda = hp_lambda, lambda_env = hp_lambda_env,
                   anova_samples = fixed_anova, fit_sampling = fit_sampling,
                   n_reps = n_reps, .before = 1),
          here("Calanda_JSDM", "output", "results",
               paste0("exp8_sites_ci_summary_", param_tag, ".csv")))

# Model-level summary across reps
df_model_rep_summary = df_model %>%
  group_by(metric_type) %>%
  summarise(
    mean = mean(value), sd = sd(value),
    min = min(value), max = max(value),
    cv = sd(value) / abs(mean(value)),
    .groups = "drop"
  ) %>%
  mutate(alpha = hp_alpha, lambda = hp_lambda, lambda_env = hp_lambda_env,
         anova_samples = fixed_anova, fit_sampling = fit_sampling,
         n_reps = n_reps, .before = 1)

write_csv(df_model_rep_summary,
          here("Calanda_JSDM", "output", "results",
               paste0("exp8_model_summary_", param_tag, ".csv")))

cat("  Saved exp8_species_ci_summary, exp8_sites_ci_summary, exp8_model_summary CSVs\n")

# Print quick overview
cat("\n--- Model-level R2 stability ---\n")
print(df_model_rep_summary %>% select(metric_type, mean, sd, cv))

cat("\n--- Species-level SD summary (across components) ---\n")
df_species_stats %>%
  group_by(component) %>%
  summarise(
    mean_sd = mean(sd, na.rm = TRUE),
    median_sd = median(sd, na.rm = TRUE),
    max_sd = max(sd, na.rm = TRUE),
    pct_sd_gt_01 = mean(sd > 0.01, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  print()

cat("\n=== Post Experiment 8 analysis complete ===\n")
