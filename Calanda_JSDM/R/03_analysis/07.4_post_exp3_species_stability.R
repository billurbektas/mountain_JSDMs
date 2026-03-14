# ==============================================================================
# Script: 14_post_exp3_species_stability.R
# Purpose: Identify which species drive the instability in anova variance
#          partitioning. Same fitted model (sampling = 5000), 4 anova runs
#          (10k, 20k, 30k, 50k samples). Per-species variation across runs
#          is pure MC noise. Relate instability to species prevalence.
#
# Inputs:
#   - output/results/anova_saturation/exp3_model_fit5000.rds (for Y matrix)
#   - output/results/anova_saturation/exp3_anova_*.rds (4 VP runs)
#
# Outputs:
#   - plot/exp3_species_sd_vs_prevalence.pdf
#   - plot/exp3_species_range_vs_prevalence.pdf
#   - plot/exp3_species_pairwise_scatter.pdf
#   - plot/exp3_stable_vs_unstable_violins.pdf
#   - output/results/exp3_species_stability.csv
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
cat("\n=== Step 1: Loading data ===\n")

# Get Y matrix for prevalence (detect param_tag from model file)
model_files = list.files(
  here("Calanda_JSDM", "output", "results", "anova_saturation"),
  pattern = "^exp3_model_fit5000.*\\.rds$", full.names = TRUE
)
model = readRDS(model_files[1])
Y = model$data$Y
n_presences = colSums(Y)
n_sites = nrow(Y)
prevalence = n_presences / n_sites

# Load all anova runs
vp_files = list.files(
  here("Calanda_JSDM", "output", "results", "anova_saturation"),
  pattern = "^exp3_anova_.*\\.rds$", full.names = TRUE
)

species_arrays = list()  # list of 352 x 4 matrices
anova_sizes = c()

for (f in vp_files) {
  run = readRDS(f)
  ns = run$anova_samples
  anova_sizes = c(anova_sizes, ns)
  species_arrays[[as.character(ns)]] = run$partition$species  # 352 x 4
}

anova_sizes = sort(anova_sizes)
cat(sprintf("  %d species, %d anova runs (%s)\n",
            ncol(Y), length(anova_sizes), paste(anova_sizes, collapse = ", ")))

# ==============================================================================
# STEP 2: PER-SPECIES VARIABILITY ACROSS RUNS
# ==============================================================================
cat("\n=== Step 2: Computing per-species variability ===\n")

# For each component (env, spa, codist, r2), stack values across runs
components = c("env", "spa", "codist", "r2")
comp_labels = c(env = "Environment", spa = "Spatial", codist = "Biotic", r2 = "Overall")

stability_rows = list()

for (sp_idx in seq_len(ncol(Y))) {
  sp_name = colnames(Y)[sp_idx]

  for (comp in components) {
    values = sapply(as.character(anova_sizes), function(ns) {
      species_arrays[[ns]][sp_idx, comp]
    })

    stability_rows = c(stability_rows, list(tibble(
      species     = sp_name,
      species_idx = sp_idx,
      component   = comp_labels[comp],
      n_presences = n_presences[sp_idx],
      prevalence  = prevalence[sp_idx],
      mean_value  = mean(values),
      sd_value    = sd(values),
      range_value = diff(range(values)),
      cv_value    = ifelse(mean(values) != 0, sd(values) / abs(mean(values)), NA_real_),
      min_value   = min(values),
      max_value   = max(values)
    )))
  }
}

df_stability = bind_rows(stability_rows) %>%
  mutate(component = factor(component, levels = names(unit_colors)))

cat(sprintf("  Stability table: %d rows (%d species x %d components)\n",
            nrow(df_stability), n_distinct(df_stability$species), n_distinct(df_stability$component)))

# Summary stats
cat("\n  SD across runs (median per component):\n")
df_stability %>%
  group_by(component) %>%
  summarise(
    median_sd    = median(sd_value, na.rm = TRUE),
    mean_sd      = mean(sd_value, na.rm = TRUE),
    max_sd       = max(sd_value, na.rm = TRUE),
    median_range = median(range_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()

# ==============================================================================
# STEP 3: SD vs PREVALENCE SCATTER
# ==============================================================================
cat("\n=== Step 3: SD vs prevalence scatter ===\n")

p_sd_prev = ggplot(df_stability,
                   aes(x = n_presences, y = sd_value, color = component)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 0.8, alpha = 0.15) +
  scale_color_manual(values = unit_colors) +
  facet_wrap(~component, scales = "free_y") +
  labs(
    title = "Per-species SD of variance partitioning across 4 ANOVA runs",
    subtitle = "Same fitted model, anova samples = 10k / 20k / 30k / 50k",
    x = "Number of presences (out of 576 sites)",
    y = "SD across runs",
    color = "Component"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))

pdf(here("Calanda_JSDM", "plot", "detailed_modeling_experiment_figures", "exp3_species_sd_vs_prevalence.pdf"),
    width = 12, height = 8)
print(p_sd_prev)
dev.off()
cat("  Saved exp3_species_sd_vs_prevalence.pdf\n")

# ==============================================================================
# STEP 4: RANGE vs PREVALENCE SCATTER
# ==============================================================================
cat("\n=== Step 4: Range vs prevalence scatter ===\n")

p_range_prev = ggplot(df_stability,
                      aes(x = n_presences, y = range_value, color = component)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 0.8, alpha = 0.15) +
  scale_color_manual(values = unit_colors) +
  facet_wrap(~component, scales = "free_y") +
  labs(
    title = "Per-species range (max - min) of variance partitioning across 4 ANOVA runs",
    subtitle = "Same fitted model, anova samples = 10k / 20k / 30k / 50k",
    x = "Number of presences (out of 576 sites)",
    y = "Range across runs",
    color = "Component"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))

pdf(here("Calanda_JSDM", "plot", "detailed_modeling_experiment_figures", "exp3_species_range_vs_prevalence.pdf"),
    width = 12, height = 8)
print(p_range_prev)
dev.off()
cat("  Saved exp3_species_range_vs_prevalence.pdf\n")

# ==============================================================================
# STEP 5: PAIRWISE SCATTER (10k vs 50k) PER COMPONENT
# ==============================================================================
cat("\n=== Step 5: Pairwise scatter (10k vs 50k) ===\n")

min_ns = as.character(min(anova_sizes))
max_ns = as.character(max(anova_sizes))

df_pairwise = tibble(
  species     = colnames(Y),
  n_presences = as.integer(n_presences),
  env_lo      = species_arrays[[min_ns]][, "env"],
  env_hi      = species_arrays[[max_ns]][, "env"],
  spa_lo      = species_arrays[[min_ns]][, "spa"],
  spa_hi      = species_arrays[[max_ns]][, "spa"],
  codist_lo   = species_arrays[[min_ns]][, "codist"],
  codist_hi   = species_arrays[[max_ns]][, "codist"],
  r2_lo       = species_arrays[[min_ns]][, "r2"],
  r2_hi       = species_arrays[[max_ns]][, "r2"]
)

make_scatter = function(df, col_lo, col_hi, comp_label, comp_color) {
  r = cor(df[[col_lo]], df[[col_hi]], use = "complete.obs")
  ggplot(df, aes(x = .data[[col_lo]], y = .data[[col_hi]], color = n_presences)) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.3, linetype = "dashed", color = "grey50") +
    geom_point(alpha = 0.5, size = 1.5) +
    scale_color_viridis_c(option = "C", direction = -1) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
             label = sprintf("r = %.3f", r), size = 4, fontface = "bold") +
    labs(
      title = comp_label,
      x = paste0("ANOVA samples = ", format(as.integer(min_ns), big.mark = ",")),
      y = paste0("ANOVA samples = ", format(as.integer(max_ns), big.mark = ",")),
      color = "Presences"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "right")
}

p_env    = make_scatter(df_pairwise, "env_lo", "env_hi", "Environment", color_env)
p_spa    = make_scatter(df_pairwise, "spa_lo", "spa_hi", "Spatial", color_spa)
p_codist = make_scatter(df_pairwise, "codist_lo", "codist_hi", "Biotic", color_codist)
p_r2     = make_scatter(df_pairwise, "r2_lo", "r2_hi", "Overall R\u00B2", color_overall)

p_scatter = (p_r2 | p_env) / (p_spa | p_codist) +
  plot_annotation(
    title = sprintf("Experiment 3: Species-level partition (%s vs %s ANOVA samples)",
                    format(as.integer(min_ns), big.mark = ","),
                    format(as.integer(max_ns), big.mark = ",")),
    subtitle = "Color = number of presences (darker = rarer). Same fitted model.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10)
    )
  )

pdf(here("Calanda_JSDM", "plot", "detailed_modeling_experiment_figures", "exp3_species_pairwise_scatter.pdf"),
    width = 12, height = 10)
print(p_scatter)
dev.off()
cat("  Saved exp3_species_pairwise_scatter.pdf\n")

# ==============================================================================
# STEP 6: STABLE vs UNSTABLE SPECIES (by prevalence bins)
# ==============================================================================
cat("\n=== Step 6: Stable vs unstable species ===\n")

# Bin species by prevalence
df_stability = df_stability %>%
  mutate(
    prevalence_bin = cut(n_presences,
      breaks = c(0, 10, 30, 100, 200, Inf),
      labels = c("1-10", "11-30", "31-100", "101-200", ">200"),
      include.lowest = TRUE
    )
  )

# Count species per bin
cat("  Species per prevalence bin:\n")
df_stability %>%
  filter(component == "Overall") %>%
  count(prevalence_bin) %>%
  print()

# Violin of SD per prevalence bin
p_violin_sd = ggplot(df_stability,
                     aes(x = prevalence_bin, y = sd_value, fill = component)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  geom_violin(alpha = 0.3, scale = "width",
              position = position_dodge(width = 0.7),
              draw_quantiles = 0.5, linewidth = 0.2) +
  scale_fill_manual(values = unit_colors) +
  labs(
    title = "MC noise (SD across 4 ANOVA runs) by species prevalence",
    subtitle = "Same fitted model. Rarer species expected to be noisier",
    x = "Number of presences", y = "SD across runs",
    fill = "Component"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"))

# Violin of range per prevalence bin
p_violin_range = ggplot(df_stability,
                        aes(x = prevalence_bin, y = range_value, fill = component)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  geom_violin(alpha = 0.3, scale = "width",
              position = position_dodge(width = 0.7),
              draw_quantiles = 0.5, linewidth = 0.2) +
  scale_fill_manual(values = unit_colors) +
  labs(
    title = "MC noise (range across 4 ANOVA runs) by species prevalence",
    x = "Number of presences", y = "Range across runs",
    fill = "Component"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"))

p_combined = p_violin_sd / p_violin_range +
  plot_annotation(
    title = "Experiment 3: Species stability by prevalence",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

pdf(here("Calanda_JSDM", "plot", "detailed_modeling_experiment_figures", "exp3_stable_vs_unstable_violins.pdf"),
    width = 12, height = 12)
print(p_combined)
dev.off()
cat("  Saved exp3_stable_vs_unstable_violins.pdf\n")

# ==============================================================================
# STEP 7: CORRELATION BY PREVALENCE BIN
# ==============================================================================
cat("\n=== Step 7: Correlation by prevalence bin ===\n")

prev_bins = list(
  "1-10"    = which(n_presences >= 1  & n_presences <= 10),
  "11-30"   = which(n_presences >= 11 & n_presences <= 30),
  "31-100"  = which(n_presences >= 31 & n_presences <= 100),
  "101-200" = which(n_presences >= 101 & n_presences <= 200),
  ">200"    = which(n_presences > 200)
)

cor_by_bin = list()
for (bin_name in names(prev_bins)) {
  idx = prev_bins[[bin_name]]
  if (length(idx) < 3) next

  for (comp in components) {
    vals_lo = species_arrays[[min_ns]][idx, comp]
    vals_hi = species_arrays[[max_ns]][idx, comp]
    r = cor(vals_lo, vals_hi, use = "complete.obs")

    cor_by_bin = c(cor_by_bin, list(tibble(
      prevalence_bin = bin_name,
      component      = comp_labels[comp],
      n_species      = length(idx),
      correlation    = r
    )))
  }
}

df_cor_by_bin = bind_rows(cor_by_bin) %>%
  mutate(
    component = factor(component, levels = names(unit_colors)),
    prevalence_bin = factor(prevalence_bin, levels = names(prev_bins))
  )

cat("\n  Correlation (10k vs 50k) by prevalence bin:\n")
df_cor_by_bin %>%
  pivot_wider(names_from = component, values_from = correlation) %>%
  print()

p_cor_bin = ggplot(df_cor_by_bin,
                   aes(x = prevalence_bin, y = correlation, fill = component)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  geom_hline(yintercept = 1, linewidth = 0.3, linetype = "dashed", color = "grey60") +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.8) +
  geom_text(aes(label = sprintf("%.2f", correlation)),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 2.8) +
  scale_fill_manual(values = unit_colors) +
  labs(
    title = sprintf("Correlation between %s and %s ANOVA samples, by species prevalence",
                    format(as.integer(min_ns), big.mark = ","),
                    format(as.integer(max_ns), big.mark = ",")),
    subtitle = "Higher prevalence species expected to be more stable",
    x = "Number of presences", y = "Pearson r",
    fill = "Component"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

pdf(here("Calanda_JSDM", "plot", "detailed_modeling_experiment_figures", "exp3_correlation_by_prevalence.pdf"),
    width = 10, height = 6)
print(p_cor_bin)
dev.off()
cat("  Saved exp3_correlation_by_prevalence.pdf\n")

# ==============================================================================
# STEP 8: SAVE SUMMARY
# ==============================================================================
cat("\n=== Step 8: Saving summary ===\n")

write_csv(df_stability,
          here("Calanda_JSDM", "output", "results", "exp3_species_stability.csv"))
cat("  Saved exp3_species_stability.csv\n")

cat("\n=== Post Experiment 3 species stability analysis complete ===\n")
