# ==============================================================================
# Script: 08_variance_partitioning.R
# Purpose: Variance partitioning from 10-fold CV results.
#          - Venn diagram: median McFadden R² anova components + CI table
#          - Ternary plots: species (filtered by AUC >= 0.7) and sites
#            (filtered by logloss <= log(2)), colored by altitude
#          - Summary dataframes with mean, CI, weight = 1/CI_width²
#          - Reports both McFadden R² (anova) and Tjur R² (discrimination)
#
# Inputs:
#   - output/results/cv/fold_*_<param_tag>.rds  (10 fold results)
#   - output/results/exp6_species_cv_<param_tag>.csv
#   - output/results/exp6_sites_cv_<param_tag>.csv
#   - output/data_calanda_jsdm_*.rds  (for Y matrix and site info)
#
# Outputs:
#   - plot/venn_cv.pdf
#   - plot/ternary_species_cv.pdf
#   - plot/ternary_sites_cv.pdf
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

# --- Colors ---
color_env    = "#81caf3"
color_spa    = "#d00000"
color_codist = "#00bd89"

# ==============================================================================
# LOAD DATA
# ==============================================================================
cat("\n=== Variance partitioning (10-fold CV) ===\n")

# Detect param_tag from fold files
fold_files = sort(list.files(
  here("Calanda_JSDM", "output", "results", "cv"),
  pattern = "^fold_.*\\.rds$", full.names = TRUE
))
param_tag = regmatches(basename(fold_files[1]),
                       regexpr("a[0-9.]+_l[0-9.]+_le[0-9.]+(?=\\.)",
                               basename(fold_files[1]), perl = TRUE))
k = length(fold_files)
fold_data_list = lapply(fold_files, readRDS)

cat(sprintf("  Loaded %d folds | param_tag = %s\n", k, param_tag))

# Species and site CV metrics
species_cv = read_csv(here("Calanda_JSDM", "output", "results",
                           paste0("exp6_species_cv_", param_tag, ".csv")),
                      show_col_types = FALSE)
site_cv = read_csv(here("Calanda_JSDM", "output", "results",
                        paste0("exp6_sites_cv_", param_tag, ".csv")),
                   show_col_types = FALSE)

# Load Y matrix and environmental data
data_files = sort(list.files(
  here("Calanda_JSDM", "output"),
  pattern = "^data_calanda_jsdm_[0-9].*\\.rds$", full.names = TRUE
), decreasing = TRUE)
data_calanda = readRDS(data_files[1])
Y = data_calanda$Y
sp_names = colnames(Y)
n_species = ncol(Y)
n_sites = nrow(Y)

# ==============================================================================
# 1. MODEL-LEVEL R² — McFadden (anova) and Tjur (discrimination)
# ==============================================================================
cat("\n--- Model-level R² across folds ---\n")

# 1a) McFadden R² from anova (7 components + Full)
anova_model_names = c("F_A", "F_B", "F_AB", "F_S", "F_AS", "F_BS", "F_ABS", "Full")
anova_labels = c("Environment", "Biotic", "Env x Bio", "Spatial",
                 "Env x Spa", "Bio x Spa", "Shared", "Full")

anova_r2_list = list()
tjur_r2_vec = numeric(k)

for (fi in seq_len(k)) {
  fd = fold_data_list[[fi]]
  anova_res = fd$partition$anova$results
  r2_mcf = setNames(anova_res$`R2 McFadden`, anova_res$models)
  anova_r2_list[[fi]] = tibble(
    fold = fi,
    model = anova_model_names,
    label = anova_labels,
    mcfadden_r2 = as.numeric(r2_mcf[anova_model_names])
  )
  tjur_r2_vec[fi] = fd$partition$R2
}

df_anova_r2 = bind_rows(anova_r2_list)

# Summary: median + 95% CI for each anova component
df_anova_summary = df_anova_r2 %>%
  group_by(model, label) %>%
  summarise(
    median_r2 = median(mcfadden_r2, na.rm = TRUE),
    mean_r2   = mean(mcfadden_r2, na.rm = TRUE),
    sd_r2     = sd(mcfadden_r2, na.rm = TRUE),
    ci_lo     = mean_r2 - qt(0.975, n() - 1) * sd_r2 / sqrt(n()),
    ci_hi     = mean_r2 + qt(0.975, n() - 1) * sd_r2 / sqrt(n()),
    .groups   = "drop"
  )

# Tjur R² summary
tjur_summary = tibble(
  metric = "Tjur R2 (discrimination)",
  mean_val = mean(tjur_r2_vec),
  median_val = median(tjur_r2_vec),
  sd_val = sd(tjur_r2_vec),
  ci_lo = mean_val - qt(0.975, k - 1) * sd_val / sqrt(k),
  ci_hi = mean_val + qt(0.975, k - 1) * sd_val / sqrt(k)
)

cat("\n  McFadden R² (anova decomposition) — median [mean +/- 95% CI]:\n")
for (i in seq_len(nrow(df_anova_summary))) {
  row = df_anova_summary[i, ]
  cat(sprintf("    %-15s: median = %.4f  [%.4f +/- %.4f, %.4f]\n",
              row$label, row$median_r2, row$mean_r2, row$ci_lo, row$ci_hi))
}

cat(sprintf("\n  Tjur R² (discrimination): mean = %.4f, median = %.4f [%.4f, %.4f]\n",
            tjur_summary$mean_val, tjur_summary$median_val,
            tjur_summary$ci_lo, tjur_summary$ci_hi))

# ==============================================================================
# 2. VENN DIAGRAM — median McFadden R² + CI table
# ==============================================================================
cat("\n--- Venn diagram ---\n")

# Build a fake anova object with median R² to pass to plot_anova_custom
median_r2 = df_anova_summary %>% select(model, median_r2)
venn_results = data.frame(
  models = c(anova_model_names, "Saturated", "Null"),
  `ll` = NA,
  `Residual deviance` = NA,
  `Deviance` = NA,
  `R2 Nagelkerke` = NA,
  `R2 McFadden` = c(median_r2$median_r2, NA, 0),
  check.names = FALSE
)
venn_obj = list(results = venn_results, spatial = TRUE)

# CI table for annotation
ci_table = df_anova_summary %>%
  filter(model != "Full") %>%
  mutate(ci_text = sprintf("%s: %.3f [%.3f, %.3f]",
                           label, median_r2, ci_lo, ci_hi)) %>%
  pull(ci_text)
full_row = df_anova_summary %>% filter(model == "Full")
full_text = sprintf("Full model: %.3f [%.3f, %.3f]",
                    full_row$median_r2, full_row$ci_lo, full_row$ci_hi)
tjur_text = sprintf("Tjur R2: %.3f [%.3f, %.3f]",
                    tjur_summary$median_val, tjur_summary$ci_lo, tjur_summary$ci_hi)

pdf(here("Calanda_JSDM", "plot", "venn_cv.pdf"), width = 12, height = 8)
layout(matrix(c(1, 2), nrow = 1), widths = c(2, 1.2))
# Left: Venn
par(mar = c(2, 2, 3, 0))
plot_anova_custom(venn_obj, cols = c(color_env, color_codist, color_spa))
title(main = sprintf("Variance Partitioning (median across %d folds)", k),
      cex.main = 1.2)
# Right: CI table
par(mar = c(2, 0, 3, 2))
plot(NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
title(main = "McFadden R\u00B2: median [95% CI]", cex.main = 1)
all_lines = c(full_text, "", ci_table, "", tjur_text)
n_lines = length(all_lines)
for (li in seq_along(all_lines)) {
  y_pos = 0.95 - (li - 1) * 0.08
  font_face = if (li == 1 || li == n_lines) 2 else 1  # bold for full + tjur
  text(0.05, y_pos, all_lines[li], adj = c(0, 1), cex = 0.85, font = font_face)
}
dev.off()
cat("  Saved venn_cv.pdf\n")

# ==============================================================================
# 3. SPECIES-LEVEL VP — aggregate across folds, filter by AUC
# ==============================================================================
cat("\n--- Species VP (filtered by AUC >= 0.7) ---\n")

sp_vp_list = list()
for (fi in seq_len(k)) {
  fd = fold_data_list[[fi]]
  sp = fd$partition$species
  sp_vp_list[[fi]] = tibble(
    fold = fi,
    species_idx = seq_len(nrow(sp)),
    r2 = sp[, "r2"],
    env = sp[, "env"],
    spa = sp[, "spa"],
    codist = sp[, "codist"]
  )
}
df_sp_vp = bind_rows(sp_vp_list)

# Mean + CI per species per component
df_sp_summary = df_sp_vp %>%
  pivot_longer(cols = c(r2, env, spa, codist), names_to = "component", values_to = "value") %>%
  group_by(species_idx, component) %>%
  summarise(
    mean_val = mean(value, na.rm = TRUE),
    sd_val   = sd(value, na.rm = TRUE),
    ci_lo    = mean_val - qt(0.975, n() - 1) * sd_val / sqrt(n()),
    ci_hi    = mean_val + qt(0.975, n() - 1) * sd_val / sqrt(n()),
    ci_width = ci_hi - ci_lo,
    .groups  = "drop"
  ) %>%
  mutate(
    species_name = sp_names[species_idx],
    weight = 1 / pmax(ci_width, 0.001)^2
  )

# Join AUC
df_sp_summary = df_sp_summary %>%
  left_join(tibble(species_idx = seq_len(n_species), auc = species_cv$test_auc),
            by = "species_idx")

# Filter: AUC >= 0.7
sp_keep = species_cv$test_auc >= 0.7 & !is.na(species_cv$test_auc)
n_sp_keep = sum(sp_keep)
n_sp_drop = n_species - n_sp_keep
cat(sprintf("  %d / %d species kept (AUC >= 0.7), %d removed\n",
            n_sp_keep, n_species, n_sp_drop))

df_sp_filtered = df_sp_summary %>% filter(sp_keep[species_idx])

# Wide format for ternary
df_sp_tern = df_sp_filtered %>%
  filter(component %in% c("env", "spa", "codist")) %>%
  select(species_idx, species_name, component, mean_val) %>%
  pivot_wider(names_from = component, values_from = mean_val)

# ==============================================================================
# 4. SITE-LEVEL VP — aggregate across folds, filter by logloss
# ==============================================================================
cat("\n--- Site VP (filtered by logloss <= log(2)) ---\n")

si_vp_list = list()
for (fi in seq_len(k)) {
  fd = fold_data_list[[fi]]
  si = fd$partition$sites
  si_vp_list[[fi]] = tibble(
    fold = fi,
    site_idx = fd$train_idx,
    r2 = si[, "r2"],
    env = si[, "env"],
    spa = si[, "spa"],
    codist = si[, "codist"]
  )
}
df_si_vp = bind_rows(si_vp_list)

df_si_summary = df_si_vp %>%
  pivot_longer(cols = c(r2, env, spa, codist), names_to = "component", values_to = "value") %>%
  group_by(site_idx, component) %>%
  summarise(
    mean_val = mean(value, na.rm = TRUE),
    sd_val   = sd(value, na.rm = TRUE),
    ci_lo    = mean_val - qt(0.975, n() - 1) * sd_val / sqrt(n()),
    ci_hi    = mean_val + qt(0.975, n() - 1) * sd_val / sqrt(n()),
    ci_width = ci_hi - ci_lo,
    .groups  = "drop"
  ) %>%
  mutate(weight = 1 / pmax(ci_width, 0.001)^2)

# Join logloss
df_si_summary = df_si_summary %>%
  left_join(tibble(site_idx = seq_len(n_sites), logloss = site_cv$logloss),
            by = "site_idx")

# Filter: logloss <= log(2)
si_keep = site_cv$logloss <= log(2) & !is.na(site_cv$logloss)
n_si_keep = sum(si_keep)
n_si_drop = n_sites - n_si_keep
cat(sprintf("  %d / %d sites kept (logloss <= log(2) = %.3f), %d removed\n",
            n_si_keep, n_sites, log(2), n_si_drop))

df_si_filtered = df_si_summary %>% filter(si_keep[site_idx])

# Wide format for ternary
df_si_tern = df_si_filtered %>%
  filter(component %in% c("env", "spa", "codist")) %>%
  select(site_idx, component, mean_val) %>%
  pivot_wider(names_from = component, values_from = mean_val)

# ==============================================================================
# 5. PREPARE WIDE VP DATA FOR SCATTER PLOTS
# ==============================================================================
cat("\n--- Preparing VP scatter data ---\n")

# Species: add r2 to wide format
df_sp_r2 = df_sp_filtered %>%
  filter(component == "r2") %>%
  select(species_idx, r2 = mean_val)
df_sp_tern = df_sp_tern %>%
  left_join(df_sp_r2, by = "species_idx")

# Sites: add r2 to wide format
df_si_r2 = df_si_filtered %>%
  filter(component == "r2") %>%
  select(site_idx, r2 = mean_val)
df_si_tern = df_si_tern %>%
  left_join(df_si_r2, by = "site_idx")

# ==============================================================================
# 6. SAVE SUMMARY DATAFRAMES
# ==============================================================================
cat("\n--- Saving summary dataframes ---\n")

write_csv(df_sp_summary,
          here("Calanda_JSDM", "output", "results",
               paste0("vp_species_summary_", param_tag, ".csv")))
write_csv(df_si_summary,
          here("Calanda_JSDM", "output", "results",
               paste0("vp_sites_summary_", param_tag, ".csv")))
write_csv(df_anova_summary,
          here("Calanda_JSDM", "output", "results",
               paste0("vp_anova_summary_", param_tag, ".csv")))

cat(sprintf("  Saved vp_species_summary_%s.csv (%d species x %d components)\n",
            param_tag, n_species, 4))
cat(sprintf("  Saved vp_sites_summary_%s.csv (%d sites x %d components)\n",
            param_tag, n_sites, 4))
cat(sprintf("  Saved vp_anova_summary_%s.csv\n", param_tag))

# ==============================================================================
# 7. CONSOLIDATED LANDSCAPE PDF
# ==============================================================================
cat("\n--- Building consolidated landscape PDF ---\n")

# --- Row 1b: VP component violin plots (species + sites side-by-side) ---
# Prepare long data for violin: species
df_sp_violin = df_sp_filtered %>%
  filter(component %in% c("env", "spa", "codist", "r2")) %>%
  mutate(
    type = "Species",
    component = factor(component,
                       levels = c("env", "codist", "spa", "r2"),
                       labels = c("Environment", "Biotic", "Spatial", "R^2"))
  )

# Sites
df_si_violin = df_si_filtered %>%
  filter(component %in% c("env", "spa", "codist", "r2")) %>%
  mutate(
    type = "Sites",
    component = factor(component,
                       levels = c("env", "codist", "spa", "r2"),
                       labels = c("Environment", "Biotic", "Spatial", "R^2"))
  )

df_violin = bind_rows(df_sp_violin, df_si_violin) %>%
  mutate(type = factor(type, levels = c("Species", "Sites")))

# Compute medians for annotation
df_medians = df_violin %>%
  group_by(component, type) %>%
  summarise(median_val = median(mean_val, na.rm = TRUE), .groups = "drop")

# Component violin colors
comp_colors = c("Environment" = color_env, "Biotic" = color_codist,
                "Spatial" = color_spa, "R^2" = "grey50")

# Violin fill colors per component
violin_fill = c("Environment" = color_env, "Biotic" = color_codist,
                "Spatial" = color_spa, "R^2" = "grey60")
df_violin = df_violin %>%
  mutate(comp_fill = violin_fill[as.character(component)])

p_violins = ggplot(df_violin, aes(x = type, y = mean_val)) +
  geom_violin(aes(fill = component), alpha = 0.2, color = NA) +
  geom_jitter(aes(color = ci_width), width = 0.15, size = 1, alpha = 0.5) +
  geom_crossbar(data = df_medians,
                aes(x = type, y = median_val, ymin = median_val, ymax = median_val),
                width = 0.5, linewidth = 0.6, color = "black") +
  geom_text(data = df_medians,
            aes(x = type, y = median_val,
                label = sprintf("%.3f", median_val)),
            vjust = -0.8, hjust = -0.9, size = 3.5, fontface = "bold") +
  scale_color_gradient(low = "grey10", high = "grey90", name = "CI width") +
  scale_fill_manual(values = violin_fill, guide = "none") +
  facet_wrap(~ component, nrow = 1, scales = "free_y",
             labeller = label_parsed) +
  labs(x = NULL, y = "Variance proportion") +
  theme_bw(base_size = 18) +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 18),
    strip.text = element_text(size = 17, face = "bold"),
    legend.position = "right"
  )

# --- Row 2: Scatter plots with correlations ---
# Helper to add correlation annotations
build_vp_scatter_combined = function(df, label_species = FALSE) {
  # Pearson correlations
  cor_eb = cor.test(df$env, df$codist, method = "pearson")
  cor_es = cor.test(df$env, df$spa, method = "pearson")

  cor_label = paste0(
    "r(Env, Bio) = ", sprintf("%.2f", cor_eb$estimate),
    "\nr(Env, Spa) = ", sprintf("%.2f", cor_es$estimate)
  )

  p = ggplot(df, aes(x = env, y = codist, color = spa, size = r2)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey60") +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.3, linetype = "dashed", color = "grey50") +
    geom_point(alpha = 0.7) +
    scale_color_gradient(low = "grey85", high = color_spa, name = "Spatial") +
    scale_size_continuous(range = c(1, 6), name = expression(R^2)) +
    annotate("text",
             x = 0.02,
             y = max(df$codist, na.rm = TRUE),
             label = cor_label, hjust = 0, vjust = 1, size = 4.5,
             fontface = "italic") +
    labs(x = "Environment", y = "Biotic") +
    theme_bw(base_size = 18) +
    theme(
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      legend.position = "right"
    )

  if (label_species && "species_name" %in% names(df)) {
    labeled = bind_rows(
      df %>% slice_max(env, n = 3),
      df %>% slice_max(spa, n = 3),
      df %>% slice_max(codist, n = 3),
      df %>% slice_min(r2, n = 2),
      df %>% slice_max(r2, n = 2)
    ) %>% distinct(species_name, .keep_all = TRUE)
    p = p + geom_text_repel(data = labeled, aes(label = species_name),
                            size = 3.5, fontface = "italic",
                            max.overlaps = 15, segment.color = "grey50",
                            show.legend = FALSE)
  }
  p
}

p_scatter_sp = build_vp_scatter_combined(df_sp_tern, label_species = TRUE)
p_scatter_si = build_vp_scatter_combined(df_si_tern)

# --- Compose full figure ---
# Row 1: Venn (base R) captured as grob + violins
venn_grob = wrap_elements(full = ~ {
  par(mar = c(1, 1, 2, 1))
  plot_anova_custom(venn_obj, cols = c(color_env, color_codist, color_spa))
  # Add CI annotation below
  full_r2 = df_anova_summary %>% filter(model == "Full")
  mtext(sprintf("Full McFadden R\u00B2: %.3f [%.3f, %.3f]",
                full_r2$median_r2, full_r2$ci_lo, full_r2$ci_hi),
        side = 1, line = 0, cex = 0.9, font = 2)
  mtext(sprintf("Tjur R\u00B2: %.3f [%.3f, %.3f]",
                tjur_summary$median_val, tjur_summary$ci_lo, tjur_summary$ci_hi),
        side = 1, line = 1, cex = 0.85)
})

# Set row 1 widths: Venn narrower, violins wider
row1 = venn_grob + p_violins + plot_layout(widths = c(1, 2.5))
row2 = p_scatter_sp + p_scatter_si + plot_layout(widths = c(1, 1))
layout_combined = row1 / row2 + plot_layout(heights = c(0.8, 1))

pdf(here("Calanda_JSDM", "plot", "vp_combined.pdf"), width = 22, height = 14)
print(layout_combined)
dev.off()
cat("  Saved vp_combined.pdf\n")

cat("\n=== Variance partitioning complete ===\n")
