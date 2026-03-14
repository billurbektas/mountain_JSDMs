# ==============================================================================
# Script: 12_post_exp6_cv.R
# Purpose: Post-analysis of sjSDM Experiment 6 - k-fold cross-validation
#          1) Predictive performance: species AUC violin + site log-loss violin
#          2) VP stability across folds: full model R2, species-level and
#             site-level E/S/C/R2 means with 95% CIs across 10 folds
#
# Inputs:
#   - output/results/exp6_species_cv_<param_tag>.csv
#   - output/results/exp6_sites_cv_<param_tag>.csv
#   - output/results/cv/fold_*_<param_tag>.rds (predictions + partition)
#   - data_calanda_jsdm_*.rds (for Y matrix)
#   - param_tag format: a<alpha>_l<lambda>_le<lambda_env>
#
# Outputs:
#   - plot/exp6_prediction_violins_<param_tag>.pdf
#   - plot/exp6_vp_model_stability_<param_tag>.pdf
#   - plot/exp6_vp_species_stability_<param_tag>.pdf
#   - plot/exp6_vp_sites_stability_<param_tag>.pdf
#   - output/results/exp6_cv_summary_<param_tag>.csv
#
# Requires: R/00_setup/functions_calanda.R
# ==============================================================================

library(tidyverse)
library(conflicted)
library(here)
library(patchwork)
library(pROC)

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
cat("\n=== Step 1: Loading exp6 CV results ===\n")

# Detect param_tag from species CV file
species_cv_files = list.files(
  here("Calanda_JSDM", "output", "results"),
  pattern = "^exp6_species_cv_.*\\.csv$", full.names = TRUE
)

tag_match = regmatches(basename(species_cv_files[1]),
                       regexpr("a[0-9.]+_l[0-9.]+_le[0-9.]+(?=\\.)", basename(species_cv_files[1]), perl = TRUE))
if (length(tag_match) == 0 || tag_match == "") {
  stop("Could not detect param_tag from filenames.")
}
param_tag = tag_match
cat(sprintf("  Detected param_tag: %s\n", param_tag))

# Load species and site CSVs
species_cv = read_csv(here("Calanda_JSDM", "output", "results",
                           paste0("exp6_species_cv_", param_tag, ".csv")),
                      show_col_types = FALSE)
site_cv = read_csv(here("Calanda_JSDM", "output", "results",
                        paste0("exp6_sites_cv_", param_tag, ".csv")),
                   show_col_types = FALSE)

cat(sprintf("  Species: %d (test AUC available: %d)\n",
            nrow(species_cv), sum(!is.na(species_cv$test_auc))))
cat(sprintf("  Sites: %d (log-loss available: %d)\n",
            nrow(site_cv), sum(!is.na(site_cv$logloss))))

# Load per-fold RDS files
fold_files = list.files(
  here("Calanda_JSDM", "output", "results", "cv"),
  pattern = paste0("^fold_.*_", param_tag, "\\.rds$"), full.names = TRUE
)
fold_files = sort(fold_files)
k = length(fold_files)
cat(sprintf("  Folds: %d\n", k))

fold_data_list = lapply(fold_files, readRDS)

# Extract hyperparameters
hp_alpha      = fold_data_list[[1]]$alpha
hp_lambda     = fold_data_list[[1]]$lambda
hp_lambda_env = fold_data_list[[1]]$lambda_env
hp_fit_sampling = fold_data_list[[1]]$fit_sampling
hp_subtitle = sprintf("alpha = %s, lambda = %s, lambda_env = %s | fit sampling = %s | %d-fold CV",
                       hp_alpha, hp_lambda, hp_lambda_env, hp_fit_sampling, k)
cat(sprintf("  Hyperparameters: %s\n", hp_subtitle))

# Load Y matrix
data_files = sort(list.files(
  here("Calanda_JSDM", "output"),
  pattern = "^data_calanda_jsdm_[0-9].*\\.rds$", full.names = TRUE
), decreasing = TRUE)
data_calanda = readRDS(data_files[1])
cat(sprintf("  Data file: %s\n", basename(data_files[1])))
Y = data_calanda$Y
n_species = ncol(Y)
n_sites   = nrow(Y)
cat(sprintf("  Y matrix: %d sites x %d species\n", n_sites, n_species))

# ==============================================================================
# STEP 2: PREDICTION VIOLINS (Species AUC + Site log-loss)
# ==============================================================================
cat("\n=== Step 2: Prediction violin plots ===\n")

df_auc = species_cv %>%
  filter(!is.na(test_auc)) %>%
  transmute(metric = "Species test AUC", value = test_auc)

df_ll = site_cv %>%
  filter(!is.na(logloss)) %>%
  transmute(metric = "Site log-loss", value = logloss)

# Null model log-loss: predict training prevalence per fold (fair OOS comparison)
eps = 1e-7
folds = readRDS(here("Calanda_JSDM", "output", "results", "cv", "stratified_folds.rds"))
null_logloss = numeric(nrow(Y))
for (fold_i in seq_len(k)) {
  train_idx = setdiff(seq_len(nrow(Y)), folds[[fold_i]])
  test_idx = folds[[fold_i]]
  prevalence_train = colMeans(Y[train_idx, ])
  p_null = pmin(pmax(prevalence_train, eps), 1 - eps)
  for (i in test_idx) {
    null_logloss[i] = -mean(Y[i, ] * log(p_null) + (1 - Y[i, ]) * log(1 - p_null))
  }
}
null_logloss_mean = mean(null_logloss)
cat(sprintf("  Null model (prevalence) site log-loss: mean = %.3f\n", null_logloss_mean))

# Species AUC violin
p_auc = ggplot(df_auc, aes(x = metric, y = value)) +
  geom_violin(fill = color_env, alpha = 0.4, draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 0.8, color = color_env) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey40") +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "grey30") +
  annotate("text", x = 0.55, y = 0.51, label = "random = 0.5", size = 3, color = "grey40") +
  annotate("text", x = 0.55, y = 0.71,
           label = sprintf("0.7  (%d/%d above)", sum(df_auc$value > 0.7), nrow(df_auc)),
           size = 3, color = "grey30") +
  labs(title = "A) Species test AUC", x = NULL, y = "AUC") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Site log-loss violin
p_ll = ggplot(df_ll, aes(x = metric, y = value)) +
  geom_violin(fill = color_spa, alpha = 0.4, draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter(width = 0.15, alpha = 0.15, size = 0.5, color = color_spa) +
  geom_hline(yintercept = log(2), linetype = "dotted", color = "grey40") +
  annotate("text", x = 0.55, y = log(2) + 0.01,
           label = sprintf("random = %.3f  (%d/%d below)", log(2),
                           sum(df_ll$value < log(2)), nrow(df_ll)),
           size = 3, color = "grey40") +
  geom_hline(yintercept = null_logloss_mean, linetype = "dashed", color = "grey30") +
  annotate("text", x = 0.55, y = null_logloss_mean + 0.01,
           label = sprintf("prevalence = %.3f  (%d/%d below)", null_logloss_mean,
                           sum(df_ll$value < null_logloss_mean), nrow(df_ll)),
           size = 3, color = "grey30") +
  labs(title = "B) Site log-loss", x = NULL, y = "Log-loss (lower = better)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

p_violins = (p_auc | p_ll) +
  plot_annotation(
    title = "Experiment 6: Predictive performance (out-of-fold)",
    subtitle = hp_subtitle,
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10)
    )
  )

pdf(here("Calanda_JSDM", "plot", "detailed_modeling_experiment_figures", paste0("exp6_prediction_violins_", param_tag, ".pdf")),
    width = 10, height = 7)
print(p_violins)
dev.off()
cat(sprintf("  Saved exp6_prediction_violins_%s.pdf\n", param_tag))

# ==============================================================================
# STEP 3: VP STABILITY — FULL MODEL R2 ACROSS FOLDS
# ==============================================================================
cat("\n=== Step 3: VP model-level stability across folds ===\n")

# Extract model-level anova R2 per fold
model_vp_list = list()

for (fi in seq_along(fold_data_list)) {
  fd = fold_data_list[[fi]]
  p  = fd$partition

  anova_res = p$anova$results
  r2_mcf = setNames(anova_res$`R2 McFadden`, anova_res$models)

  model_vp_list[[fi]] = tibble(
    fold = fd$fold,
    Overall     = r2_mcf["Full"],
    Environment = r2_mcf["F_A"],
    Spatial     = r2_mcf["F_S"],
    Biotic      = r2_mcf["F_B"]
  )
}

df_model_vp = bind_rows(model_vp_list) %>%
  pivot_longer(-fold, names_to = "component", values_to = "R2") %>%
  mutate(component = factor(component, levels = names(unit_colors)))

# Summary stats
df_model_summary = df_model_vp %>%
  group_by(component) %>%
  summarise(
    mean_R2 = mean(R2, na.rm = TRUE),
    sd_R2   = sd(R2, na.rm = TRUE),
    ci_lo   = mean_R2 - qt(0.975, k - 1) * sd_R2 / sqrt(k),
    ci_hi   = mean_R2 + qt(0.975, k - 1) * sd_R2 / sqrt(k),
    .groups = "drop"
  )

p_model_vp = ggplot(df_model_vp, aes(x = component, y = R2, color = component)) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
  geom_crossbar(data = df_model_summary,
                aes(x = component, y = mean_R2, ymin = ci_lo, ymax = ci_hi, fill = component),
                alpha = 0.2, width = 0.4, linewidth = 0.4) +
  scale_color_manual(values = unit_colors) +
  scale_fill_manual(values = unit_colors) +
  geom_text(data = df_model_summary,
            aes(x = component, y = ci_hi,
                label = sprintf("%.4f\n[%.4f]", mean_R2, sd_R2)),
            vjust = -0.3, size = 3, show.legend = FALSE) +
  labs(
    title = "Experiment 6: Model-level McFadden R2 across folds",
    subtitle = paste0(hp_subtitle, "\nDots = per-fold R2, bars = mean +/- 95% CI"),
    x = NULL, y = "McFadden R2"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none") +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(15, 5.5, 5.5, 5.5))

pdf(here("Calanda_JSDM", "plot", "detailed_modeling_experiment_figures", paste0("exp6_vp_model_stability_", param_tag, ".pdf")),
    width = 8, height = 7)
print(p_model_vp)
dev.off()
cat(sprintf("  Saved exp6_vp_model_stability_%s.pdf\n", param_tag))

# ==============================================================================
# STEP 4: VP STABILITY — SPECIES-LEVEL MEANS ACROSS FOLDS
# ==============================================================================
cat("\n=== Step 4: VP species-level stability across folds ===\n")

# For each fold, compute mean species-level E/S/C/R2
sp_vp_list = list()

for (fi in seq_along(fold_data_list)) {
  fd = fold_data_list[[fi]]
  sp = fd$partition$species

  sp_vp_list[[fi]] = tibble(
    fold = fd$fold,
    Overall     = mean(sp$r2, na.rm = TRUE),
    Environment = mean(sp$env, na.rm = TRUE),
    Spatial     = mean(sp$spa, na.rm = TRUE),
    Biotic      = mean(sp$codist, na.rm = TRUE)
  )
}

df_sp_vp = bind_rows(sp_vp_list) %>%
  pivot_longer(-fold, names_to = "component", values_to = "mean_value") %>%
  mutate(component = factor(component, levels = names(unit_colors)))

df_sp_summary = df_sp_vp %>%
  group_by(component) %>%
  summarise(
    grand_mean = mean(mean_value, na.rm = TRUE),
    sd_val     = sd(mean_value, na.rm = TRUE),
    ci_lo      = grand_mean - qt(0.975, k - 1) * sd_val / sqrt(k),
    ci_hi      = grand_mean + qt(0.975, k - 1) * sd_val / sqrt(k),
    .groups = "drop"
  )

# A) Crossbar plot of means across folds
p_sp_means = ggplot(df_sp_vp, aes(x = component, y = mean_value, color = component)) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
  geom_crossbar(data = df_sp_summary,
                aes(x = component, y = grand_mean, ymin = ci_lo, ymax = ci_hi, fill = component),
                alpha = 0.2, width = 0.4, linewidth = 0.4) +
  scale_color_manual(values = unit_colors) +
  scale_fill_manual(values = unit_colors) +
  geom_text(data = df_sp_summary,
            aes(x = component, y = ci_hi,
                label = sprintf("%.4f\n[%.4f]", grand_mean, sd_val)),
            vjust = -0.3, size = 3, show.legend = FALSE) +
  labs(
    title = "A) Mean species VP across folds",
    subtitle = "Dots = per-fold mean, bars = grand mean +/- 95% CI",
    x = NULL, y = "Mean species-level value"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none") +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(15, 5.5, 5.5, 5.5))

# B) Per-species SD across folds (how variable is each species?)
sp_per_species_list = list()
for (fi in seq_along(fold_data_list)) {
  fd = fold_data_list[[fi]]
  sp = fd$partition$species
  for (col in c("r2", "env", "spa", "codist")) {
    comp_label = recode(col, r2 = "Overall", env = "Environment",
                        spa = "Spatial", codist = "Biotic")
    sp_per_species_list = c(sp_per_species_list, list(tibble(
      fold = fd$fold, component = comp_label,
      species_idx = seq_len(nrow(sp)),
      value = sp[, col]
    )))
  }
}

df_sp_per_species = bind_rows(sp_per_species_list) %>%
  mutate(component = factor(component, levels = names(unit_colors)))

df_sp_species_sd = df_sp_per_species %>%
  group_by(component, species_idx) %>%
  summarise(sd_across_folds = sd(value, na.rm = TRUE), .groups = "drop")

p_sp_sd = ggplot(df_sp_species_sd, aes(x = sd_across_folds, fill = component)) +
  geom_histogram(bins = 30, alpha = 0.7, color = "white") +
  facet_wrap(~ component, scales = "free_y") +
  scale_fill_manual(values = unit_colors) +
  labs(
    title = "B) Per-species SD of VP values across folds",
    subtitle = "How much does each species' partition change across 10 training sets?",
    x = "SD across folds", y = "Count"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

p_sp_stability = p_sp_means / p_sp_sd +
  plot_annotation(
    title = "Experiment 6: Species-level VP stability across folds",
    subtitle = hp_subtitle,
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10)
    )
  )

pdf(here("Calanda_JSDM", "plot", "detailed_modeling_experiment_figures", paste0("exp6_vp_species_stability_", param_tag, ".pdf")),
    width = 12, height = 12)
print(p_sp_stability)
dev.off()
cat(sprintf("  Saved exp6_vp_species_stability_%s.pdf\n", param_tag))

# --- C) Per-species ranked mean + CI across folds ---
cat("  Building per-species VP error bar plot...\n")

# Get species names from Y
sp_names = colnames(Y)

# Count how many folds each species has presences in the training set
n_folds_with_presence = matrix(NA, nrow = k, ncol = n_species)
for (fi in seq_along(fold_data_list)) {
  fd = fold_data_list[[fi]]
  n_folds_with_presence[fi, ] = colSums(Y[fd$train_idx, ]) > 0
}
sp_fold_count = colSums(n_folds_with_presence)

df_sp_species_ci = df_sp_per_species %>%
  group_by(component, species_idx) %>%
  summarise(
    mean_val = mean(value, na.rm = TRUE),
    sd_val   = sd(value, na.rm = TRUE),
    ci_lo    = mean_val - qt(0.975, n() - 1) * sd_val / sqrt(n()),
    ci_hi    = mean_val + qt(0.975, n() - 1) * sd_val / sqrt(n()),
    .groups  = "drop"
  ) %>%
  mutate(
    species_name = sp_names[species_idx],
    n_folds_present = sp_fold_count[species_idx]
  )

cat(sprintf("  Species in all %d folds: %d / %d\n",
            k, sum(sp_fold_count == k), n_species))
cat(sprintf("  Species missing from >= 1 fold: %d\n",
            sum(sp_fold_count < k)))

# Rank species by mean value within each component
df_sp_species_ci = df_sp_species_ci %>%
  group_by(component) %>%
  mutate(rank = rank(mean_val, ties.method = "first")) %>%
  ungroup()

p_sp_ranked = ggplot(df_sp_species_ci,
                     aes(x = rank, y = mean_val, color = component)) +
  geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi, alpha = n_folds_present / k),
                  size = 0.15, linewidth = 0.2) +
  facet_wrap(~ component, scales = "free_y") +
  scale_color_manual(values = unit_colors) +
  scale_alpha_continuous(range = c(0.3, 1), name = paste0("Folds with\npresence (/ ", k, ")")) +
  labs(
    title = "Experiment 6: Per-species VP values (mean +/- 95% CI across folds)",
    subtitle = paste0("Species ranked by mean value | opacity = fraction of folds with presences\n", hp_subtitle),
    x = "Species (ranked)", y = "VP value"
  ) +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(face = "bold"))

pdf(here("Calanda_JSDM", "plot", "detailed_modeling_experiment_figures", paste0("exp6_vp_species_ranked_", param_tag, ".pdf")),
    width = 14, height = 10)
print(p_sp_ranked)
dev.off()
cat(sprintf("  Saved exp6_vp_species_ranked_%s.pdf\n", param_tag))

# --- D) Per-species ranked PROPORTIONS (mean + CI across folds) ---
cat("  Building per-species VP proportion plots...\n")

# Compute proportions per species per fold: E/(E+S+C), S/(E+S+C), C/(E+S+C)
sp_prop_list = list()
sp_sum_vs_r2 = list()
for (fi in seq_along(fold_data_list)) {
  fd = fold_data_list[[fi]]
  sp = fd$partition$species
  total = sp$env + sp$spa + sp$codist
  sp_sum_vs_r2[[fi]] = tibble(
    fold = fd$fold,
    species_idx = seq_len(nrow(sp)),
    r2 = sp$r2,
    sum_esc = total
  )
  for (comp in c("env", "spa", "codist")) {
    comp_label = recode(comp, env = "Environment", spa = "Spatial", codist = "Biotic")
    sp_prop_list = c(sp_prop_list, list(tibble(
      fold = fd$fold, component = comp_label,
      species_idx = seq_len(nrow(sp)),
      proportion = sp[[comp]] / total
    )))
  }
}

# Sanity check: sum(E+S+C) vs R2
df_sp_sum_check = bind_rows(sp_sum_vs_r2)
sp_sum_r2_cor = cor(df_sp_sum_check$sum_esc, df_sp_sum_check$r2, use = "complete.obs")
sp_sum_r2_diff = mean(abs(df_sp_sum_check$sum_esc - df_sp_sum_check$r2), na.rm = TRUE)
cat(sprintf("  Sanity check (species): cor(sum_ESC, R2) = %.4f, mean|diff| = %.6f\n",
            sp_sum_r2_cor, sp_sum_r2_diff))

prop_colors = c("Environment" = color_env, "Spatial" = color_spa, "Biotic" = color_codist)

df_sp_prop = bind_rows(sp_prop_list) %>%
  mutate(component = factor(component, levels = c("Environment", "Spatial", "Biotic")))

df_sp_prop_ci = df_sp_prop %>%
  group_by(component, species_idx) %>%
  summarise(
    mean_val = mean(proportion, na.rm = TRUE),
    sd_val   = sd(proportion, na.rm = TRUE),
    ci_lo    = mean_val - qt(0.975, n() - 1) * sd_val / sqrt(n()),
    ci_hi    = mean_val + qt(0.975, n() - 1) * sd_val / sqrt(n()),
    .groups  = "drop"
  ) %>%
  mutate(species_name = sp_names[species_idx]) %>%
  group_by(component) %>%
  mutate(rank = rank(mean_val, ties.method = "first")) %>%
  ungroup()

p_sp_prop_ranked = ggplot(df_sp_prop_ci,
                          aes(x = rank, y = mean_val, color = component)) +
  geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi),
                  size = 0.15, linewidth = 0.2, alpha = 0.7) +
  facet_wrap(~ component, scales = "free_y") +
  scale_color_manual(values = prop_colors) +
  labs(
    title = "Experiment 6: Per-species VP proportions (mean +/- 95% CI across folds)",
    subtitle = paste0("Proportion = component / (E+S+C) | Species ranked by mean proportion\n",
                      sprintf("Sanity check: cor(E+S+C, R2) = %.4f, mean|diff| = %.6f\n", sp_sum_r2_cor, sp_sum_r2_diff),
                      hp_subtitle),
    x = "Species (ranked)", y = "Proportion of explained variance"
  ) +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(face = "bold"))

pdf(here("Calanda_JSDM", "plot", "detailed_modeling_experiment_figures", paste0("exp6_vp_species_proportions_", param_tag, ".pdf")),
    width = 14, height = 10)
print(p_sp_prop_ranked)
dev.off()
cat(sprintf("  Saved exp6_vp_species_proportions_%s.pdf\n", param_tag))

# ==============================================================================
# STEP 5: VP STABILITY — SITE-LEVEL MEANS ACROSS FOLDS
# ==============================================================================
cat("\n=== Step 5: VP site-level stability across folds ===\n")

# For each fold, compute mean site-level E/S/C/R2
si_vp_list = list()

for (fi in seq_along(fold_data_list)) {
  fd = fold_data_list[[fi]]
  si = fd$partition$sites

  si_vp_list[[fi]] = tibble(
    fold = fd$fold,
    Overall     = mean(si$r2, na.rm = TRUE),
    Environment = mean(si$env, na.rm = TRUE),
    Spatial     = mean(si$spa, na.rm = TRUE),
    Biotic      = mean(si$codist, na.rm = TRUE)
  )
}

df_si_vp = bind_rows(si_vp_list) %>%
  pivot_longer(-fold, names_to = "component", values_to = "mean_value") %>%
  mutate(component = factor(component, levels = names(unit_colors)))

df_si_summary = df_si_vp %>%
  group_by(component) %>%
  summarise(
    grand_mean = mean(mean_value, na.rm = TRUE),
    sd_val     = sd(mean_value, na.rm = TRUE),
    ci_lo      = grand_mean - qt(0.975, k - 1) * sd_val / sqrt(k),
    ci_hi      = grand_mean + qt(0.975, k - 1) * sd_val / sqrt(k),
    .groups = "drop"
  )

# A) Crossbar plot of means across folds
p_si_means = ggplot(df_si_vp, aes(x = component, y = mean_value, color = component)) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
  geom_crossbar(data = df_si_summary,
                aes(x = component, y = grand_mean, ymin = ci_lo, ymax = ci_hi, fill = component),
                alpha = 0.2, width = 0.4, linewidth = 0.4) +
  scale_color_manual(values = unit_colors) +
  scale_fill_manual(values = unit_colors) +
  geom_text(data = df_si_summary,
            aes(x = component, y = ci_hi,
                label = sprintf("%.4f\n[%.4f]", grand_mean, sd_val)),
            vjust = -0.3, size = 3, show.legend = FALSE) +
  labs(
    title = "A) Mean site VP across folds",
    subtitle = "Dots = per-fold mean, bars = grand mean +/- 95% CI",
    x = NULL, y = "Mean site-level value"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none") +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(15, 5.5, 5.5, 5.5))

# B) Per-site SD across folds
si_per_site_list = list()
for (fi in seq_along(fold_data_list)) {
  fd = fold_data_list[[fi]]
  si = fd$partition$sites
  for (col in c("r2", "env", "spa", "codist")) {
    comp_label = recode(col, r2 = "Overall", env = "Environment",
                        spa = "Spatial", codist = "Biotic")
    si_per_site_list = c(si_per_site_list, list(tibble(
      fold = fd$fold, component = comp_label,
      site_idx = seq_len(nrow(si)),
      value = si[, col]
    )))
  }
}

df_si_per_site = bind_rows(si_per_site_list) %>%
  mutate(component = factor(component, levels = names(unit_colors)))

df_si_site_sd = df_si_per_site %>%
  group_by(component, site_idx) %>%
  summarise(sd_across_folds = sd(value, na.rm = TRUE), .groups = "drop")

p_si_sd = ggplot(df_si_site_sd, aes(x = sd_across_folds, fill = component)) +
  geom_histogram(bins = 30, alpha = 0.7, color = "white") +
  facet_wrap(~ component, scales = "free_y") +
  scale_fill_manual(values = unit_colors) +
  labs(
    title = "B) Per-site SD of VP values across folds",
    subtitle = "How much does each site's partition change across 10 training sets?",
    x = "SD across folds", y = "Count"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

p_si_stability = p_si_means / p_si_sd +
  plot_annotation(
    title = "Experiment 6: Site-level VP stability across folds",
    subtitle = hp_subtitle,
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10)
    )
  )

pdf(here("Calanda_JSDM", "plot", "detailed_modeling_experiment_figures", paste0("exp6_vp_sites_stability_", param_tag, ".pdf")),
    width = 12, height = 12)
print(p_si_stability)
dev.off()
cat(sprintf("  Saved exp6_vp_sites_stability_%s.pdf\n", param_tag))

# --- C) Per-site ranked mean + CI across folds ---
cat("  Building per-site VP error bar plot...\n")

# Count how many folds each site appears in the training set
si_fold_count = integer(n_sites)
for (fi in seq_along(fold_data_list)) {
  fd = fold_data_list[[fi]]
  si_fold_count[fd$train_idx] = si_fold_count[fd$train_idx] + 1L
}

df_si_site_ci = df_si_per_site %>%
  group_by(component, site_idx) %>%
  summarise(
    mean_val = mean(value, na.rm = TRUE),
    sd_val   = sd(value, na.rm = TRUE),
    ci_lo    = mean_val - qt(0.975, n() - 1) * sd_val / sqrt(n()),
    ci_hi    = mean_val + qt(0.975, n() - 1) * sd_val / sqrt(n()),
    .groups  = "drop"
  ) %>%
  mutate(n_folds_train = si_fold_count[site_idx])

cat(sprintf("  Sites in all %d training folds: %d / %d (expected: %d in 9 folds)\n",
            k, sum(si_fold_count == k), n_sites, n_sites))

# Rank sites by mean value within each component
df_si_site_ci = df_si_site_ci %>%
  group_by(component) %>%
  mutate(rank = rank(mean_val, ties.method = "first")) %>%
  ungroup()

p_si_ranked = ggplot(df_si_site_ci,
                     aes(x = rank, y = mean_val, color = component)) +
  geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi, alpha = n_folds_train / k),
                  size = 0.08, linewidth = 0.15) +
  facet_wrap(~ component, scales = "free_y") +
  scale_color_manual(values = unit_colors) +
  scale_alpha_continuous(range = c(0.3, 1), name = paste0("Training folds\n(/ ", k, ")")) +
  labs(
    title = "Experiment 6: Per-site VP values (mean +/- 95% CI across folds)",
    subtitle = paste0("Sites ranked by mean value | opacity = fraction of folds in training set\n", hp_subtitle),
    x = "Site (ranked)", y = "VP value"
  ) +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(face = "bold"))

pdf(here("Calanda_JSDM", "plot", "detailed_modeling_experiment_figures", paste0("exp6_vp_sites_ranked_", param_tag, ".pdf")),
    width = 14, height = 10)
print(p_si_ranked)
dev.off()
cat(sprintf("  Saved exp6_vp_sites_ranked_%s.pdf\n", param_tag))

# --- D) Per-site ranked PROPORTIONS (mean + CI across folds) ---
cat("  Building per-site VP proportion plots...\n")

# Compute proportions per site per fold: E/(E+S+C), S/(E+S+C), C/(E+S+C)
si_prop_list = list()
si_sum_vs_r2 = list()
for (fi in seq_along(fold_data_list)) {
  fd = fold_data_list[[fi]]
  si = fd$partition$sites
  total = si$env + si$spa + si$codist
  si_sum_vs_r2[[fi]] = tibble(
    fold = fd$fold,
    site_idx = seq_len(nrow(si)),
    r2 = si$r2,
    sum_esc = total
  )
  for (comp in c("env", "spa", "codist")) {
    comp_label = recode(comp, env = "Environment", spa = "Spatial", codist = "Biotic")
    si_prop_list = c(si_prop_list, list(tibble(
      fold = fd$fold, component = comp_label,
      site_idx = seq_len(nrow(si)),
      proportion = si[[comp]] / total
    )))
  }
}

# Sanity check: sum(E+S+C) vs R2
df_si_sum_check = bind_rows(si_sum_vs_r2)
si_sum_r2_cor = cor(df_si_sum_check$sum_esc, df_si_sum_check$r2, use = "complete.obs")
si_sum_r2_diff = mean(abs(df_si_sum_check$sum_esc - df_si_sum_check$r2), na.rm = TRUE)
cat(sprintf("  Sanity check (sites): cor(sum_ESC, R2) = %.4f, mean|diff| = %.6f\n",
            si_sum_r2_cor, si_sum_r2_diff))

df_si_prop = bind_rows(si_prop_list) %>%
  mutate(component = factor(component, levels = c("Environment", "Spatial", "Biotic")))

df_si_prop_ci = df_si_prop %>%
  group_by(component, site_idx) %>%
  summarise(
    mean_val = mean(proportion, na.rm = TRUE),
    sd_val   = sd(proportion, na.rm = TRUE),
    ci_lo    = mean_val - qt(0.975, n() - 1) * sd_val / sqrt(n()),
    ci_hi    = mean_val + qt(0.975, n() - 1) * sd_val / sqrt(n()),
    .groups  = "drop"
  ) %>%
  group_by(component) %>%
  mutate(rank = rank(mean_val, ties.method = "first")) %>%
  ungroup()

p_si_prop_ranked = ggplot(df_si_prop_ci,
                          aes(x = rank, y = mean_val, color = component)) +
  geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi),
                  size = 0.08, linewidth = 0.15, alpha = 0.7) +
  facet_wrap(~ component, scales = "free_y") +
  scale_color_manual(values = prop_colors) +
  labs(
    title = "Experiment 6: Per-site VP proportions (mean +/- 95% CI across folds)",
    subtitle = paste0("Proportion = component / (E+S+C) | Sites ranked by mean proportion\n",
                      sprintf("Sanity check: cor(E+S+C, R2) = %.4f, mean|diff| = %.6f\n", si_sum_r2_cor, si_sum_r2_diff),
                      hp_subtitle),
    x = "Site (ranked)", y = "Proportion of explained variance"
  ) +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(face = "bold"))

pdf(here("Calanda_JSDM", "plot", "detailed_modeling_experiment_figures", paste0("exp6_vp_sites_proportions_", param_tag, ".pdf")),
    width = 14, height = 10)
print(p_si_prop_ranked)
dev.off()
cat(sprintf("  Saved exp6_vp_sites_proportions_%s.pdf\n", param_tag))

# ==============================================================================
# STEP 6: SUMMARY TABLE
# ==============================================================================
cat("\n=== Step 6: Summary table ===\n")

df_summary = bind_cols(
  tibble(
    alpha = hp_alpha, lambda = hp_lambda, lambda_env = hp_lambda_env,
    fit_sampling = hp_fit_sampling, k_folds = k,
    n_species = n_species, n_sites = n_sites,
    # Prediction
    species_auc_mean   = mean(species_cv$test_auc, na.rm = TRUE),
    species_auc_median = median(species_cv$test_auc, na.rm = TRUE),
    species_auc_sd     = sd(species_cv$test_auc, na.rm = TRUE),
    species_auc_gt07   = sum(species_cv$test_auc > 0.7, na.rm = TRUE),
    site_logloss_mean   = mean(site_cv$logloss, na.rm = TRUE),
    site_logloss_median = median(site_cv$logloss, na.rm = TRUE),
    site_logloss_sd     = sd(site_cv$logloss, na.rm = TRUE)
  ),
  # Model-level VP means + SDs
  df_model_summary %>%
    pivot_wider(names_from = component,
                values_from = c(mean_R2, sd_R2),
                names_glue = "model_{component}_{.value}"),
  # Species-level VP means + SDs
  df_sp_summary %>%
    pivot_wider(names_from = component,
                values_from = c(grand_mean, sd_val),
                names_glue = "species_{component}_{.value}"),
  # Site-level VP means + SDs
  df_si_summary %>%
    pivot_wider(names_from = component,
                values_from = c(grand_mean, sd_val),
                names_glue = "sites_{component}_{.value}")
)

summary_file = paste0("exp6_cv_summary_", param_tag, ".csv")
write_csv(df_summary,
          here("Calanda_JSDM", "output", "results", summary_file))
cat(sprintf("  Saved %s\n", summary_file))

cat("\n--- Key metrics ---\n")
cat(sprintf("  Species test AUC: mean=%.3f, median=%.3f (%.0f%% > 0.7)\n",
            df_summary$species_auc_mean, df_summary$species_auc_median,
            100 * df_summary$species_auc_gt07 / n_species))
cat(sprintf("  Site log-loss: mean=%.3f, median=%.3f\n",
            df_summary$site_logloss_mean, df_summary$site_logloss_median))
cat("\n  Model-level VP (mean +/- SD across folds):\n")
for (comp in names(unit_colors)) {
  row = df_model_summary %>% filter(component == comp)
  cat(sprintf("    %s: %.4f +/- %.4f\n", comp, row$mean_R2, row$sd_R2))
}
cat("\n  Species-level VP means (grand mean +/- SD across folds):\n")
for (comp in names(unit_colors)) {
  row = df_sp_summary %>% filter(component == comp)
  cat(sprintf("    %s: %.4f +/- %.4f\n", comp, row$grand_mean, row$sd_val))
}
cat("\n  Site-level VP means (grand mean +/- SD across folds):\n")
for (comp in names(unit_colors)) {
  row = df_si_summary %>% filter(component == comp)
  cat(sprintf("    %s: %.4f +/- %.4f\n", comp, row$grand_mean, row$sd_val))
}

cat("\n=== Post Experiment 6 analysis complete ===\n")
