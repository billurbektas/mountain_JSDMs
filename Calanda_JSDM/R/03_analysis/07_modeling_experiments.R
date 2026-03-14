# ==============================================================================
# Script: final_experiment_figures.R
# Purpose: Consolidated final figures for all model experiments (1-6).
#          Produces publication-ready PDFs from cached experiment results.
#
# Figures:
#   1) Exp 1 — Convergence (named models) + Overall anova R2
#   2) Exp 2 — Convergence (named models) + Overall anova R2
#   3) Final selected model convergence
#   4) Exp 3 — Species/site correlations + violin plots with linked medians
#   5) Exp 5 — Anova R2 differences (drop-one)
#   6) Exp 6 — Prediction violins, species/site rank (colored by AUC/logloss),
#              species/site VP proportions
#   7) Merged timing plot (Exp 3 + 4)
#   8) Merged convergence (Exp 1 selected configs)
#
# Requires: R/00_setup/functions_calanda.R
# ==============================================================================

library(tidyverse)
library(conflicted)
library(here)
library(patchwork)
library(ggrepel)
library(pROC)

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::select)

source(here("Calanda_JSDM", "R", "00_setup", "functions_calanda.R"))

# --- Colors ---
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

prop_colors = c("Environment" = color_env, "Spatial" = color_spa, "Biotic" = color_codist)

anova_colors = c(
  "Overall"   = "black",
  "Environment" = color_env,
  "Spatial"     = color_spa,
  "Biotic"      = color_codist,
  "Env x Bio" = "#5aa17f",
  "Env x Spa" = "#c06060",
  "Bio x Spa" = "#8b6508",
  "Shared"    = "grey50"
)

dir.create(here("Calanda_JSDM", "plot"), showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# HELPER: Load exp1 runs
# ==============================================================================
load_exp1 = function() {
  run_files = sort(list.files(
    here("Calanda_JSDM", "output", "results", "runs"),
    pattern = "\\.rds$", full.names = TRUE
  ))
  history_list = list(); metrics_list = list()
  for (f in run_files) {
    run = readRDS(f)
    cfg = run$config
    rid = cfg$run_id; sf = cfg$spatial_form
    a = cfg$alpha_common; l = cfg$lambda_common
    history_list[[rid]] = tibble(
      run_id = rid, spatial_form = sf, alpha = a, lambda = l,
      iteration = seq_len(nrow(run$model$history)),
      loss = as.numeric(run$model$history)
    )
    anova_res = run$partition$anova$results
    r2_mcf = setNames(anova_res$`R2 McFadden`, anova_res$models)
    metrics_list[[rid]] = tibble(
      run_id = rid, spatial_form = sf, alpha = a, lambda = l,
      metric_type = c("overall_R2", "env_R2", "spa_R2", "bio_R2",
                      "env_bio_R2", "env_spa_R2", "bio_spa_R2", "shared_R2"),
      value = c(r2_mcf["Full"], r2_mcf["F_A"], r2_mcf["F_S"], r2_mcf["F_B"],
                r2_mcf["F_AB"], r2_mcf["F_AS"], r2_mcf["F_BS"], r2_mcf["F_ABS"])
    )
  }
  list(
    history = bind_rows(history_list),
    metrics = bind_rows(metrics_list)
  )
}

# ==============================================================================
# HELPER: Load exp2 runs (+ inject exp1 anchors)
# ==============================================================================
load_exp2 = function() {
  run_files = sort(list.files(
    here("Calanda_JSDM", "output", "results", "decoupled"),
    pattern = "\\.rds$", full.names = TRUE
  ))
  history_list = list(); metrics_list = list()
  for (f in run_files) {
    run = readRDS(f)
    cfg = run$config
    rid = cfg$run_id; v = cfg$vary
    lv = cfg$lambda_vary; a = cfg$alpha_common
    history_list[[rid]] = tibble(
      run_id = rid, vary = v, lambda_vary = lv, alpha = a,
      iteration = seq_len(nrow(run$model$history)),
      loss = as.numeric(run$model$history)
    )
    anova_res = run$partition$anova$results
    r2_mcf = setNames(anova_res$`R2 McFadden`, anova_res$models)
    metrics_list[[rid]] = tibble(
      run_id = rid, vary = v, lambda_vary = lv, alpha = a,
      metric_type = c("overall_R2", "env_R2", "spa_R2", "bio_R2",
                      "env_bio_R2", "env_spa_R2", "bio_spa_R2", "shared_R2"),
      value = c(r2_mcf["Full"], r2_mcf["F_A"], r2_mcf["F_S"], r2_mcf["F_B"],
                r2_mcf["F_AB"], r2_mcf["F_AS"], r2_mcf["F_BS"], r2_mcf["F_ABS"])
    )
  }
  df_history = bind_rows(history_list)
  df_metrics = bind_rows(metrics_list)

  # Inject exp1 anchors (lambda = 0.01)
  for (alpha_val in c(0, 0.5, 1)) {
    exp1_file = here("Calanda_JSDM", "output", "results", "runs",
                     paste0("LIN_XY_XY_a", alpha_val, "_l0.01.rds"))
    if (file.exists(exp1_file)) {
      exp1 = readRDS(exp1_file)
      for (v in c("lambda_env", "lambda_sp", "lambda_bio")) {
        rid = paste0("exp1_anchor_", v, "_a", alpha_val)
        df_history = bind_rows(df_history, tibble(
          run_id = rid, vary = v, lambda_vary = 0.01, alpha = alpha_val,
          iteration = seq_len(nrow(exp1$model$history)),
          loss = as.numeric(exp1$model$history)
        ))
        anova_res = exp1$partition$anova$results
        r2_mcf = setNames(anova_res$`R2 McFadden`, anova_res$models)
        df_metrics = bind_rows(df_metrics, tibble(
          run_id = rid, vary = v, lambda_vary = 0.01, alpha = alpha_val,
          metric_type = c("overall_R2", "env_R2", "spa_R2", "bio_R2",
                          "env_bio_R2", "env_spa_R2", "bio_spa_R2", "shared_R2"),
          value = c(r2_mcf["Full"], r2_mcf["F_A"], r2_mcf["F_S"], r2_mcf["F_B"],
                    r2_mcf["F_AB"], r2_mcf["F_AS"], r2_mcf["F_BS"], r2_mcf["F_ABS"])
        ))
      }
    }
  }
  list(history = df_history, metrics = df_metrics)
}

# ==============================================================================
# FIGURE 1: Experiment 1 — Convergence + Anova R2
# ==============================================================================
cat("\n=== Figure 1: Experiment 1 ===\n")
exp1 = load_exp1()

# Model label: spatial_form | alpha | lambda
exp1$history = exp1$history %>%
  mutate(model_label = paste0(spatial_form, " | a=", alpha, " | l=", lambda))
exp1$metrics = exp1$metrics %>%
  mutate(model_label = paste0(spatial_form, " | a=", alpha, " | l=", lambda))

# 1a) Anova R2 (convergence moved to merged figure)
df1_anova = exp1$metrics %>%
  mutate(
    component = recode(metric_type,
      overall_R2 = "Overall", env_R2 = "Environment",
      spa_R2 = "Spatial", bio_R2 = "Biotic",
      env_bio_R2 = "Env x Bio", env_spa_R2 = "Env x Spa",
      bio_spa_R2 = "Bio x Spa", shared_R2 = "Shared"),
    component = factor(component, levels = names(anova_colors)),
    lambda_f = factor(lambda),
    alpha_lab = paste0("alpha = ", alpha)
  )

df1_labels = df1_anova %>% filter(lambda == max(lambda))

p1_anova = ggplot(df1_anova,
                  aes(x = lambda_f, y = value, color = component, group = component)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  geom_text_repel(data = df1_labels, aes(label = component),
                  nudge_x = 0.3, direction = "y", hjust = 0,
                  size = 2.5, segment.size = 0.2, show.legend = FALSE) +
  scale_color_manual(values = anova_colors) +
  facet_grid(spatial_form ~ alpha_lab) +
  labs(title = "Experiment 1: Model-level McFadden R\u00B2",
       x = expression(lambda), y = "McFadden R\u00B2", color = "Component") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold")) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5.5, 40, 5.5, 5.5))

pdf(here("Calanda_JSDM", "plot", "fig_exp1_anova_R2.pdf"), width = 14, height = 8)
print(p1_anova)
dev.off()
cat("  Saved fig_exp1_anova_R2.pdf\n")

# ==============================================================================
# FIGURE 2: Experiment 2 — Convergence + Anova R2
# ==============================================================================
cat("\n=== Figure 2: Experiment 2 ===\n")
exp2 = load_exp2()

vary_labels = c(
  lambda_env = "Vary lambda[env]",
  lambda_sp  = "Vary lambda[spa]",
  lambda_bio = "Vary lambda[bio]"
)

df2_anova = exp2$metrics %>%
  mutate(
    component = recode(metric_type,
      overall_R2 = "Overall", env_R2 = "Environment",
      spa_R2 = "Spatial", bio_R2 = "Biotic",
      env_bio_R2 = "Env x Bio", env_spa_R2 = "Env x Spa",
      bio_spa_R2 = "Bio x Spa", shared_R2 = "Shared"),
    component = factor(component, levels = names(anova_colors)),
    lambda_f = factor(lambda_vary),
    alpha_lab = paste0("alpha = ", alpha)
  )

df2_labels = df2_anova %>% filter(lambda_vary == max(lambda_vary))

p2_anova = ggplot(df2_anova,
                  aes(x = lambda_f, y = value, color = component, group = component)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  geom_text_repel(data = df2_labels, aes(label = component),
                  nudge_x = 0.3, direction = "y", hjust = 0,
                  size = 2.5, segment.size = 0.2, show.legend = FALSE) +
  scale_color_manual(values = anova_colors) +
  facet_grid(alpha_lab ~ vary, labeller = labeller(vary = vary_labels)) +
  labs(title = "Experiment 2: Model-level McFadden R\u00B2 (decoupled lambda)",
       subtitle = "Others anchored at lambda = 0.01",
       x = expression(lambda[varied]), y = "McFadden R\u00B2", color = "Component") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold")) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5.5, 40, 5.5, 5.5))

pdf(here("Calanda_JSDM", "plot", "fig_exp2_anova_R2.pdf"), width = 14, height = 10)
print(p2_anova)
dev.off()
cat("  Saved fig_exp2_anova_R2.pdf\n")

# (Selected model convergence included in merged convergence figure below)

# ==============================================================================
# FIGURE 4: Experiment 3 — Correlations + Violin plots with linked medians
# ==============================================================================
cat("\n=== Figure 4: Experiment 3 ===\n")

# Load anova saturation data
all_anova_files = list.files(
  here("Calanda_JSDM", "output", "results", "anova_saturation"),
  pattern = "^exp3_anova_.*\\.rds$", full.names = TRUE
)
tag_match = regmatches(basename(all_anova_files[1]),
                       regexpr("a[0-9.]+_l[0-9.]+_le[0-9.]+(?=\\.)", basename(all_anova_files[1]), perl = TRUE))
param_tag = tag_match
vp_files = all_anova_files[grepl(param_tag, basename(all_anova_files), fixed = TRUE)]

exp3_data = list()
for (f in vp_files) {
  run = readRDS(f)
  exp3_data[[as.character(run$anova_samples)]] = run
}
anova_sizes = sort(as.integer(names(exp3_data)))
ref_size = max(anova_sizes)
hp = exp3_data[[as.character(anova_sizes[1])]]
hp_subtitle_3 = sprintf("alpha = %s, lambda = %s, lambda_env = %s | fit sampling = 5000",
                         hp$alpha, hp$lambda, hp$lambda_env)

# 4a) Correlations (species + sites, no MAD)
ref_species = exp3_data[[as.character(ref_size)]]$partition$species
ref_sites   = exp3_data[[as.character(ref_size)]]$partition$sites

corr_list = list()
for (ns in anova_sizes) {
  p_sp = exp3_data[[as.character(ns)]]$partition$species
  p_si = exp3_data[[as.character(ns)]]$partition$sites
  for (col in c("r2", "env", "spa", "codist")) {
    comp_label = recode(col, r2 = "Overall", env = "Environment",
                        spa = "Spatial", codist = "Biotic")
    corr_list = c(corr_list, list(tibble(
      anova_samples = ns, component = comp_label,
      species_r = cor(p_sp[, col], ref_species[, col], use = "complete.obs"),
      sites_r   = cor(p_si[, col], ref_sites[, col], use = "complete.obs")
    )))
  }
}
df_corr = bind_rows(corr_list) %>%
  mutate(component = factor(component, levels = names(unit_colors)))

p3_sp_cor = ggplot(df_corr, aes(x = anova_samples, y = species_r, color = component)) +
  geom_hline(yintercept = 1, linewidth = 0.3, linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = 0.99, linewidth = 0.3, linetype = "dotted", color = "grey50") +
  geom_line(linewidth = 0.8) + geom_point(size = 2.5) +
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(values = unit_colors) +
  labs(title = "A) Species-level correlation with reference",
       x = "ANOVA samples (log scale)", y = "Pearson r", color = "Component") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

p3_si_cor = ggplot(df_corr, aes(x = anova_samples, y = sites_r, color = component)) +
  geom_hline(yintercept = 1, linewidth = 0.3, linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = 0.99, linewidth = 0.3, linetype = "dotted", color = "grey50") +
  geom_line(linewidth = 0.8) + geom_point(size = 2.5) +
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(values = unit_colors) +
  labs(title = "B) Site-level correlation with reference",
       x = "ANOVA samples (log scale)", y = "Pearson r", color = "Component") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

p3_corr = (p3_sp_cor | p3_si_cor) +
  plot_annotation(
    title = sprintf("Experiment 3: Correlation with reference (%s samples)", format(ref_size, big.mark = ",")),
    subtitle = hp_subtitle_3,
    theme = theme(plot.title = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(size = 10))
  )

pdf(here("Calanda_JSDM", "plot", "fig_exp3_correlations.pdf"), width = 14, height = 7)
print(p3_corr)
dev.off()
cat("  Saved fig_exp3_correlations.pdf\n")

# 4b) Violin plots of DIFFERENCES from reference, with linked min/max/median
build_diff_violin = function(exp3_data, anova_sizes, ref_size, unit_type = "species") {
  ref_dat = exp3_data[[as.character(ref_size)]]$partition[[unit_type]]
  # Exclude ref from comparison sizes
  comp_sizes = anova_sizes[anova_sizes != ref_size]

  df_list = list()
  for (ns in comp_sizes) {
    dat = exp3_data[[as.character(ns)]]$partition[[unit_type]]
    for (col in c("env", "spa", "codist", "r2")) {
      comp_label = recode(col, r2 = "Overall", env = "Environment",
                          spa = "Spatial", codist = "Biotic")
      df_list = c(df_list, list(tibble(
        anova_samples = ns, component = comp_label,
        idx = seq_len(nrow(dat)),
        diff = abs(dat[, col] - ref_dat[, col])
      )))
    }
  }
  df = bind_rows(df_list) %>%
    mutate(
      component = factor(component, levels = names(unit_colors)),
      anova_samples_f = factor(anova_samples, levels = comp_sizes,
                               labels = scales::comma(comp_sizes))
    )

  # Compute min/max/median per component per anova_samples
  df_stats = df %>%
    group_by(component, anova_samples, anova_samples_f) %>%
    summarise(
      min_val = min(diff, na.rm = TRUE),
      max_val = max(diff, na.rm = TRUE),
      median_val = median(diff, na.rm = TRUE),
      .groups = "drop"
    )

  title_label = ifelse(unit_type == "species", "Species", "Site")

  p = ggplot(df, aes(x = anova_samples_f, y = diff, fill = component)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
    geom_violin(alpha = 0.3, scale = "width", draw_quantiles = 0.5, linewidth = 0.2) +
    # Linked medians
    geom_line(data = df_stats,
              aes(x = anova_samples_f, y = median_val, color = component, group = component),
              linewidth = 0.6, linetype = "solid", inherit.aes = FALSE) +
    geom_point(data = df_stats,
               aes(x = anova_samples_f, y = median_val, color = component),
               size = 1.5, inherit.aes = FALSE) +
    geom_text(data = df_stats,
              aes(x = anova_samples_f, y = median_val, label = sprintf("%.4f", median_val)),
              vjust = -0.8, size = 2.2, inherit.aes = FALSE) +
    # Linked min
    geom_line(data = df_stats,
              aes(x = anova_samples_f, y = min_val, color = component, group = component),
              linewidth = 0.4, linetype = "dotted", inherit.aes = FALSE) +
    geom_text(data = df_stats,
              aes(x = anova_samples_f, y = min_val, label = sprintf("%.3f", min_val)),
              vjust = 1.5, size = 1.8, color = "grey40", inherit.aes = FALSE) +
    # Linked max
    geom_line(data = df_stats,
              aes(x = anova_samples_f, y = max_val, color = component, group = component),
              linewidth = 0.4, linetype = "dotted", inherit.aes = FALSE) +
    geom_text(data = df_stats,
              aes(x = anova_samples_f, y = max_val, label = sprintf("%.3f", max_val)),
              vjust = -0.8, size = 1.8, color = "grey40", inherit.aes = FALSE) +
    scale_fill_manual(values = unit_colors) +
    scale_color_manual(values = unit_colors) +
    facet_wrap(~ component, scales = "free_y") +
    labs(title = sprintf("%s-level: absolute difference from reference (%s samples)",
                         title_label, format(ref_size, big.mark = ",")),
         x = "ANOVA samples", y = "Absolute absolute difference from reference", fill = "Component") +
    theme_bw(base_size = 10) +
    theme(legend.position = "none", strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  p
}

p3_sp_diff = build_diff_violin(exp3_data, anova_sizes, ref_size, "species")
p3_si_diff = build_diff_violin(exp3_data, anova_sizes, ref_size, "sites")

p3_violins = p3_sp_diff / p3_si_diff +
  plot_annotation(
    title = sprintf("Experiment 3: VP absolute difference from reference (%s samples)",
                    format(ref_size, big.mark = ",")),
    subtitle = paste0(hp_subtitle_3, "\nSolid = median, dotted = min/max"),
    theme = theme(plot.title = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(size = 10))
  )

pdf(here("Calanda_JSDM", "plot", "fig_exp3_violins.pdf"), width = 16, height = 14)
print(p3_violins)
dev.off()
cat("  Saved fig_exp3_violins.pdf\n")

# ==============================================================================
# FIGURE 5: Experiment 5 — Drop-one anova diff (reuse existing plot logic)
# ==============================================================================
cat("\n=== Figure 5: Experiment 5 ===\n")

# Load dropone results
dropone_files = sort(list.files(
  here("Calanda_JSDM", "output", "results", "dropone"),
  pattern = "^dropone_.*\\.rds$", full.names = TRUE
))
tag5 = regmatches(basename(dropone_files[1]),
                  regexpr("a[0-9.]+_l[0-9.]+_le[0-9.]+(?=\\.)", basename(dropone_files[1]), perl = TRUE))

base_file = dropone_files[grepl("dropone_BASE_", dropone_files)]
base_run = readRDS(base_file)
base_anova = base_run$partition$anova$results
base_r2 = setNames(base_anova$`R2 McFadden`, base_anova$models)

drop_files = dropone_files[!grepl("dropone_BASE_", dropone_files)]
drop_list = list()
for (f in drop_files) {
  run = readRDS(f)
  pred = run$dropped
  anova_res = run$partition$anova$results
  r2_mcf = setNames(anova_res$`R2 McFadden`, anova_res$models)
  drop_list = c(drop_list, list(tibble(
    predictor = pred,
    Overall = r2_mcf["Full"] - base_r2["Full"],
    Environment = r2_mcf["F_A"] - base_r2["F_A"],
    Spatial = r2_mcf["F_S"] - base_r2["F_S"],
    Biotic = r2_mcf["F_B"] - base_r2["F_B"]
  )))
}

df_drop = bind_rows(drop_list) %>%
  pivot_longer(cols = c(Overall, Environment, Spatial, Biotic),
               names_to = "component", values_to = "diff") %>%
  mutate(component = factor(component, levels = names(unit_colors)))

# Order predictors by Overall diff
pred_order = df_drop %>%
  filter(component == "Overall") %>%
  arrange(diff) %>%
  pull(predictor)
df_drop = df_drop %>% mutate(predictor = factor(predictor, levels = pred_order))

p5 = ggplot(df_drop, aes(x = predictor, y = diff, fill = component)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.8) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  scale_fill_manual(values = unit_colors) +
  labs(title = "Experiment 5: Drop-one predictor — change in McFadden R\u00B2",
       subtitle = sprintf("Base R\u00B2 = %.4f | %s", base_r2["Full"],
                          sprintf("alpha = %s, lambda = %s, lambda_env = %s",
                                  base_run$alpha, base_run$lambda, base_run$lambda_env)),
       x = "Dropped predictor", y = expression(Delta ~ "McFadden" ~ R^2), fill = "Component") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))

pdf(here("Calanda_JSDM", "plot", "fig_exp5_dropone.pdf"), width = 12, height = 7)
print(p5)
dev.off()
cat("  Saved fig_exp5_dropone.pdf\n")

# ==============================================================================
# FIGURE 6: Experiment 6 — Prediction violins + ranked VP colored by AUC/logloss
# ==============================================================================
cat("\n=== Figure 6: Experiment 6 ===\n")

# Load CV data
species_cv = read_csv(here("Calanda_JSDM", "output", "results",
                           paste0("exp6_species_cv_", param_tag, ".csv")),
                      show_col_types = FALSE)
site_cv = read_csv(here("Calanda_JSDM", "output", "results",
                        paste0("exp6_sites_cv_", param_tag, ".csv")),
                   show_col_types = FALSE)
fold_files = sort(list.files(
  here("Calanda_JSDM", "output", "results", "cv"),
  pattern = paste0("^fold_.*_", param_tag, "\\.rds$"), full.names = TRUE
))
k = length(fold_files)
fold_data_list = lapply(fold_files, readRDS)

hp6 = fold_data_list[[1]]
hp_subtitle_6 = sprintf("alpha = %s, lambda = %s, lambda_env = %s | fit sampling = %s | %d-fold CV",
                         hp6$alpha, hp6$lambda, hp6$lambda_env, hp6$fit_sampling, k)

# Load Y matrix
data_files = sort(list.files(
  here("Calanda_JSDM", "output"),
  pattern = "^data_calanda_jsdm_[0-9].*\\.rds$", full.names = TRUE
), decreasing = TRUE)
data_calanda = readRDS(data_files[1])
Y = data_calanda$Y
n_species = ncol(Y); n_sites = nrow(Y)
sp_names = colnames(Y)

# 6a) Prediction violins (with null + random lines, counts)
df_auc = species_cv %>% filter(!is.na(test_auc)) %>%
  transmute(metric = "Species test AUC", value = test_auc)
df_ll = site_cv %>% filter(!is.na(logloss)) %>%
  transmute(metric = "Site log-loss", value = logloss)

eps = 1e-7
folds = readRDS(here("Calanda_JSDM", "output", "results", "cv", "stratified_folds.rds"))
null_logloss = numeric(nrow(Y))
for (fold_i in seq_len(k)) {
  train_idx = setdiff(seq_len(nrow(Y)), folds[[fold_i]])
  test_idx = folds[[fold_i]]
  prev_train = colMeans(Y[train_idx, ])
  p_null = pmin(pmax(prev_train, eps), 1 - eps)
  for (i in test_idx) {
    null_logloss[i] = -mean(Y[i, ] * log(p_null) + (1 - Y[i, ]) * log(1 - p_null))
  }
}
null_ll_mean = mean(null_logloss)

p6_auc = ggplot(df_auc, aes(x = metric, y = value)) +
  geom_violin(fill = color_env, alpha = 0.4, draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 0.8, color = color_env) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey40") +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "grey30") +
  annotate("text", x = 1, y = 0.51, label = "random = 0.5", size = 3, color = "grey40", hjust = -0.1) +
  annotate("text", x = 1, y = 0.71,
           label = sprintf("0.7  (%d/%d above)", sum(df_auc$value > 0.7), nrow(df_auc)),
           size = 3, color = "grey30", hjust = -0.1) +
  labs(title = "A) Species test AUC", x = NULL, y = "AUC") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

p6_ll = ggplot(df_ll, aes(x = metric, y = value)) +
  geom_violin(fill = color_spa, alpha = 0.4, draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter(width = 0.15, alpha = 0.15, size = 0.5, color = color_spa) +
  geom_hline(yintercept = log(2), linetype = "dotted", color = "grey40") +
  annotate("text", x = 1, y = log(2) + 0.01,
           label = sprintf("random = %.3f  (%d/%d below)", log(2),
                           sum(df_ll$value < log(2)), nrow(df_ll)),
           size = 3, color = "grey40", hjust = -0.1) +
  geom_hline(yintercept = null_ll_mean, linetype = "dashed", color = "grey30") +
  annotate("text", x = 0.95, y = null_ll_mean + 0.01,
           label = sprintf("prevalence = %.3f  (%d/%d below)", null_ll_mean,
                           sum(df_ll$value < null_ll_mean), nrow(df_ll)),
           size = 3, color = "grey30", hjust = 0) +
  labs(title = "B) Site log-loss", x = NULL, y = "Log-loss (lower = better)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

p6_violins = (p6_auc | p6_ll) +
  plot_annotation(
    title = "Experiment 6: Predictive performance (out-of-fold)",
    subtitle = hp_subtitle_6,
    theme = theme(plot.title = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(size = 10))
  )

pdf(here("Calanda_JSDM", "plot", "fig_exp6_prediction_violins.pdf"), width = 10, height = 7)
print(p6_violins)
dev.off()
cat("  Saved fig_exp6_prediction_violins.pdf\n")

# 6b) Per-species ranked VP colored by AUC
sp_per_species_list = list()
for (fi in seq_along(fold_data_list)) {
  fd = fold_data_list[[fi]]
  sp = fd$partition$species
  for (col in c("r2", "env", "spa", "codist")) {
    comp_label = recode(col, r2 = "Overall", env = "Environment",
                        spa = "Spatial", codist = "Biotic")
    sp_per_species_list = c(sp_per_species_list, list(tibble(
      fold = fd$fold, component = comp_label,
      species_idx = seq_len(nrow(sp)), value = sp[, col]
    )))
  }
}
df_sp_per_species = bind_rows(sp_per_species_list) %>%
  mutate(component = factor(component, levels = names(unit_colors)))

df_sp_ci = df_sp_per_species %>%
  group_by(component, species_idx) %>%
  summarise(
    mean_val = mean(value, na.rm = TRUE),
    sd_val   = sd(value, na.rm = TRUE),
    ci_lo    = mean_val - qt(0.975, n() - 1) * sd_val / sqrt(n()),
    ci_hi    = mean_val + qt(0.975, n() - 1) * sd_val / sqrt(n()),
    .groups  = "drop"
  ) %>%
  mutate(species_name = sp_names[species_idx])

# Join AUC per species
sp_auc = tibble(species_idx = seq_len(n_species), auc = species_cv$test_auc)
df_sp_ci = df_sp_ci %>% left_join(sp_auc, by = "species_idx")

df_sp_ci = df_sp_ci %>%
  group_by(component) %>%
  mutate(rank = rank(mean_val, ties.method = "first")) %>%
  ungroup()

df_sp_ci = df_sp_ci %>% mutate(low_auc = !is.na(auc) & auc < 0.7)

n_low_auc = n_distinct(df_sp_ci$species_idx[df_sp_ci$low_auc])
n_ok_auc  = n_distinct(df_sp_ci$species_idx[!df_sp_ci$low_auc])

p6_sp_ranked = ggplot(df_sp_ci, aes(x = rank, y = mean_val, color = auc)) +
  geom_pointrange(data = filter(df_sp_ci, low_auc),
                  aes(ymin = ci_lo, ymax = ci_hi), size = 0.08, linewidth = 0.15,
                  color = "black", alpha = 1) +
  geom_pointrange(data = filter(df_sp_ci, !low_auc),
                  aes(ymin = ci_lo, ymax = ci_hi), size = 0.08, linewidth = 0.15, alpha = 0.5) +
  facet_wrap(~ component, scales = "free_y") +
  scale_color_viridis_c(option = "plasma", name = "Species AUC") +
  labs(title = "Experiment 6: Per-species VP values (mean \u00B1 95% CI across folds)",
       subtitle = paste0(sprintf("Black = AUC < 0.7 (%d species) | Colored = AUC >= 0.7 (%d species)\n",
                                  n_low_auc, n_ok_auc), hp_subtitle_6),
       x = "Species (ranked)", y = "VP value") +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(face = "bold"))

pdf(here("Calanda_JSDM", "plot", "fig_exp6_species_ranked.pdf"), width = 14, height = 10)
print(p6_sp_ranked)
dev.off()
cat("  Saved fig_exp6_species_ranked.pdf\n")

# 6c) Per-site ranked VP colored by logloss
si_per_site_list = list()
for (fi in seq_along(fold_data_list)) {
  fd = fold_data_list[[fi]]
  si = fd$partition$sites
  for (col in c("r2", "env", "spa", "codist")) {
    comp_label = recode(col, r2 = "Overall", env = "Environment",
                        spa = "Spatial", codist = "Biotic")
    si_per_site_list = c(si_per_site_list, list(tibble(
      fold = fd$fold, component = comp_label,
      site_idx = seq_len(nrow(si)), value = si[, col]
    )))
  }
}
df_si_per_site = bind_rows(si_per_site_list) %>%
  mutate(component = factor(component, levels = names(unit_colors)))

df_si_ci = df_si_per_site %>%
  group_by(component, site_idx) %>%
  summarise(
    mean_val = mean(value, na.rm = TRUE),
    sd_val   = sd(value, na.rm = TRUE),
    ci_lo    = mean_val - qt(0.975, n() - 1) * sd_val / sqrt(n()),
    ci_hi    = mean_val + qt(0.975, n() - 1) * sd_val / sqrt(n()),
    .groups  = "drop"
  )

# Join logloss per site
si_ll = tibble(site_idx = seq_len(n_sites), logloss = site_cv$logloss)
df_si_ci = df_si_ci %>% left_join(si_ll, by = "site_idx")

df_si_ci = df_si_ci %>%
  group_by(component) %>%
  mutate(rank = rank(mean_val, ties.method = "first")) %>%
  ungroup()

df_si_ci = df_si_ci %>% mutate(high_ll = !is.na(logloss) & logloss > log(2))

n_high_ll = n_distinct(df_si_ci$site_idx[df_si_ci$high_ll])
n_ok_ll   = n_distinct(df_si_ci$site_idx[!df_si_ci$high_ll])

p6_si_ranked = ggplot(df_si_ci, aes(x = rank, y = mean_val, color = logloss)) +
  geom_pointrange(data = filter(df_si_ci, high_ll),
                  aes(ymin = ci_lo, ymax = ci_hi), size = 0.04, linewidth = 0.1,
                  color = "black", alpha = 1) +
  geom_pointrange(data = filter(df_si_ci, !high_ll),
                  aes(ymin = ci_lo, ymax = ci_hi), size = 0.04, linewidth = 0.1, alpha = 0.5) +
  facet_wrap(~ component, scales = "free_y") +
  scale_color_viridis_c(option = "inferno", direction = -1, name = "Site log-loss") +
  labs(title = "Experiment 6: Per-site VP values (mean \u00B1 95% CI across folds)",
       subtitle = paste0(sprintf("Black = log-loss > log(2) (%d sites) | Colored = log-loss <= log(2) (%d sites)\n",
                                  n_high_ll, n_ok_ll), hp_subtitle_6),
       x = "Site (ranked)", y = "VP value") +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(face = "bold"))

pdf(here("Calanda_JSDM", "plot", "fig_exp6_sites_ranked.pdf"), width = 14, height = 10)
print(p6_si_ranked)
dev.off()
cat("  Saved fig_exp6_sites_ranked.pdf\n")

# (Proportion plots removed per user request)

# ==============================================================================
# FIGURE 7: Merged timing plot (Exp 3 + Exp 4)
# ==============================================================================
cat("\n=== Figure 7: Merged timing ===\n")

# Exp3 timing (anova sampling)
df_timing_3 = tibble(
  anova_samples = anova_sizes,
  vp_time_min = sapply(anova_sizes, function(ns) exp3_data[[as.character(ns)]]$vp_time_min)
)

# Exp4 timing (fit sampling)
fit_files = sort(list.files(
  here("Calanda_JSDM", "output", "results", "fit_saturation"),
  pattern = paste0("^exp4_anova_fit.*_", param_tag, "\\.rds$"), full.names = TRUE
))
df_timing_4 = tibble(
  fit_sampling = sapply(fit_files, function(f) readRDS(f)$fit_sampling),
  vp_time_min = sapply(fit_files, function(f) readRDS(f)$vp_time_min)
)

color_exp3 = "#1f77b4"
color_exp4 = "#ff7f0e"

df_timing_merged = bind_rows(
  df_timing_3 %>% transmute(samples = anova_samples, vp_time_min, experiment = "Exp 3: ANOVA sampling"),
  df_timing_4 %>% transmute(samples = fit_sampling, vp_time_min, experiment = "Exp 4: Fit sampling")
)

p_timing = ggplot(df_timing_merged, aes(x = samples, y = vp_time_min, color = experiment)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  geom_text(aes(label = paste0(round(vp_time_min, 1), " min")), vjust = -1, size = 3,
            show.legend = FALSE) +
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(values = c("Exp 3: ANOVA sampling" = color_exp3,
                                "Exp 4: Fit sampling" = color_exp4)) +
  labs(title = "Variance partitioning wall time",
       subtitle = hp_subtitle_3,
       x = "Samples (log scale)", y = "VP time (min)", color = NULL) +
  theme_bw(base_size = 11) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(10, 5.5, 5.5, 5.5), legend.position = "bottom")

pdf(here("Calanda_JSDM", "plot", "fig_merged_timing.pdf"), width = 14, height = 6)
print(p_timing)
dev.off()
cat("  Saved fig_merged_timing.pdf\n")

# ==============================================================================
# FIGURE 8: Merged convergence (all Exp 1 + Exp 2 models, selected in black)
# ==============================================================================
cat("\n=== Figure 8: Merged convergence ===\n")

# Combine exp1 + exp2 histories
df_conv_exp1 = exp1$history %>%
  transmute(experiment = "Exp 1", model_label, iteration, loss)
df_conv_exp2 = exp2$history %>%
  mutate(model_label = paste0(vary, " | lv=", lambda_vary, " | a=", alpha)) %>%
  transmute(experiment = "Exp 2", model_label, iteration, loss)

df_conv_all = bind_rows(df_conv_exp1, df_conv_exp2) %>%
  mutate(uid = paste0(experiment, " | ", model_label))

# Selected model: LIN_XY_XY | a=1 | l=0.01 from Exp 1
selected_uid = df_conv_all %>%
  filter(experiment == "Exp 1",
         model_label == "LIN_XY_XY | a=1 | l=0.01") %>%
  pull(uid) %>% unique()

df_conv_all = df_conv_all %>%
  mutate(is_selected = uid %in% selected_uid)

# Distinct color per model, override selected to black
all_uids = unique(df_conv_all$uid)
n_all = length(all_uids)
set.seed(42)
conv_palette = setNames(scales::hue_pal()(n_all), all_uids)
conv_palette[selected_uid] = "black"

p_merged_conv = ggplot(df_conv_all, aes(x = iteration, y = loss,
                                         group = uid, color = uid)) +
  geom_line(data = filter(df_conv_all, !is_selected),
            alpha = 0.35, linewidth = 0.3) +
  geom_line(data = filter(df_conv_all, is_selected),
            linewidth = 0.9) +
  scale_color_manual(values = conv_palette, name = "Model") +
  labs(title = "Training convergence: all Exp 1 + Exp 2 models",
       subtitle = "Selected final model in black",
       x = "Iteration", y = "Loss") +
  theme_bw(base_size = 11) +
  theme(legend.text = element_text(size = 5), legend.key.height = unit(0.35, "cm"))

pdf(here("Calanda_JSDM", "plot", "fig_merged_convergence.pdf"), width = 14, height = 8)
print(p_merged_conv)
dev.off()
cat("  Saved fig_merged_convergence.pdf\n")

cat("\n=== All final figures complete ===\n")
