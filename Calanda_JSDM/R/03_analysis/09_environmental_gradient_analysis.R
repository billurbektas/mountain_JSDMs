# ==============================================================================
# Script: 09_environmental_gradient_analysis.R
# Purpose: Community-level (site) and species-level environmental gradient
#          analysis using 10-fold CV results.
#
#   COMMUNITY (sites):
#     - Weighted regressions of VP components vs altitude
#     - Weighted stepwise AIC regressions of VP vs all 12 env variables
#
#   SPECIES:
#     - Extract species betas from final model
#     - Weighted stepwise AIC regressions of VP vs all 12 betas
#
#   COMBINED FIGURES:
#     - Climate variables (summer_temp, et_annual, fdd, disturbance, nutrient):
#       species row 1, community row 2 → env_gradient_climate.pdf
#     - Other variables (remaining 7):
#       species row 1, community row 2 → env_gradient_other.pdf
#
#   Filtering:
#     - Sites: logloss <= log(2)
#     - Species: AUC >= 0.7
#     - Weights = 1 / CI_width^2
#
# Inputs:
#   - output/results/vp_sites_summary_<param_tag>.csv
#   - output/results/vp_species_summary_<param_tag>.csv
#   - output/results/exp6_sites_cv_<param_tag>.csv
#   - output/results/exp6_species_cv_<param_tag>.csv
#   - output/data_calanda_jsdm_2026-03-06.rds
#   - output/final_model_se_<param_tag>.rds
#
# Outputs:
#   - plot/community_altitude_regressions.pdf
#   - plot/env_gradient_climate.pdf
#   - plot/env_gradient_other.pdf
#
# Requires: R/00_setup/functions_calanda.R, sjSDM
# ==============================================================================

library(tidyverse)
library(conflicted)
library(here)
library(patchwork)
library(ggrepel)
library(sjSDM)

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::select)

source(here("Calanda_JSDM", "R", "00_setup", "functions_calanda.R"))

# --- Colors & labels ---
color_env    = "#81caf3"
color_spa    = "#d00000"
color_codist = "#00bd89"
comp_colors  = c("env" = color_env, "codist" = color_codist, "spa" = color_spa)

vp_components = c("env", "codist", "spa")
vp_labels = c("env" = "Environment", "codist" = "Biotic", "spa" = "Spatial")

# Full variable names
var_full_names = c(
  "summer_temp"    = "Summer Temperature",
  "fdd"            = "Freezing Degree Days",
  "et.annual"      = "Annual Evapotranspiration",
  "et_annual"      = "Annual Evapotranspiration",
  "slope"          = "Slope",
  "rocks_cover"    = "Rock Cover",
  "trees_cover"    = "Tree Cover",
  "shrubs_cover"   = "Shrub Cover",
  "soil_depth_var" = "Soil Depth Variability",
  "tpi"            = "Topographic Position Index",
  "flowdir"        = "Flow Direction",
  "nutrient"       = "Nutrient Index",
  "disturbance"    = "Disturbance Index"
)

# The 12 model env predictors
env_cols = c("summer_temp", "fdd", "et_annual", "slope", "rocks_cover",
             "trees_cover", "shrubs_cover", "soil_depth_var",
             "tpi", "flowdir", "nutrient", "disturbance")

# Variable groups
climate_vars = c("summer_temp", "fdd", "et_annual")
other_vars = setdiff(env_cols, climate_vars)

# ==============================================================================
# HELPER: build a partial effect plot for one variable + one level
# ==============================================================================
build_partial_plot = function(evar, models, df, predictor_cols,
                               vp_components, vp_labels, comp_colors,
                               var_full_names, x_label_prefix = "") {

  pred_parts = list()
  annot_parts = list()

  for (comp in vp_components) {
    best_mod = models[[comp]]
    s = summary(best_mod)
    retained = names(coef(best_mod))[-1]

    relevant = retained[grepl(paste0("^", evar, "$|^I\\(", evar), retained)]
    if (length(relevant) == 0) {
      annot_parts[[comp]] = sprintf("%s: not retained", vp_labels[comp])
      next
    }

    pred_df = as.data.frame(lapply(predictor_cols, function(v) rep(mean(df[[v]]), 200)))
    names(pred_df) = predictor_cols
    pred_df[[evar]] = seq(min(df[[evar]]), max(df[[evar]]), length.out = 200)
    pred_df$y = predict(best_mod, pred_df)
    pred_df$component = comp
    pred_parts[[comp]] = pred_df

    coef_info = as.data.frame(s$coefficients) %>%
      rownames_to_column("term") %>%
      filter(grepl(paste0("^", evar, "$|^I\\(", evar), term)) %>%
      mutate(
        short = ifelse(grepl("\\^2", term), "x\u00B2", "x"),
        sig = ifelse(`Pr(>|t|)` < 0.001, "***",
              ifelse(`Pr(>|t|)` < 0.01, "**",
              ifelse(`Pr(>|t|)` < 0.05, "*", ""))),
        label = sprintf("%s = %.4f%s", short, Estimate, sig)
      )
    annot_parts[[comp]] = sprintf("%s: R2=%.3f, %s",
                                   vp_labels[comp], s$adj.r.squared,
                                   paste(coef_info$label, collapse = ", "))
  }

  df_preds = bind_rows(pred_parts)

  # Points in long format
  df_pts = df %>%
    select(all_of(paste0(vp_components, "_mean_val")), all_of(evar)) %>%
    pivot_longer(cols = all_of(paste0(vp_components, "_mean_val")),
                 names_to = "component", values_to = "value") %>%
    mutate(component = gsub("_mean_val", "", component))

  annot_text = paste(annot_parts[vp_components], collapse = "\n")
  full_name = ifelse(evar %in% names(var_full_names), var_full_names[evar], evar)
  x_lab = if (x_label_prefix != "") paste0(x_label_prefix, full_name) else full_name

  p = ggplot() +
    geom_point(data = df_pts,
               aes(x = .data[[evar]], y = value, color = component),
               size = 0.3, alpha = 0.3) +
    geom_line(data = df_preds,
              aes(x = .data[[evar]], y = y, color = component),
              linewidth = 1.2) +
    scale_color_manual(values = comp_colors, labels = vp_labels, name = NULL) +
    annotate("text",
             x = min(df[[evar]]),
             y = max(df_pts$value, na.rm = TRUE),
             label = annot_text,
             hjust = 0, vjust = 1, size = 4, fontface = "italic") +
    labs(x = x_lab, y = "Variance proportion") +
    theme_bw(base_size = 12) +
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 9),
          legend.position = "none")

  return(p)
}

# ==============================================================================
# LOAD SHARED DATA
# ==============================================================================
cat("\n=== Environmental gradient analysis ===\n")

# Detect param_tag
fold_files = list.files(
  here("Calanda_JSDM", "output", "results", "cv"),
  pattern = "^fold_.*\\.rds$"
)
param_tag = regmatches(fold_files[1],
                       regexpr("a[0-9.]+_l[0-9.]+_le[0-9.]+(?=\\.)",
                               fold_files[1], perl = TRUE))

# Environmental data
data_calanda = readRDS(here("Calanda_JSDM", "output", "data_calanda_jsdm_2026-03-06.rds"))
X = data_calanda$X
n_sites = nrow(X)

# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  PART A: COMMUNITY-LEVEL (SITES)                                        ║
# ╚════════════════════════════════════════════════════════════════════════════╝
cat("\n========== COMMUNITY-LEVEL ==========\n")

df_si_summary = read_csv(
  here("Calanda_JSDM", "output", "results",
       paste0("vp_sites_summary_", param_tag, ".csv")),
  show_col_types = FALSE
)

site_cv = read_csv(
  here("Calanda_JSDM", "output", "results",
       paste0("exp6_sites_cv_", param_tag, ".csv")),
  show_col_types = FALSE
)

si_keep = site_cv$logloss <= log(2) & !is.na(site_cv$logloss)
keep_idx = which(si_keep)
n_keep = length(keep_idx)
cat(sprintf("  %d / %d sites kept (logloss <= log(2))\n", n_keep, n_sites))

df_wide = df_si_summary %>%
  filter(site_idx %in% keep_idx) %>%
  filter(component %in% c("env", "spa", "codist", "r2")) %>%
  select(site_idx, component, mean_val, ci_width, weight) %>%
  pivot_wider(
    names_from = component,
    values_from = c(mean_val, ci_width, weight),
    names_glue = "{component}_{.value}"
  )

df_wide = df_wide %>%
  mutate(
    altitude       = X[site_idx, "altitude"],
    slope          = X[site_idx, "slope"],
    summer_temp    = X[site_idx, "summer_temp"],
    fdd            = X[site_idx, "fdd"],
    et_annual      = X[site_idx, "et.annual"],
    soil_depth_var = X[site_idx, "soil_depth_var"],
    trees_cover    = X[site_idx, "trees_cover"],
    shrubs_cover   = X[site_idx, "shrubs_cover"],
    rocks_cover    = X[site_idx, "rocks_cover"],
    flowdir        = X[site_idx, "flowdir"],
    tpi            = X[site_idx, "tpi"],
    nutrient       = X[site_idx, "nutrient"],
    disturbance    = X[site_idx, "disturbance"]
  )

cat(sprintf("  Wide data: %d sites x %d columns\n", nrow(df_wide), ncol(df_wide)))

# ==============================================================================
# A1. ALTITUDE REGRESSIONS (weighted)
# ==============================================================================
cat("\n--- VP components vs altitude (weighted) ---\n")

altitude_models = list()
altitude_preds = list()

for (comp in vp_components) {
  y_col = paste0(comp, "_mean_val")
  w_col = paste0(comp, "_weight")

  full_formula = as.formula(paste0(y_col, " ~ altitude + I(altitude^2)"))
  full_mod = lm(full_formula, data = df_wide, weights = df_wide[[w_col]])
  best_mod = step(full_mod, direction = "backward", trace = 0)
  altitude_models[[comp]] = best_mod

  s = summary(best_mod)
  cat(sprintf("  %s: adj.R2 = %.3f, p = %.2e\n", vp_labels[comp],
              s$adj.r.squared, pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3],
                                  lower.tail = FALSE)))
  print(s$coefficients)
  cat("\n")

  pred_df = data.frame(altitude = seq(min(df_wide$altitude), max(df_wide$altitude),
                                       length.out = 200))
  pred_df$y = predict(best_mod, pred_df)
  pred_df$component = comp
  altitude_preds[[comp]] = pred_df
}

df_alt_preds = bind_rows(altitude_preds)

p_alt_list = list()
for (comp in vp_components) {
  y_col = paste0(comp, "_mean_val")
  w_col = paste0(comp, "_weight")
  preds = df_alt_preds %>% filter(component == comp)
  s = summary(altitude_models[[comp]])

  coef_df = as.data.frame(s$coefficients) %>%
    rownames_to_column("term") %>%
    filter(term != "(Intercept)") %>%
    mutate(
      short = case_when(
        grepl("I\\(", term) ~ "altitude2",
        TRUE ~ "altitude"
      ),
      sig = ifelse(`Pr(>|t|)` < 0.001, "***",
            ifelse(`Pr(>|t|)` < 0.01, "**",
            ifelse(`Pr(>|t|)` < 0.05, "*", ""))),
      label = sprintf("%s: %.4f%s", short, Estimate, sig)
    )
  effect_label = paste(coef_df$label, collapse = "\n")

  p = ggplot(df_wide, aes(x = altitude, y = .data[[y_col]])) +
    geom_point(aes(size = .data[[w_col]]), alpha = 0.4, color = comp_colors[comp]) +
    geom_line(data = preds, aes(x = altitude, y = y), linewidth = 1, color = "black") +
    scale_size_continuous(range = c(0.5, 3), guide = "none") +
    annotate("text", x = min(df_wide$altitude), y = max(df_wide[[y_col]], na.rm = TRUE),
             label = sprintf("adj.R2 = %.3f\n%s", s$adj.r.squared, effect_label),
             hjust = 0, vjust = 1, size = 4.5, fontface = "italic") +
    labs(x = "Altitude (scaled)", y = vp_labels[comp]) +
    theme_bw(base_size = 16) +
    theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14))

  p_alt_list[[comp]] = p
}

p_altitude = p_alt_list[["env"]] + p_alt_list[["codist"]] + p_alt_list[["spa"]] +
  plot_layout(nrow = 1)

pdf(here("Calanda_JSDM", "plot", "community_altitude_regressions.pdf"), width = 18, height = 6)
print(p_altitude)
dev.off()
cat("  Saved community_altitude_regressions.pdf\n")

# ==============================================================================
# A2. WEIGHTED STEPWISE: VP ~ all 12 env variables
# ==============================================================================
cat("\n--- VP vs all 12 env variables (weighted stepwise AIC) ---\n")

comm_linear = env_cols
comm_quad = paste0("I(", env_cols, "^2)")
comm_rhs = paste(c(comm_linear, comm_quad), collapse = " + ")

community_env_models = list()

for (comp in vp_components) {
  y_col = paste0(comp, "_mean_val")
  w_col = paste0(comp, "_weight")

  full_formula = as.formula(paste0(y_col, " ~ ", comm_rhs))
  full_mod = lm(full_formula, data = df_wide, weights = df_wide[[w_col]])
  best_mod = step(full_mod, direction = "backward", trace = 0)
  community_env_models[[comp]] = best_mod

  s = summary(best_mod)
  cat(sprintf("\n  %s: adj.R2 = %.3f, p = %.2e\n", vp_labels[comp],
              s$adj.r.squared, pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3],
                                  lower.tail = FALSE)))
  coef_df = as.data.frame(s$coefficients) %>%
    rownames_to_column("term") %>%
    filter(term != "(Intercept)")
  print(s$coefficients)
  cat(sprintf("  Retained %d / %d terms\n", nrow(coef_df), length(c(comm_linear, comm_quad))))
}


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  PART B: SPECIES-LEVEL                                                  ║
# ╚════════════════════════════════════════════════════════════════════════════╝
cat("\n========== SPECIES-LEVEL ==========\n")

model_file = here("Calanda_JSDM", "output",
                   paste0("final_model_se_", param_tag, ".rds"))
fm = readRDS(model_file)
environment(fm$get_model)$device = "cpu"
cat(sprintf("  Loaded final model: %s\n", basename(model_file)))

model_summary = summary(fm)
species_betas = as.data.frame(t(model_summary$coefs[-1, ]))
species_betas$species_name = fm$species
rownames(species_betas) = NULL

# Rename et.annual to et_annual to match community env_cols
if ("et.annual" %in% names(species_betas)) {
  names(species_betas)[names(species_betas) == "et.annual"] = "et_annual"
}

beta_cols = setdiff(names(species_betas), "species_name")
cat(sprintf("  Beta matrix: %d species x %d env predictors\n",
            nrow(species_betas), length(beta_cols)))

df_sp_summary = read_csv(
  here("Calanda_JSDM", "output", "results",
       paste0("vp_species_summary_", param_tag, ".csv")),
  show_col_types = FALSE
)

species_cv = read_csv(
  here("Calanda_JSDM", "output", "results",
       paste0("exp6_species_cv_", param_tag, ".csv")),
  show_col_types = FALSE
)

sp_keep = species_cv$test_auc >= 0.7 & !is.na(species_cv$test_auc)
n_sp_keep = sum(sp_keep)
cat(sprintf("  %d / %d species kept (AUC >= 0.7)\n", n_sp_keep, nrow(species_cv)))

keep_names = species_cv$species[sp_keep]

df_sp_wide = df_sp_summary %>%
  filter(species_name %in% keep_names) %>%
  filter(component %in% c("env", "spa", "codist", "r2")) %>%
  select(species_idx, species_name, component, mean_val, ci_width, weight) %>%
  pivot_wider(
    names_from = component,
    values_from = c(mean_val, ci_width, weight),
    names_glue = "{component}_{.value}"
  )

df_merged = df_sp_wide %>%
  inner_join(species_betas, by = "species_name")

cat(sprintf("  Merged data: %d species with betas + VP\n", nrow(df_merged)))

# ==============================================================================
# B1. WEIGHTED STEPWISE: VP ~ all 12 betas
# ==============================================================================
cat("\n--- VP vs species betas (weighted stepwise AIC) ---\n")

linear_terms = beta_cols
quad_terms = paste0("I(", beta_cols, "^2)")
all_terms = c(linear_terms, quad_terms)
rhs = paste(all_terms, collapse = " + ")

beta_models = list()

for (comp in vp_components) {
  y_col = paste0(comp, "_mean_val")
  w_col = paste0(comp, "_weight")

  full_formula = as.formula(paste0(y_col, " ~ ", rhs))
  full_mod = lm(full_formula, data = df_merged, weights = df_merged[[w_col]])
  best_mod = step(full_mod, direction = "backward", trace = 0)
  beta_models[[comp]] = best_mod

  s = summary(best_mod)
  cat(sprintf("\n  %s: adj.R2 = %.3f, p = %.2e\n", vp_labels[comp],
              s$adj.r.squared, pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3],
                                  lower.tail = FALSE)))
  coef_df = as.data.frame(s$coefficients) %>%
    rownames_to_column("term") %>%
    filter(term != "(Intercept)")
  print(s$coefficients)
  cat(sprintf("  Retained %d / %d terms\n", nrow(coef_df), length(all_terms)))
}


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  PART C: COMBINED FIGURES — species (row 1) + community (row 2)         ║
# ╚════════════════════════════════════════════════════════════════════════════╝
cat("\n========== COMBINED FIGURES ==========\n")

# --- Helper to build a combined PDF for a set of variables ---
build_combined_pdf = function(vars, filename, pdf_width, pdf_height) {

  sp_plots = list()
  comm_plots = list()

  for (v in vars) {
    sp_plots[[v]] = build_partial_plot(
      evar = v, models = beta_models, df = df_merged,
      predictor_cols = beta_cols,
      vp_components = vp_components, vp_labels = vp_labels,
      comp_colors = comp_colors, var_full_names = var_full_names,
      x_label_prefix = "Response to "
    )

    comm_plots[[v]] = build_partial_plot(
      evar = v, models = community_env_models, df = df_wide,
      predictor_cols = env_cols,
      vp_components = vp_components, vp_labels = vp_labels,
      comp_colors = comp_colors, var_full_names = var_full_names,
      x_label_prefix = ""
    )
  }

  n = length(vars)

  # Column headers via ggtitle on first row only
  sp_plots[[1]] = sp_plots[[1]] + ggtitle("A) Species") +
    theme(plot.title = element_text(size = 18, face = "bold", hjust = 0))
  comm_plots[[1]] = comm_plots[[1]] + ggtitle("B) Community") +
    theme(plot.title = element_text(size = 18, face = "bold", hjust = 0))

  # Interleave: sp1, comm1, sp2, comm2, ...
  paired = list()
  for (i in seq_along(vars)) {
    paired[[2 * i - 1]] = sp_plots[[i]]
    paired[[2 * i]]     = comm_plots[[i]]
  }

  # Enable legend on just one plot so collect works
  paired[[1]] = paired[[1]] +
    theme(legend.position = "bottom", legend.text = element_text(size = 12))

  p_combined = wrap_plots(paired, ncol = 2) +
    plot_layout(guides = "collect") +
    plot_annotation(theme = theme(legend.position = "bottom"))

  pdf(here("Calanda_JSDM", "plot", filename), width = pdf_width, height = pdf_height)
  print(p_combined)
  dev.off()
  cat(sprintf("  Saved %s\n", filename))
}

# --- Climate variables: 3 rows x 2 cols (species | community) ---
build_combined_pdf(climate_vars, "env_gradient_climate.pdf",
                   pdf_width = 10, pdf_height = length(climate_vars) * 4)

# --- Other variables: 9 rows x 2 cols (species | community) ---
build_combined_pdf(other_vars, "env_gradient_other.pdf",
                   pdf_width = 10, pdf_height = length(other_vars) * 4)

# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  PART D: EFFECT SIZE FOREST PLOT — all variables, species + community   ║
# ╚════════════════════════════════════════════════════════════════════════════╝
cat("\n========== EFFECT SIZE PLOT ==========\n")

# Helper: extract coefficients from a model list
extract_coefs = function(models, vp_components, vp_labels, level_label) {
  coef_list = list()
  for (comp in vp_components) {
    s = summary(models[[comp]])
    cc = as.data.frame(s$coefficients) %>%
      rownames_to_column("term") %>%
      filter(term != "(Intercept)") %>%
      mutate(
        component = comp,
        level = level_label,
        # Parse variable name and term type
        base_var = gsub("I\\(|\\^2\\)", "", term),
        term_type = ifelse(grepl("\\^2", term), "quadratic", "linear"),
        estimate = Estimate,
        se = `Std. Error`,
        p_value = `Pr(>|t|)`,
        ci_lo = estimate - 1.96 * se,
        ci_hi = estimate + 1.96 * se,
        significant = p_value < 0.05
      )
    coef_list[[comp]] = cc
  }
  bind_rows(coef_list)
}

df_coefs_comm = extract_coefs(community_env_models, vp_components, vp_labels, "Community")
df_coefs_sp   = extract_coefs(beta_models, vp_components, vp_labels, "Species")

df_coefs = bind_rows(df_coefs_comm, df_coefs_sp) %>%
  mutate(
    component = factor(component, levels = vp_components, labels = vp_labels),
    level = factor(level, levels = c("Species", "Community")),
    var_label = ifelse(base_var %in% names(var_full_names),
                       var_full_names[base_var], base_var),
    var_label = factor(var_label, levels = rev(unique(var_label))),
    term_type = factor(term_type, levels = c("linear", "quadratic"))
  )

cat(sprintf("  %d total coefficients (%d community, %d species)\n",
            nrow(df_coefs), nrow(df_coefs_comm), nrow(df_coefs_sp)))

# Dodge position for community vs species side by side
p_forest = ggplot(df_coefs,
                  aes(x = estimate, y = var_label, color = component,
                      shape = term_type, alpha = significant)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi),
                width = 0.3, linewidth = 0.6,
                position = position_dodge(width = 0.7),
                orientation = "y") +
  geom_point(size = 2.5,
             position = position_dodge(width = 0.7)) +
  scale_color_manual(values = c("Environment" = color_env,
                                "Biotic" = color_codist,
                                "Spatial" = color_spa),
                     name = NULL) +
  scale_shape_manual(values = c("linear" = 16, "quadratic" = 17),
                     name = "Term type") +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3), guide = "none") +
  facet_wrap(~ level, nrow = 1) +
  labs(x = "Coefficient estimate (+/- 1.96 SE)", y = NULL) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 13),
    legend.position = "bottom",
    legend.text = element_text(size = 11),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 13, face = "bold")
  )

pdf(here("Calanda_JSDM", "plot", "env_effect_sizes.pdf"),
    width = 16, height = max(8, length(unique(df_coefs$var_label)) * 0.6))
print(p_forest)
dev.off()
cat("  Saved env_effect_sizes.pdf\n")

cat("\n=== Environmental gradient analysis complete ===\n")
