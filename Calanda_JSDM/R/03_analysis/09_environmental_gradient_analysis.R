# ==============================================================================
# Script: 09_environmental_gradient_analysis.R
# Author: Billur Bektas
# Claude (Anthropic) was used to assist with code refactoring, validation, and documentation.
#
# Purpose: Community-level (site) and species-level environmental gradient
#          analysis using 10-fold CV results.
#
#   COMMUNITY (sites):
#     - Stepwise BIC regressions of VP vs all 12 env variables (unweighted)
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
#     - Sites: logloss <= log(2), ci_width < 0.1 (unweighted, AIC selection)
#     - Species: AUC >= 0.7, weighted by 1 / CI_width^2 (AIC selection)
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
#   - plot/env_gradient_climate.pdf
#   - plot/env_gradient_other.pdf
#
# Requires: R/00_setup/functions_calanda.R, sjSDM
# ==============================================================================

library(tidyverse)
library(conflicted)
library(here)
library(patchwork)
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
vp_labels = c("env" = "Environment", "codist" = "Species Associations", "spa" = "Space")

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

# Filter by logloss using site names (not positional indices)
site_logloss = site_cv %>% select(site, logloss)
keep_names = site_logloss %>%
  filter(logloss <= log(2), !is.na(logloss)) %>%
  pull(site)
n_keep = length(keep_names)
cat(sprintf("  %d / %d sites kept (logloss <= log(2))\n", n_keep, n_sites))

df_wide = df_si_summary %>%
  filter(site_name %in% keep_names) %>%
  filter(component %in% c("env", "spa", "codist", "r2")) %>%
  select(site_name, site_idx, component, mean_val, ci_width) %>%
  pivot_wider(
    names_from = component,
    values_from = c(mean_val, ci_width),
    names_glue = "{component}_{.value}"
  )

# Remove sites with unstable VP estimates (ci_width >= 0.1 on any component)
ci_cols = grep("_ci_width$", names(df_wide), value = TRUE)
unstable = apply(df_wide[, ci_cols], 1, function(r) any(r >= 0.1, na.rm = TRUE))
if (any(unstable)) {
  cat(sprintf("  Removing %d site(s) with ci_width >= 0.1\n", sum(unstable)))
  df_wide = df_wide[!unstable, ]
}

# Join environmental data by site name (not positional index into X)
X_df = as.data.frame(X) %>%
  rownames_to_column("site_name") %>%
  select(site_name, slope, summer_temp, fdd, et_annual = et.annual,
         soil_depth_var, trees_cover, shrubs_cover, rocks_cover,
         flowdir, tpi, nutrient, disturbance)

df_wide = df_wide %>%
  left_join(X_df, by = "site_name")

cat(sprintf("  Wide data: %d sites x %d columns\n", nrow(df_wide), ncol(df_wide)))

# ==============================================================================
# A1. WEIGHTED STEPWISE: VP ~ all 12 env variables
# ==============================================================================
cat("\n--- VP vs all 12 env variables (stepwise AIC, unweighted) ---\n")

comm_linear = env_cols
comm_quad = paste0("I(", env_cols, "^2)")
comm_rhs = paste(c(comm_linear, comm_quad), collapse = " + ")

community_env_models = list()
comm_env_results = list()

for (comp in vp_components) {
  y_col = paste0(comp, "_mean_val")

  full_formula = as.formula(paste0(y_col, " ~ ", comm_rhs))
  full_mod = lm(full_formula, data = df_wide)
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

  comm_env_results[[comp]] = extract_model_results(full_mod, best_mod, "env_gradient_community", comp)
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
sp_beta_results = list()

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

  sp_beta_results[[comp]] = extract_model_results(full_mod, best_mod, "env_gradient_species", comp)
}


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  PART C: COMBINED FIGURES — species (row 1) + community (row 2)         ║
# ╚════════════════════════════════════════════════════════════════════════════╝
cat("\n========== COMBINED FIGURES ==========\n")

# --- Climate variables: 3 rows x 2 cols (species | community) ---
build_combined_pdf(climate_vars, "env_gradient_climate.pdf",
                   pdf_width = 10, pdf_height = length(climate_vars) * 4)

# --- Land cover variables ---
cover_vars = c("rocks_cover", "trees_cover", "shrubs_cover", "soil_depth_var")
build_combined_pdf(cover_vars, "env_gradient_cover.pdf",
                   pdf_width = 10, pdf_height = length(cover_vars) * 4)

# --- Topography variables ---
topo_vars = c("slope", "tpi", "flowdir")
build_combined_pdf(topo_vars, "env_gradient_topography.pdf",
                   pdf_width = 10, pdf_height = length(topo_vars) * 4)

# --- Biotic index variables ---
index_vars = c("nutrient", "disturbance")
build_combined_pdf(index_vars, "env_gradient_indices.pdf",
                   pdf_width = 10, pdf_height = length(index_vars) * 4)

# ==============================================================================
# SAVE MODEL RESULTS AS CSV
# ==============================================================================
cat("\n--- Saving model results ---\n")

all_results = c(comm_env_results, sp_beta_results)

# Coefficients
all_coefs = bind_rows(lapply(all_results, function(r) r$coefficients))
write_csv(all_coefs, here("Calanda_JSDM", "output",
                           paste0("env_gradient_coefficients_", param_tag, ".csv")))
cat(sprintf("  Saved env_gradient_coefficients_%s.csv (%d rows)\n", param_tag, nrow(all_coefs)))

# AIC steps
all_aic = bind_rows(lapply(all_results, function(r) r$bic_steps))
write_csv(all_aic, here("Calanda_JSDM", "output",
                         paste0("env_gradient_aic_steps_", param_tag, ".csv")))
cat(sprintf("  Saved env_gradient_aic_steps_%s.csv (%d rows)\n", param_tag, nrow(all_aic)))

# Model summaries
all_summaries = bind_rows(lapply(all_results, function(r) r$model_summary))
write_csv(all_summaries, here("Calanda_JSDM", "output",
                               paste0("env_gradient_model_summary_", param_tag, ".csv")))
cat(sprintf("  Saved env_gradient_model_summary_%s.csv (%d rows)\n", param_tag, nrow(all_summaries)))

cat("\n=== Environmental gradient analysis complete ===\n")
