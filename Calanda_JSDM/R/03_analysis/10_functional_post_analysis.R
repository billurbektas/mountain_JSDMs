# ==============================================================================
# Script: 10_functional_post_analysis.R
# Author: Billur Bektas
# Claude (Anthropic) was used to assist with code refactoring, validation, and documentation.
#
# Purpose: Functional trait post-JSDM analysis.
#
#   5 core traits: vegetative height, LNC, LA, seed mass, LDMC
#
#   SPECIES-LEVEL:
#     - Trait medians (log-transformed) + kCVs as predictors
#     - Weighted stepwise AIC regressions of VP vs all trait predictors
#     - Median and kCV shown in same facet (different line types)
#
#   COMMUNITY-LEVEL:
#     - Unweighted community trait means + variances as predictors
#     - Weighted stepwise AIC regressions of VP vs community trait predictors
#     - Mean and variance shown in same facet (different line types)
#
#   Filtering:
#     - Sites: logloss <= log(2), ci_width < 0.1 (unweighted, AIC selection)
#     - Species: AUC >= 0.7, weighted by 1 / CI_width^2 (AIC selection)
#     - To switch to BIC: change selection_k and selection_k_comm to log(n)
#
# Inputs:
#   - output/traits_medians_imputed.csv
#   - output/traits_kcv_imputed.csv
#   - output/community_traits_unweighted_imputed.csv
#   - output/results/vp_species_summary_<param_tag>.csv
#   - output/results/vp_sites_summary_<param_tag>.csv
#   - output/results/exp6_species_cv_<param_tag>.csv
#   - output/results/exp6_sites_cv_<param_tag>.csv
#   - output/data_calanda_jsdm_2026-03-06.rds
#
# Outputs:
#   - plot/functional_gradient.pdf
#
# Requires: R/00_setup/functions_calanda.R
# ==============================================================================

library(tidyverse)
library(conflicted)
library(here)
library(ggpubr)

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

# Core traits
core_traits = c("vegetative_height", "LNC", "seed_mass", "LDMC")

# Woody species to exclude from species-level analysis
woody_species = c(
  "Fagus sylvatica", "Picea abies", "Sorbus aucuparia", "Acer pseudoplatanus",
  "Larix decidua", "Fraxinus excelsior", "Pinus sylvestris", "Corylus avellana",
  "Cornus sanguinea", "Ligustrum vulgare", "Lonicera xylosteum",
  "Rubus idaeus", "Rubus fruticosus", "Rubus caesius"
)

# Display names for each core trait
trait_display = c(
  "vegetative_height" = "Vegetative Height",
  "LNC"               = "Leaf N Content",
  "seed_mass"         = "Seed Mass",
  "LDMC"              = "Leaf Dry Matter Content"
)

# ==============================================================================
# LOAD SHARED DATA
# ==============================================================================
cat("\n=== Functional trait post-JSDM analysis ===\n")

# Detect param_tag
fold_files = list.files(
  here("Calanda_JSDM", "output", "results", "cv"),
  pattern = "^fold_.*\\.rds$"
)
param_tag = regmatches(fold_files[1],
                       regexpr("a[0-9.]+_l[0-9.]+_le[0-9.]+(?=\\.)",
                               fold_files[1], perl = TRUE))

# Site index data
data_calanda = readRDS(here("Calanda_JSDM", "output", "data_calanda_jsdm_2026-03-06.rds"))
n_sites = nrow(data_calanda$X)

# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  PART A: SPECIES-LEVEL — VP vs trait medians + kCVs                     ║
# ╚════════════════════════════════════════════════════════════════════════════╝
cat("\n========== SPECIES-LEVEL ==========\n")

# --- Load species VP ---
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
keep_names = species_cv$species[sp_keep]
keep_names = setdiff(keep_names, woody_species)
n_sp_keep = length(keep_names)
cat(sprintf("  %d / %d species kept (AUC >= 0.7, excluding %d woody)\n",
            n_sp_keep, nrow(species_cv), sum(species_cv$species[sp_keep] %in% woody_species)))

df_sp_wide = df_sp_summary %>%
  filter(species_name %in% keep_names) %>%
  filter(component %in% c("env", "spa", "codist", "r2")) %>%
  select(species_idx, species_name, component, mean_val, ci_width, weight) %>%
  pivot_wider(
    names_from = component,
    values_from = c(mean_val, ci_width, weight),
    names_glue = "{component}_{.value}"
  )

# --- Load species traits ---
sp_medians = read_csv(here("Calanda_JSDM", "output", "traits_medians_imputed.csv"),
                       show_col_types = FALSE) %>%
  rename(species_name = species_TNRS)

sp_kcv = read_csv(here("Calanda_JSDM", "output", "traits_kcv_imputed.csv"),
                   show_col_types = FALSE) %>%
  rename(species_name = species_TNRS)

# Log-transform median traits
log_traits_median = paste0("Median_", core_traits)
sp_medians = sp_medians %>%
  mutate(across(all_of(log_traits_median), ~ log(.x)))

# Merge medians + kCVs
sp_traits = sp_medians %>%
  left_join(sp_kcv, by = "species_name")

# Only keep the 5 core traits (median + kCV) + per-trait distinctiveness
sp_median_cols = paste0("Median_", core_traits)
sp_kcv_cols = paste0("kCV_", core_traits)
sp_distinct_cols = paste0("distinct_", core_traits)

# Load per-trait functional distinctiveness
sp_distinct = read_csv(
  here("Calanda_JSDM", "output", "species_functional_distinctiveness.csv"),
  show_col_types = FALSE
) %>%
  rename(species_name = species)

sp_trait_cols = c(sp_median_cols, sp_kcv_cols, sp_distinct_cols)

# Merge with VP + distinctiveness
df_sp_merged = df_sp_wide %>%
  inner_join(sp_traits, by = "species_name") %>%
  left_join(sp_distinct, by = "species_name") %>%
  drop_na(all_of(sp_trait_cols))

cat(sprintf("  Merged: %d species with traits + distinctiveness + VP\n", nrow(df_sp_merged)))

# Standardize trait predictors
df_sp_merged = df_sp_merged %>%
  mutate(across(all_of(sp_trait_cols), ~ as.numeric(scale(.x))))

# --- Weighted stepwise AIC: VP ~ all traits ---
cat("\n--- VP vs species traits (weighted stepwise AIC) ---\n")

# Model selection: set to log(n) for BIC, 2 for AIC
selection_k = 2  # AIC; change to log(nrow(df_sp_merged)) for BIC

# Linear terms only
sp_rhs = paste(sp_trait_cols, collapse = " + ")

sp_trait_models = list()
sp_trait_results = list()

for (comp in vp_components) {
  y_col = paste0(comp, "_mean_val")
  w_col = paste0(comp, "_weight")

  full_formula = as.formula(paste0(y_col, " ~ ", sp_rhs))
  full_mod = lm(full_formula, data = df_sp_merged, weights = df_sp_merged[[w_col]])
  best_mod = step(full_mod, direction = "backward", trace = 0, k = selection_k)
  sp_trait_models[[comp]] = best_mod

  s = summary(best_mod)
  cat(sprintf("\n  %s: adj.R2 = %.3f, p = %.2e\n", vp_labels[comp],
              s$adj.r.squared, pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3],
                                  lower.tail = FALSE)))
  coef_df = as.data.frame(s$coefficients) %>%
    rownames_to_column("term") %>%
    filter(term != "(Intercept)")
  print(s$coefficients)
  cat(sprintf("  Retained %d / %d terms\n", nrow(coef_df), length(sp_trait_cols)))

  sp_trait_results[[comp]] = extract_model_results(full_mod, best_mod, "functional_species", comp)
}

# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  PART B: COMMUNITY-LEVEL — VP vs unweighted trait means + variances     ║
# ╚════════════════════════════════════════════════════════════════════════════╝
cat("\n========== COMMUNITY-LEVEL ==========\n")

# --- Load site VP ---
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
keep_names = site_cv %>%
  filter(logloss <= log(2), !is.na(logloss)) %>%
  pull(site)
cat(sprintf("  %d / %d sites kept (logloss <= log(2))\n", length(keep_names), n_sites))

df_si_wide = df_si_summary %>%
  filter(site_name %in% keep_names) %>%
  filter(component %in% c("env", "spa", "codist", "r2")) %>%
  select(site_name, site_idx, component, mean_val, ci_width) %>%
  pivot_wider(
    names_from = component,
    values_from = c(mean_val, ci_width),
    names_glue = "{component}_{.value}"
  )

# Remove sites with unstable VP estimates (ci_width >= 0.1 on any component)
ci_cols = grep("_ci_width$", names(df_si_wide), value = TRUE)
unstable = apply(df_si_wide[, ci_cols], 1, function(r) any(r >= 0.1, na.rm = TRUE))
if (any(unstable)) {
  cat(sprintf("  Removing %d site(s) with ci_width >= 0.1\n", sum(unstable)))
  df_si_wide = df_si_wide[!unstable, ]
}

# --- Load community traits (unweighted) ---
comm_traits = read_csv(
  here("Calanda_JSDM", "output", "community_traits_unweighted_imputed.csv"),
  show_col_types = FALSE
)

# Only keep the 5 core traits (mean + var)
comm_mean_cols = paste0("Median_", core_traits, "_mean")
comm_var_cols = paste0("Median_", core_traits, "_var")
comm_trait_cols = c(comm_mean_cols, comm_var_cols)

# Merge site VP with community traits by site name (not positional index)
df_comm_merged = df_si_wide %>%
  rename(plot_id_releve = site_name) %>%
  inner_join(comm_traits %>% select(plot_id_releve, all_of(comm_trait_cols)),
             by = "plot_id_releve") %>%
  drop_na(all_of(comm_trait_cols))

cat(sprintf("  Merged: %d sites with community traits + VP\n", nrow(df_comm_merged)))

# Standardize community trait predictors
df_comm_merged = df_comm_merged %>%
  mutate(across(all_of(comm_trait_cols), ~ as.numeric(scale(.x))))

# --- Stepwise AIC: VP ~ all community traits (unweighted) ---
cat("\n--- VP vs community traits (stepwise AIC, unweighted) ---\n")

# Model selection: set to log(n) for BIC, 2 for AIC
selection_k_comm = 2  # AIC; change to log(nrow(df_comm_merged)) for BIC

# Linear terms only
comm_rhs = paste(comm_trait_cols, collapse = " + ")

comm_trait_models = list()
comm_trait_results = list()

for (comp in vp_components) {
  y_col = paste0(comp, "_mean_val")

  full_formula = as.formula(paste0(y_col, " ~ ", comm_rhs))
  full_mod = lm(full_formula, data = df_comm_merged)
  best_mod = step(full_mod, direction = "backward", trace = 0, k = selection_k_comm)
  comm_trait_models[[comp]] = best_mod

  s = summary(best_mod)
  cat(sprintf("\n  %s: adj.R2 = %.3f, p = %.2e\n", vp_labels[comp],
              s$adj.r.squared, pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3],
                                  lower.tail = FALSE)))
  coef_df = as.data.frame(s$coefficients) %>%
    rownames_to_column("term") %>%
    filter(term != "(Intercept)")
  print(s$coefficients)
  cat(sprintf("  Retained %d / %d terms\n", nrow(coef_df), length(comm_trait_cols)))

  comm_trait_results[[comp]] = extract_model_results(full_mod, best_mod, "functional_community", comp)
}

# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  PART C: PAIRED PARTIAL EFFECT FIGURES                                  ║
# ╚════════════════════════════════════════════════════════════════════════════╝
cat("\n========== PARTIAL EFFECT FIGURES ==========\n")

# Short labels for annotations (full labels used in legends)
annot_short = c("env" = "Environment", "codist" = "S. Associations", "spa" = "Space")

# Split traits: main figure vs other
main_traits = c("vegetative_height", "LNC", "seed_mass")
other_traits = c("LDMC")

# --- Helper to build paired figure for a set of traits ---
build_functional_figure = function(traits_subset) {
  paired_rows = list()
  first_row = TRUE
  for (tr in traits_subset) {
    med_col = paste0("Median_", tr)
    kcv_col = paste0("kCV_", tr)
    mean_col = paste0("Median_", tr, "_mean")
    var_col = paste0("Median_", tr, "_var")

    p_sp = build_paired_plot(
      evar1 = med_col, evar2 = kcv_col,
      label1 = "Median", label2 = "Plasticity",
      trait_name = trait_display[tr],
      models = sp_trait_models, df = df_sp_merged,
      predictor_cols = sp_trait_cols,
      vp_components = vp_components, vp_labels = vp_labels,
      comp_colors = comp_colors,
      evar3 = paste0("distinct_", tr), label3 = "Distinctiveness",
      annot_labels = annot_short
    )

    p_comm = build_paired_plot(
      evar1 = mean_col, evar2 = var_col,
      label1 = "Mean", label2 = "Variance",
      trait_name = trait_display[tr],
      models = comm_trait_models, df = df_comm_merged,
      predictor_cols = comm_trait_cols,
      vp_components = vp_components, vp_labels = vp_labels,
      comp_colors = comp_colors,
      annot_labels = annot_short
    )

    if (first_row) {
      p_sp = p_sp + ggtitle("A) Species") +
        theme(plot.title = element_text(size = 18, face = "bold", hjust = 0))
      p_comm = p_comm + ggtitle("B) Community") +
        theme(plot.title = element_text(size = 18, face = "bold", hjust = 0))
      first_row = FALSE
    }

    paired_rows[[length(paired_rows) + 1]] = p_sp
    paired_rows[[length(paired_rows) + 1]] = p_comm
  }

  ggarrange(
    plotlist = paired_rows, ncol = 2, nrow = length(traits_subset),
    common.legend = TRUE, legend = "bottom"
  )
}

# --- Main figure ---
p_main = build_functional_figure(main_traits)
pdf(here("Calanda_JSDM", "plot", "functional_gradient.pdf"),
    width = 12, height = length(main_traits) * 4)
print(p_main)
dev.off()
cat("  Saved functional_gradient.pdf\n")

# --- Other figure (LDMC) ---
p_other = build_functional_figure(other_traits)
pdf(here("Calanda_JSDM", "plot", "functional_gradient_other.pdf"),
    width = 12, height = length(other_traits) * 4)
print(p_other)
dev.off()
cat("  Saved functional_gradient_other.pdf\n")

# ==============================================================================
# SAVE MODEL RESULTS AS CSV
# ==============================================================================
cat("\n--- Saving model results ---\n")

all_results = c(sp_trait_results, comm_trait_results)

# Coefficients
all_coefs = bind_rows(lapply(all_results, function(r) r$coefficients))
write_csv(all_coefs, here("Calanda_JSDM", "output",
                           paste0("functional_coefficients_", param_tag, ".csv")))
cat(sprintf("  Saved functional_coefficients_%s.csv (%d rows)\n", param_tag, nrow(all_coefs)))

# AIC steps
all_aic = bind_rows(lapply(all_results, function(r) r$bic_steps))
write_csv(all_aic, here("Calanda_JSDM", "output",
                         paste0("functional_aic_steps_", param_tag, ".csv")))
cat(sprintf("  Saved functional_aic_steps_%s.csv (%d rows)\n", param_tag, nrow(all_aic)))

# Model summaries
all_summaries = bind_rows(lapply(all_results, function(r) r$model_summary))
write_csv(all_summaries, here("Calanda_JSDM", "output",
                               paste0("functional_model_summary_", param_tag, ".csv")))
cat(sprintf("  Saved functional_model_summary_%s.csv (%d rows)\n", param_tag, nrow(all_summaries)))

cat("\n=== Functional trait post-JSDM analysis complete ===\n")
