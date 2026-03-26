# Calanda JSDM — Workflow Details

See the [main README](../../README.md) for an overview. This document contains detailed inputs and outputs for each script.

---

## Step 1: Data preparation

### 01_prepare_TRY_traits.R
Fetches trait data from TRY database and ecological indicators/dispersal from floraveg.eu.

| Inputs | Outputs |
|--------|---------|
| `data/vegetation/2024_CAPHE_SpeDis_CleanData_20240214.csv` | `output/try_traits_individual.csv` |
| `data/try/try_*.txt` | `output/indicators.csv` |
| | `output/dispersal.csv` |
| | `output/try_trait_units.csv` |

### 02_prepare_field_trait_data.R
Processes field-collected trait measurements (leaf area, LMA, LDMC, CN) with quality validation.

| Inputs | Outputs |
|--------|---------|
| `data/traits/2025_TRAITS_CleanData_20251228.csv` | `output/field_traits_merged.csv` |
| `data/traits/leaf_area.csv` | `output/field_traits_clean.csv` |
| | `output/field_trait_units.csv` |

### 03_merge_and_assess_traits.R
Merges TRY + field traits, computes species medians and kCV, imputes missing values via random forest.

| Inputs | Outputs |
|--------|---------|
| `output/try_traits_individual.csv` | `output/all_traits_individual.csv` |
| `output/field_traits_clean.csv` | `output/traits_medians_imputed.csv` |
| `output/indicators.csv` | `output/traits_kcv_imputed.csv` |
| `output/dispersal.csv` | `output/traits_medians_imputation_flags.csv` |
| | `output/traits_kcv_imputation_flags.csv` |
| | `output/species_imputation_summary.csv` |
| | `plot/pca_species_means.pdf` |
| | `plot/pca_species_kcv.pdf` |
| | `plot/trait_distributions_raw_vs_log.pdf` |
| | `plot/imputation_distributions.pdf` |

### 04_prepare_climate_vegetation_data.R
Builds environmental (X) and species (Y) matrices from remote sensing and vegetation surveys; computes community traits and functional distinctiveness.

| Inputs | Outputs |
|--------|---------|
| `data/vegetation/2024_CAPHE_SpeDis_CleanData_20240214.csv` | `output/data_calanda_jsdm_<date>.rds` |
| `data/vegetation/veg.coord.csv` | `output/starter_data_<date>.RData` |
| `data/ecostress/*.csv` | `output/veg_clim.csv` |
| `data/modis/*.csv` | `output/veg_tree.csv` |
| `data/wekeo/` (Copernicus snow) | `output/community_traits_imputed.csv` |
| `data/mask/*.shp` | `output/community_traits_original.csv` |
| `output/traits_medians_imputed.csv` | `output/community_traits_unweighted_imputed.csv` |
| `output/traits_medians_imputation_flags.csv` | `output/community_traits_unweighted_original.csv` |
| `output/traits_kcv_imputation_flags.csv` | `output/species_functional_distinctiveness.csv` |
| | `output/species_imputation_summary.csv` |
| | `output/cwm_comparison.csv` |
| | `plot/cwm_imputed_vs_original.pdf` |
| | `plot/climate_calanda_surveys.pdf` |

---

## Step 2: Model fitting

### 05_test_cpu_jsdm.R
Local CPU test of the sjSDM model structure (low iterations, for debugging). No outputs.

| Inputs | Outputs |
|--------|---------|
| `output/data_calanda_jsdm_<date>.rds` | -- |

### 06_sjsdm_experiments_colab.ipynb
Full model fitting + 6 experiments on Google Colab GPU (A100). Includes runtime anova bug fix.

| Inputs | Outputs |
|--------|---------|
| `output/data_calanda_jsdm_<date>.rds` (uploaded) | `output/results/cv/fold_*.rds` |
| | `output/results/exp6_species_cv_*.csv` |
| | `output/results/exp6_sites_cv_*.csv` |
| | `output/final_model_se_*.rds` |
| | `output/results/runs/`, `decoupled/`, `anova_saturation/`, `fit_saturation/`, `dropone/` |

---

## Step 3: Model experiment post-analysis

### 07.1–07.7_post_exp*.R
Individual post-analysis for each of the 6 experiments.

| Inputs | Outputs |
|--------|---------|
| `output/results/runs/`, `decoupled/`, etc. | `output/results/exp*_*.csv` |
| | `plot/detailed_modeling_experiment_figures/*.pdf` |

### 07_modeling_experiments.R
Consolidated publication-ready figures from all experiments.

| Inputs | Outputs |
|--------|---------|
| `output/results/exp*_*.csv` | `plot/fig_exp1_anova_R2.pdf` |
| `output/results/cv/fold_*.rds` | `plot/fig_exp2_anova_R2.pdf` |
| `output/results/dropone/*.rds` | `plot/fig_exp3_correlations.pdf` |
| | `plot/fig_exp3_violins.pdf` |
| | `plot/fig_exp5_dropone.pdf` |
| | `plot/fig_exp6_prediction_violins.pdf` |
| | `plot/fig_exp6_species_ranked.pdf` |
| | `plot/fig_exp6_sites_ranked.pdf` |
| | `plot/fig_merged_convergence.pdf` |
| | `plot/fig_merged_timing.pdf` |

---

## Step 4: Post-model analysis

### 08_variance_partitioning.R
Aggregates species and site VP across 10-fold CV, computes mean R² with 95% CIs, produces Venn diagram and scatter plots.

| Inputs | Outputs |
|--------|---------|
| `output/results/cv/fold_*.rds` | `output/results/vp_species_summary_*.csv` |
| `output/results/exp6_species_cv_*.csv` | `output/results/vp_sites_summary_*.csv` |
| `output/results/exp6_sites_cv_*.csv` | `output/results/vp_anova_summary_*.csv` |
| `output/data_calanda_jsdm_*.rds` | `plot/vp_combined.pdf` |

### 09_environmental_gradient_analysis.R
AIC stepwise regressions of VP components against 12 environmental predictors (community, unweighted) and species betas (species, weighted).

| Inputs | Outputs |
|--------|---------|
| `output/results/vp_sites_summary_*.csv` | `output/env_gradient_coefficients_*.csv` |
| `output/results/vp_species_summary_*.csv` | `output/env_gradient_aic_steps_*.csv` |
| `output/results/exp6_sites_cv_*.csv` | `output/env_gradient_model_summary_*.csv` |
| `output/results/exp6_species_cv_*.csv` | `plot/env_gradient_climate.pdf` |
| `output/data_calanda_jsdm_*.rds` | `plot/env_gradient_cover.pdf` |
| `output/final_model_se_*.rds` | `plot/env_gradient_topography.pdf` |
| | `plot/env_gradient_indices.pdf` |

### 10_functional_post_analysis.R
AIC stepwise regressions of VP against species trait medians, kCV, and distinctiveness (weighted); community trait means and variances (unweighted).

| Inputs | Outputs |
|--------|---------|
| `output/results/vp_species_summary_*.csv` | `output/functional_coefficients_*.csv` |
| `output/results/vp_sites_summary_*.csv` | `output/functional_aic_steps_*.csv` |
| `output/results/exp6_species_cv_*.csv` | `output/functional_model_summary_*.csv` |
| `output/results/exp6_sites_cv_*.csv` | `plot/functional_gradient.pdf` |
| `output/data_calanda_jsdm_*.rds` | `plot/functional_gradient_other.pdf` |
| `output/traits_medians_imputed.csv` | |
| `output/traits_kcv_imputed.csv` | |
| `output/community_traits_unweighted_imputed.csv` | |
| `output/species_functional_distinctiveness.csv` | |

---

## Step 5: Visualization

### 11_map.R
Study area map with greyscale DEM relief, communities colored by RGB climate triangle (summer temp = R, ET = G, |FDD| = B), open/closed habitat shapes, species richness sizes.

| Inputs | Outputs |
|--------|---------|
| `output/veg_clim.csv` | `plot/map_calanda.pdf` |
| `output/data_calanda_jsdm_*.rds` | |
| `output/community_traits_unweighted_imputed.csv` | |
| `output/calanda_mask.shp` | |
| `output/metrics/dem.tif` | |
| `data/vegetation/2024_CAPHE_SpeDis_CleanData_20240214.csv` | |

---

## Step 6: Validation (optional)

### validation_checks.R
Pipeline integrity checks: data alignment, join correctness, VIF, residual normality, heteroscedasticity.

| Inputs | Outputs |
|--------|---------|
| All output files | Console diagnostics only |

### pipeline_attrition.R
Tracks species and site attrition from raw data through all filtering stages with reasons.

| Inputs | Outputs |
|--------|---------|
| `data/vegetation/2024_CAPHE_SpeDis_CleanData_20240214.csv` | `output/species_attrition.csv` |
| All output files | `output/site_attrition.csv` |
| | `output/attrition_summary.csv` |

### vp_diagnostic_figures.R
Diagnostic figures comparing proportional vs discard variance allocation and raw shared fraction distributions.

| Inputs | Outputs |
|--------|---------|
| `output/results/cv/fold_*.rds` | `plot/vp_diagnostic_discard_vs_proportional.pdf` |
| | `plot/vp_diagnostic_discard_vs_proportional_allfolds.pdf` |
| | `plot/vp_diagnostic_raw_shared_fractions.pdf` |
| | `plot/vp_diagnostic_raw_shared_fractions_allfolds.pdf` |

---

## Key methodological notes

### Variance partitioning levels
The Venn diagram (landscape) and violins (per-unit) use different aggregation and do not match numerically. See `docs/variance_partitioning_notes.md` for a detailed explanation.

### Site-level regressions
Community-level regressions in scripts 09 and 10 are **unweighted** (sites with ci_width >= 0.1 are excluded instead). Species-level regressions remain weighted by 1/ci_width².

### Name-based joins
All species and site matching uses name-based joins (not positional indices) to prevent silent data misalignment across folds.

---

*Last updated: 2026-03-25*
