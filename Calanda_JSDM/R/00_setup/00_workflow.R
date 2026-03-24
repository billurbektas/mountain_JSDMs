# ==============================================================================
# Script: 00_workflow.R
# Author: Billur Bektas
# Claude (Anthropic) was used to assist with code refactoring, validation, and documentation.
#
# Purpose: Coordinator script — sources all pipeline scripts in order.
#          Each script is self-contained (loads its own libraries and data).
#
# Usage:  source(here("Calanda_JSDM", "R", "00_setup", "00_workflow.R"))
#
# Toggle each step with TRUE/FALSE below.
# ==============================================================================

library(here)

# Configuration ----
run_data_prep     = FALSE   # Step 1: Data preparation (TRY traits, field traits, merge, env+veg)
run_model         = FALSE  # Step 2: Model fitting (GPU required, run on Colab)
run_experiments   = FALSE  # Step 3: Model experiment post-analysis (07.x scripts)
run_analysis      = TRUE   # Step 4: Variance partitioning + gradient analyses (08-10)
run_visualization = FALSE   # Step 5: Study area maps and RGB visualizations (11-12)
run_validation    = FALSE  # Step 6: Pipeline validation and attrition tracking

# ==============================================================================
# 1. Data preparation
# ==============================================================================
if (run_data_prep) {
  cat("\n========== STEP 1: Data preparation ==========\n")

  # 1a. Fetch TRY traits + ecological indicators + dispersal from databases
  # Inputs:  data/vegetation/2024_CAPHE_SpeDis_CleanData_20240214.csv, data/try/try_*.txt
  # Outputs: output/try_traits_individual.csv, output/indicators.csv, output/dispersal.csv
  source(here("Calanda_JSDM", "R", "01_data_prep", "01_prepare_TRY_traits.R"))

  # 1b. Process field-collected trait measurements (leaf area, LMA, LDMC, CN)
  # Inputs:  data/traits/2025_TRAITS_CleanData_20251228.csv, data/traits/leaf_area.csv
  # Outputs: output/field_traits_merged.csv, output/field_traits_clean.csv
  source(here("Calanda_JSDM", "R", "01_data_prep", "02_prepare_field_trait_data.R"))

  # 1c. Merge TRY + field traits, compute species medians and kCV, impute missing values
  # Inputs:  output/try_traits_individual.csv, output/field_traits_clean.csv, output/indicators.csv, output/dispersal.csv
  # Outputs: output/traits_medians_imputed.csv, output/traits_kcv_imputed.csv, output/all_traits_individual.csv
  source(here("Calanda_JSDM", "R", "01_data_prep", "03_merge_and_assess_traits.R"))

  # 1d. Build environmental (X) and species (Y) matrices, community traits, and distinctiveness
  # Inputs:  data/vegetation/*, data/ecostress/*, data/modis/*, data/wekeo/*, output/traits_medians_imputed.csv
  # Outputs: output/data_calanda_jsdm_<date>.rds, output/veg_clim.csv, output/community_traits_*.csv, output/species_functional_distinctiveness.csv
  source(here("Calanda_JSDM", "R", "01_data_prep", "04_prepare_climate_vegetation_data.R"))
}

# ==============================================================================
# 2. Model fitting (local CPU test; real runs on Google Colab GPU)
# ==============================================================================
if (run_model) {
  cat("\n========== STEP 2: Model fitting (CPU test) ==========\n")

  # Local CPU test of the sjSDM model structure (low iterations)
  # Inputs:  output/data_calanda_jsdm_<date>.rds
  # Outputs: none (test only)
  source(here("Calanda_JSDM", "R", "02_model", "05_test_cpu_jsdm.R"))
}
# Full model fitting + 6 experiments run via Google Colab notebook:
#   R/02_model/06_sjsdm_experiments_colab.ipynb
# See README for details on Colab execution and the anova bug fix.

# ==============================================================================
# 3. Model experiments post-analysis (07.x scripts)
# ==============================================================================
if (run_experiments) {
  cat("\n========== STEP 3: Model experiments post-analysis ==========\n")

  if (dir.exists(here("Calanda_JSDM", "output", "results", "runs"))) {
    source(here("Calanda_JSDM", "R", "03_analysis", "07.1_post_exp1_grid.R"))
  }
  if (dir.exists(here("Calanda_JSDM", "output", "results", "decoupled"))) {
    source(here("Calanda_JSDM", "R", "03_analysis", "07.2_post_exp2_decoupled.R"))
  }
  if (dir.exists(here("Calanda_JSDM", "output", "results", "anova_saturation"))) {
    source(here("Calanda_JSDM", "R", "03_analysis", "07.3_post_exp3_anova_saturation.R"))
    source(here("Calanda_JSDM", "R", "03_analysis", "07.4_post_exp3_species_stability.R"))
  }
  if (dir.exists(here("Calanda_JSDM", "output", "results", "fit_saturation"))) {
    source(here("Calanda_JSDM", "R", "03_analysis", "07.5_post_exp4_fit_saturation.R"))
  }
  if (dir.exists(here("Calanda_JSDM", "output", "results", "dropone"))) {
    source(here("Calanda_JSDM", "R", "03_analysis", "07.6_post_exp5_dropone.R"))
  }
  if (dir.exists(here("Calanda_JSDM", "output", "results", "cv"))) {
    source(here("Calanda_JSDM", "R", "03_analysis", "07.7_post_exp6_cv.R"))
  }

  # Consolidated experiment figures from all experiments
  # Inputs:  output/results/exp*_*.csv, output/results/cv/fold_*.rds
  # Outputs: plot/exp_*.pdf
  source(here("Calanda_JSDM", "R", "03_analysis", "07_modeling_experiments.R"))
}

# ==============================================================================
# 4. Post-model analysis
# ==============================================================================
if (run_analysis) {
  cat("\n========== STEP 4: Post-model analysis ==========\n")

  # Variance partitioning from 10-fold CV: Venn diagram, species/site VP summaries
  # Inputs:  output/results/cv/fold_*.rds, output/results/exp6_*_cv_*.csv, output/data_calanda_jsdm_*.rds
  # Outputs: output/results/vp_species_summary_*.csv, output/results/vp_sites_summary_*.csv, plot/vp_combined.pdf
  source(here("Calanda_JSDM", "R", "03_analysis", "08_variance_partitioning.R"))

  # Environmental gradient analysis: VP vs env predictors and species betas
  # Inputs:  output/results/vp_*_summary_*.csv, output/results/exp6_*_cv_*.csv, output/data_calanda_jsdm_*.rds, output/final_model_se_*.rds
  # Outputs: plot/env_gradient_climate.pdf, plot/env_gradient_other.pdf
  source(here("Calanda_JSDM", "R", "03_analysis", "09_environmental_gradient_analysis.R"))

  # Functional trait analysis: VP vs species traits and community traits
  # Inputs:  output/results/vp_*_summary_*.csv, output/results/exp6_*_cv_*.csv, output/traits_*.csv, output/community_traits_unweighted_imputed.csv, output/species_functional_distinctiveness.csv
  # Outputs: plot/functional_gradient.pdf
  source(here("Calanda_JSDM", "R", "03_analysis", "10_functional_post_analysis.R"))
}

# ==============================================================================
# 5. Visualization
# ==============================================================================
if (run_visualization) {
  cat("\n========== STEP 5: Visualization ==========\n")

  # Study area map with vegetation survey points
  # Inputs:  output/starter_data_*.RData, output/calanda_mask.shp
  # Outputs: plot/map_calanda.pdf
  source(here("Calanda_JSDM", "R", "04_visualization", "11_map.R"))

  # 3D RGB rayshader map of variance components on DEM
  # Inputs:  output/starter_data_*.RData, results/res_sjsdm_calanda.rds, output/calanda_mask.shp, output/metrics/dem.tif
  # Outputs: plot/calanda_rgb_3d.png, plot/calanda_rgb_2d_map.pdf, plot/rgb_legend.pdf
  source(here("Calanda_JSDM", "R", "04_visualization", "12_map_rgb_results.R"))
}

# ==============================================================================
# 6. Validation (optional)
# ==============================================================================
if (run_validation) {
  cat("\n========== STEP 6: Validation ==========\n")

  # Pipeline integrity checks: data alignment, join correctness, regression assumptions
  # Inputs:  all output files
  # Outputs: console diagnostics only
  source(here("Calanda_JSDM", "R", "05_validation", "validation_checks.R"))

  # Species and site attrition tracking from raw data through all filters
  # Inputs:  data/vegetation/*, all output files
  # Outputs: output/species_attrition.csv, output/site_attrition.csv, output/attrition_summary.csv
  source(here("Calanda_JSDM", "R", "05_validation", "pipeline_attrition.R"))
}

cat("\n========== WORKFLOW COMPLETE ==========\n")
