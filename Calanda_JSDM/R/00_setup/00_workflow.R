# ==============================================================================
# Script: 00_workflow.R
# Purpose: Coordinator script — sources all pipeline scripts in order.
#          Each script is self-contained (loads its own libraries and data).
#
# Usage:  source(here("Calanda_JSDM", "R", "00_setup", "00_workflow.R"))
#
# Note:   Step 2 (model fitting) requires GPU and is slow. By default it is
#         skipped because authoritative results live in results_from_Max/.
#         Set run_model = TRUE below to re-fit.
# ==============================================================================

library(here)

# Configuration ----
run_model = FALSE   # Set TRUE to re-fit the sjSDM model (GPU required)

# ==============================================================================
# 1. Data preparation
# ==============================================================================
cat("\n========== STEP 1: Data preparation ==========\n")

# 1a. Fetch TRY traits + indicators + dispersal from databases
source(here("Calanda_JSDM", "R", "01_data_prep", "01_prepare_TRY_traits.R"))

# 1b. Process field-collected trait measurements
source(here("Calanda_JSDM", "R", "01_data_prep", "02_prepare_field_trait_data.R"))

# 1c. Merge TRY + field traits, species summaries, impute, and assess coverage
#     (produces traits.csv needed by step 1d; coverage assessment runs only if
#      starter_data from a previous run of step 1d is available)
source(here("Calanda_JSDM", "R", "01_data_prep", "03_merge_and_assess_traits.R"))

# 1d. Prepare environmental data (uses traits.csv from step 1c)
source(here("Calanda_JSDM", "R", "01_data_prep", "04_prepare_climate_vegetation_data.R"))

# ==============================================================================
# 2. Model fitting (Test locally on CPU; real runs on Google Colab GPU)
# ==============================================================================
if (run_model) {
  cat("\n========== STEP 2: Model fitting ==========\n")
  source(here("Calanda_JSDM", "R", "02_model", "05_test_cpu_jsdm.R"))
}
# Run the model experiments via Google Colab notebook: R/02_model/06_sjsdm_experiments_colab.ipynb

# ==============================================================================
# 3. Model experiments post-analysis
# ==============================================================================
if (dir.exists(here("Calanda_JSDM", "output", "results"))) {
  cat("\n========== STEP 3: Model experiments post analysis ==========\n")
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
  # Consolidated final experiment figures
  source(here("Calanda_JSDM", "R", "03_analysis", "07_modeling_experiments.R"))
}

# ==============================================================================
# 4. Post-model analysis
# ==============================================================================
cat("\n========== STEP 4: Post-model analysis ==========\n")
source(here("Calanda_JSDM", "R", "03_analysis", "08_variance_partitioning.R"))
source(here("Calanda_JSDM", "R", "03_analysis", "09_environmental_gradient_analysis.R"))

# ==============================================================================
# 5. Visualization
# ==============================================================================
cat("\n========== STEP 5: Visualization ==========\n")
source(here("Calanda_JSDM", "R", "04_visualization", "11_map.R"))
source(here("Calanda_JSDM", "R", "04_visualization", "12_map_rgb_results.R"))

cat("\n========== WORKFLOW COMPLETE ==========\n")
