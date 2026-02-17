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
# 2. Model fitting (optional — uses results_from_Max/ by default)
# ==============================================================================
if (run_model) {
  cat("\n========== STEP 2: Model fitting ==========\n")
  source(here("Calanda_JSDM", "R", "02_model", "04_jsdm.R"))
}

# ==============================================================================
# 3. Post-model analysis
# ==============================================================================
cat("\n========== STEP 3: Analysis ==========\n")
source(here("Calanda_JSDM", "R", "03_analysis", "05_variance_partitioning.R"))
source(here("Calanda_JSDM", "R", "03_analysis", "06_community_postJSDM.R"))
source(here("Calanda_JSDM", "R", "03_analysis", "07_species_postJSDM.R"))

# ==============================================================================
# 4. Visualization
# ==============================================================================
cat("\n========== STEP 4: Visualization ==========\n")
source(here("Calanda_JSDM", "R", "04_visualization", "08_map.R"))
source(here("Calanda_JSDM", "R", "04_visualization", "09_map_rgb_results.R"))

cat("\n========== WORKFLOW COMPLETE ==========\n")
