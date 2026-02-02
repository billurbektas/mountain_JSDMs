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
source(here("Calanda_JSDM", "R", "01_data_prep", "01_prepare_data.R"))
source(here("Calanda_JSDM", "R", "01_data_prep", "01_prepare_trait_data.R"))
source(here("Calanda_JSDM", "R", "01_data_prep", "assess_trait_coverage.R"))

# ==============================================================================
# 2. Model fitting (optional — uses results_from_Max/ by default)
# ==============================================================================
if (run_model) {
  cat("\n========== STEP 2: Model fitting ==========\n")
  source(here("Calanda_JSDM", "R", "02_model", "02_jsdm.R"))
}

# ==============================================================================
# 3. Post-model analysis
# ==============================================================================
cat("\n========== STEP 3: Analysis ==========\n")
source(here("Calanda_JSDM", "R", "03_analysis", "03_variance_partitioning.R"))
source(here("Calanda_JSDM", "R", "03_analysis", "05_community_postJSDM.R"))
source(here("Calanda_JSDM", "R", "03_analysis", "06_species_postJSDM.R"))

# ==============================================================================
# 4. Visualization
# ==============================================================================
cat("\n========== STEP 4: Visualization ==========\n")
source(here("Calanda_JSDM", "R", "04_visualization", "map.R"))
source(here("Calanda_JSDM", "R", "04_visualization", "map_rgb_results.R"))

cat("\n========== WORKFLOW COMPLETE ==========\n")
