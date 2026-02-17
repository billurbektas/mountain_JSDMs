# ==============================================================================
# Script: 02_prepare_field_trait_data.R
# Purpose: Process field-collected trait measurements into clean individual-level
#          data for merging with TRY traits in step 03.
#
# Description:
#   Loads pre-cleaned trait data (with CN already integrated), merges with
#   scanned leaf area measurements, calculates derived traits (LMA, LDMC),
#   validates data quality (dry vs fresh weight, C content outliers), and
#   creates a clean dataset with problematic values set to NA.
#
# Input files:
#   - data/traits/2025_TRAITS_CleanData_20251228.csv (pre-cleaned trait data)
#   - data/traits/leaf_area.csv (scanned leaf areas)
#
# Output files:
#   - output/field_traits_merged.csv (all traits merged, with error flags)
#   - output/field_traits_clean.csv (clean version, errors set to NA)
#
# ==============================================================================

# Load required libraries ----
library(tidyverse)
library(here)

### Merge all field measured traits data ----

# Read traits_clean data
traits_clean = read_csv(here("Calanda_JSDM", "data", "traits", "2025_TRAITS_CleanData_20251228.csv"))

# Apply additional cleaning steps based on the original cleaning script ----
cat("\n--- Additional data validation ---\n")

## 1. Validate fresh vs dry weights (dry should always be < fresh)
traits_clean = traits_clean %>%
  mutate(
    biomass_weight_diff = biomass_fresh_weight - biomass_dry_weight,
    flowering_unit_weight_diff = flowering_unit_fresh_weight - flowering_unit_dry_weight,
    leaves_weight_diff = leaves_fresh_weight - leaves_dry_weight
  )

# Flag samples with negative differences (measurement errors)
weight_error_ids = traits_clean %>%
  filter(
    (biomass_weight_diff < 0 & !is.na(biomass_weight_diff)) |
    (flowering_unit_weight_diff < 0 & !is.na(flowering_unit_weight_diff)) |
    (leaves_weight_diff < 0 & !is.na(leaves_weight_diff))
  ) %>%
  pull(plant_id)

if(length(weight_error_ids) > 0) {
  cat("Found", length(weight_error_ids), "samples with dry > fresh weight (measurement errors)\n")
  traits_clean = traits_clean %>%
    mutate(weight_measurement_error = plant_id %in% weight_error_ids)
} else {
  cat("No weight measurement errors found\n")
  traits_clean = traits_clean %>%
    mutate(weight_measurement_error = FALSE)
}
# Identify outlier CN values 
outlier_c_ids =
  traits_clean %>%
  filter(C_content_corr > 70 | C_content_corr < 20) %>% # Above/below these numbers probably measurement error.
  pull(plant_id)

if(length(outlier_c_ids) > 0) {
  cat("Found", length(outlier_c_ids), "samples with outlier C values (measurement errors)\n")
  traits_clean = traits_clean %>%
    mutate(c_measurement_error = plant_id %in% outlier_c_ids)
} else {
  cat("No C measurement errors found\n")
  traits_clean = traits_clean %>%
    mutate(c_measurement_error = FALSE)
}

# Remove temporary diff columns
traits_clean = traits_clean %>%
  select(-biomass_weight_diff, -flowering_unit_weight_diff, -leaves_weight_diff)

# Read leaf_area data
leaf_area = read_csv(here("Calanda_JSDM", "data", "traits", "leaf_area.csv"))

# Check dimensions
cat("\n--- Merging datasets ---\n")
cat("traits_clean dimensions:", nrow(traits_clean), "rows x", ncol(traits_clean), "columns\n")
cat("leaf_area dimensions:", nrow(leaf_area), "rows x", ncol(leaf_area), "columns\n")

# Merge all three datasets by plant_id
all_traits = traits_clean %>%
  left_join(leaf_area %>%
              select(plant_id, n_scans, n_leaves, sum_leaf_cm2, mean_leaf_cm2, LMA),
            by = "plant_id")

# Calculate LDMC (Leaf Dry Matter Content) ----
cat("\n--- Calculating LDMC ---\n")

# LDMC = leaf's oven-dry mass / water-saturated fresh mass (mg/g)
# Units: mg/g (so multiply by 1000 to convert from g/g to mg/g)
all_traits = all_traits %>%
  mutate(
    LDMC = ifelse(!is.na(leaves_dry_weight) & !is.na(leaves_fresh_weight) & leaves_fresh_weight > 0,
                  (leaves_dry_weight / leaves_fresh_weight) * 1000,  # Convert to mg/g
                  NA_real_)
  )

# Select final trait variables for analysis ----
cat("\n--- Selecting final trait variables ---\n")

final_traits = all_traits %>%
  select(
    # Identifiers
    plant_id,
    plant_species,
    individual_nr,
    site,
    lat,
    lon,
    date_field,
    collector,

    # Plant height measurements
    vegetative_height,
    reproductive_height,
    vegetative_height_stretched,
    reproductive_height_stretched,

    # Biomass measurements
    leaves_fresh_weight,
    leaves_dry_weight,
    biomass_fresh_weight,
    biomass_dry_weight,

    # Leaf traits
    LMA,              # Leaf Mass per Area (kg/m²)
    LDMC,             # Leaf Dry Matter Content (mg/g)
    C_content_corr,              # Leaf Nitrogen Content (corrected %)
    N_content_corr,              # Leaf Carbon Content (corrected %)
    n_leaves,         # Number of leaves scanned
    sum_leaf_cm2,     # Total leaf area

    # Quality flags
    weight_measurement_error,    # Weight measurement error flag
    c_measurement_error          # Carbon content measurement error flag
  )

cat("\nQuality flags:\n")
cat("  - Weight measurement errors:", sum(final_traits$weight_measurement_error, na.rm=TRUE), "\n")
cat("  - Carbon content measurement outliers:", sum(final_traits$c_measurement_error, na.rm=TRUE), "\n")

# Create cleaned dataset (remove problematic values, not entire rows) ----
cat("\n--- Creating error-free dataset ---\n")

final_traits_clean = final_traits %>%
  mutate(
    # Set weight-dependent traits to NA if there's a weight measurement error
    LMA = ifelse(weight_measurement_error, NA_real_, LMA),
    LDMC = ifelse(weight_measurement_error, NA_real_, LDMC),
    leaves_fresh_weight = ifelse(weight_measurement_error, NA_real_, leaves_fresh_weight),
    leaves_dry_weight = ifelse(weight_measurement_error, NA_real_, leaves_dry_weight),

    # Set nutrient traits to NA if they are flagged as outliers (keep=FALSE)
    C_content_corr = ifelse(c_measurement_error, NA_real_, C_content_corr)
  )

cat("Removed problematic values:\n")

# Save all datasets
write_csv(all_traits, here("Calanda_JSDM", "output", "field_traits_merged.csv"))
write_csv(final_traits_clean, here("Calanda_JSDM", "output", "field_traits_clean.csv"))

# Save unit reference for downstream harmonization ----
field_units = tibble(
  trait = c("vegetative_height", "reproductive_height", "LMA", "LDMC",
            "N_content_corr", "C_content_corr", "sum_leaf_cm2"),
  unit = c("mm", "mm", "kg/m2", "mg/g", "% dry mass", "% dry mass", "cm2")
)
write_csv(field_units, here("Calanda_JSDM", "output", "field_trait_units.csv"))
cat("Saved field trait unit reference to output/field_trait_units.csv\n")
