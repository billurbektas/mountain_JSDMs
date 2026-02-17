# ==============================================================================
# Script: 03_merge_traits.R
# Purpose: Merge TRY and field trait data, calculate species summaries
#
# Description:
#   Combines individual-level trait observations from:
#   1. TRY database (output/try_traits_individual.csv from 01b)
#   2. Field measurements (output/final_traits_clean.csv from 02)
#
#   Then calculates species-level means and variances from the combined data,
#   merges with ecological indicators and dispersal traits, and imputes
#   missing values.
#
# Input files:
#   - output/try_traits_individual.csv (from 01b_fetch_try_traits.R)
#   - output/final_traits_clean.csv (from 02_prepare_trait_data.R)
#   - output/indicators.csv (from 01b_fetch_try_traits.R)
#   - output/dispersal.csv (from 01b_fetch_try_traits.R)
#
# Output files:
#   - output/all_traits_individual.csv (combined individual-level data)
#   - output/traits_raw.csv (species summaries before imputation)
#   - output/traits.csv (final imputed traits for JSDM pipeline)
#
# Requires: R/00_setup/functions_calanda.R (for impute_functional_traits)
# ==============================================================================

library(tidyverse)
library(here)
library(gt)

source(here("Calanda_JSDM", "R", "00_setup", "functions_calanda.R"))

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================
cat("\n=== Loading trait data ===\n")

# Load TRY individual-level data
try_individual = read_csv(
  here("Calanda_JSDM", "output", "try_traits_individual.csv"),
  show_col_types = FALSE
)
cat("TRY data:", nrow(try_individual), "observations\n")

# Load field individual-level data
field_individual = read_csv(
  here("Calanda_JSDM", "output", "final_traits_clean.csv"),
  show_col_types = FALSE
)
cat("Field data:", nrow(field_individual), "observations\n")

# Load indicators and dispersal
indicators = read_csv(
  here("Calanda_JSDM", "output", "indicators.csv"),
  show_col_types = FALSE
)

dispersal = read_csv(
  here("Calanda_JSDM", "output", "dispersal.csv"),
  show_col_types = FALSE
)

# ==============================================================================
# 2. HARMONIZE TRAIT NAMES AND RESHAPE TO LONG FORMAT
# ==============================================================================
cat("\n=== Harmonizing trait names ===\n")

# TRY data is already in long format with TraitName and Value columns
# Standardize trait names
try_long = try_individual %>%
  mutate(
    trait = case_when(
      TraitName == "N_percent" ~ "LNC",
      TraitName == "C_percent" ~ "LCC",
      TraitName == "SLA" ~ "SLA",
      TraitName == "LA" ~ "LA",
      TraitName == "LDMC" ~ "LDMC",
      TraitName == "vegetative_height" ~ "vegetative_height",
      TraitName == "reproductive_height" ~ "reproductive_height",
      TraitName == "seed_mass" ~ "seed_mass",
      TraitName == "flowering_phenology" ~ "flowering_phenology",
      TRUE ~ TraitName
    ),
    source = "TRY"
  ) %>%
  select(species_TNRS, trait, value = Value, source)

cat("TRY traits standardized:", n_distinct(try_long$trait), "traits\n")

# Field data needs to be reshaped from wide to long
# Field traits: vegetative_height, reproductive_height, LMA, LDMC, LNC, LCC
field_long = field_individual %>%
  select(
    plant_species,
    vegetative_height,
    reproductive_height,
    LMA,
    LDMC,
    LNC,
    LCC
  ) %>%
  pivot_longer(
    cols = -plant_species,
    names_to = "trait",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%
  rename(species_TNRS = plant_species) %>%
  mutate(source = "Field")

cat("Field traits standardized:", n_distinct(field_long$trait), "traits\n")

# ==============================================================================
# 3. COMBINE TRY AND FIELD DATA
# ==============================================================================
cat("\n=== Combining TRY and field data ===\n")

# Combine both sources
all_individual = bind_rows(try_long, field_long)

cat("Combined dataset:", nrow(all_individual), "observations\n")
cat("Species:", n_distinct(all_individual$species_TNRS), "\n")
cat("Traits:", n_distinct(all_individual$trait), "\n")

# Summary by source
cat("\nObservations by source:\n")
all_individual %>%
  group_by(source) %>%
  summarize(
    n_obs = n(),
    n_species = n_distinct(species_TNRS),
    .groups = "drop"
  ) %>%
  print()

# Summary by trait and source
cat("\nObservations by trait and source:\n")
all_individual %>%
  group_by(trait, source) %>%
  summarize(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = source, values_from = n, values_fill = 0) %>%
  mutate(total = TRY + Field) %>%
  arrange(desc(total)) %>%
  print()

# Save combined individual-level data
write_csv(all_individual, here("Calanda_JSDM", "output", "all_traits_individual.csv"))
cat("\nSaved combined individual data to output/all_traits_individual.csv\n")

# ==============================================================================
# 4. CALCULATE SPECIES-LEVEL SUMMARIES
# ==============================================================================
cat("\n=== Calculating species-level trait summaries ===\n")

species_traits = all_individual %>%
  group_by(species_TNRS, trait) %>%
  summarize(
    Mean = mean(value, na.rm = TRUE),
    Var = var(value, na.rm = TRUE),
    n_obs = n(),
    n_sources = n_distinct(source),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = trait,
    values_from = c(Mean, Var, n_obs, n_sources),
    names_glue = "{.value}_{trait}"
  )

cat("Species summaries calculated for", nrow(species_traits), "species\n")

# ==============================================================================
# 5. MERGE WITH INDICATORS AND DISPERSAL
# ==============================================================================
cat("\n=== Merging with indicators and dispersal ===\n")

# Merge all trait sources
traits = species_traits %>%
  left_join(dispersal %>% select(species_TNRS, dispersal), by = "species_TNRS") %>%
  left_join(indicators %>% select(species_TNRS, Light, Moisture, Nutrients, disturbance),
            by = "species_TNRS")

cat("Merged dataset:", nrow(traits), "species\n")

# Select and order columns for output
traits_ordered = traits %>%
  select(
    species_TNRS,
    # Seed and dispersal
    Mean_seed_mass,
    Var_seed_mass,
    dispersal,
    # Leaf area traits
    Mean_LA,
    Var_LA,
    Mean_SLA,
    Var_SLA,
    Mean_LMA,
    Var_LMA,
    Mean_LDMC,
    Var_LDMC,
    # Plant size
    Mean_vegetative_height,
    Var_vegetative_height,
    Mean_reproductive_height,
    Var_reproductive_height,
    # Leaf chemistry
    Mean_LNC,
    Var_LNC,
    Mean_LCC,
    Var_LCC,
    # Phenology
    Mean_flowering_phenology,
    Var_flowering_phenology,
    # Ecological indicators
    Light,
    Moisture,
    Nutrients,
    disturbance,
    # Sample sizes (optional, for QC)
    starts_with("n_obs_"),
    starts_with("n_sources_")
  )

# Save raw (un-imputed) traits
write_csv(traits_ordered, here("Calanda_JSDM", "output", "traits_raw.csv"))
cat("Saved raw traits to output/traits_raw.csv\n")

# ==============================================================================
# 6. IMPUTE MISSING TRAIT VALUES
# ==============================================================================
cat("\n=== Imputing missing trait values ===\n")

# Select only the columns needed for imputation (means and categorical)
traits_for_imputation = traits_ordered %>%
  select(
    species_TNRS,
    Mean_seed_mass,
    dispersal,
    Mean_LA,
    Mean_SLA,
    Mean_LMA,
    Mean_LDMC,
    Mean_vegetative_height,
    Mean_reproductive_height,
    Mean_LNC,
    Mean_LCC,
    Mean_flowering_phenology,
    Light,
    Moisture,
    Nutrients,
    disturbance,
    # Keep variances for output
    starts_with("Var_")
  )

# Variables to impute
variables_to_impute = c(
  "Mean_seed_mass",
  "dispersal",
  "Mean_LA",
  "Mean_SLA",
  "Mean_LMA",
  "Mean_LDMC",
  "Mean_vegetative_height",
  "Mean_reproductive_height",
  "Mean_LNC",
  "Mean_LCC",
  "Mean_flowering_phenology",
  "Light",
  "Moisture",
  "Nutrients",
  "disturbance"
)

# Show missing data before imputation
cat("\nMissing data before imputation:\n")
traits_for_imputation %>%
  summarize(across(everything(), ~ round(100 * mean(is.na(.)), 1))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "percent_missing") %>%
  filter(percent_missing > 0) %>%
  arrange(desc(percent_missing)) %>%
  print(n = 25)

# Impute missing values using random forest
imp_traits = impute_functional_traits(
  traits_for_imputation,
  variables_to_impute = variables_to_impute,
  m = 30,
  maxiter = 100,
  num_trees = 500,
  seed = 123,
  validation_fraction = 0.3
)

# Display imputation performance
cat("\nImputation performance:\n")
traits_for_imputation %>%
  summarize(across(everything(), ~ mean(is.na(.)) * 100)) %>%
  pivot_longer(cols = everything(),
               names_to = "variable",
               values_to = "percent_missing") %>%
  arrange(desc(percent_missing)) %>%
  right_join(imp_traits$performance, by = "variable") %>%
  mutate(
    percent_missing = round(percent_missing),
    r_squared = round(r_squared, 2)
  ) %>%
  select(variable, percent_missing, r_squared) %>%
  gt()

# Extract imputed data and clean column names
traits_imputed = imp_traits$imputed_data %>%
  select(species_TNRS, starts_with("Var_"), contains("_final")) %>%
  rename_with(~ str_remove(., "_final$"), ends_with("_final"))

# Save final imputed traits
write_csv(traits_imputed, here("Calanda_JSDM", "output", "traits.csv"))
cat("\nSaved imputed traits to output/traits.csv\n")

# ==============================================================================
# 7. SUMMARY
# ==============================================================================
cat("\n=== Summary ===\n")
cat("TRY individual observations:", nrow(try_individual), "\n")
cat("Field individual observations:", nrow(field_individual), "\n")
cat("Combined individual observations:", nrow(all_individual), "\n")
cat("Species in final traits:", nrow(traits_imputed), "\n")

cat("\nTrait coverage in final imputed dataset:\n")
traits_imputed %>%
  summarize(across(everything(), ~ round(100 * mean(!is.na(.)), 1))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "percent_complete") %>%
  arrange(desc(percent_complete)) %>%
  print(n = 30)

cat("\n=== Script completed successfully ===\n")
