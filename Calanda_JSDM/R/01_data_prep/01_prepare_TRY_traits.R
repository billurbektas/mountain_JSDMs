# ==============================================================================
# Script: 01_prepare_TRY_traits.R
# Author: Billur Bektas
# Claude (Anthropic) was used to assist with code refactoring, validation, and documentation.
#
# Purpose: Fetch TRY trait data and ecological indicators using tidyTRY
#
# Description:
#   Uses the tidyTRY package to:
#   1. Extract species list from raw vegetation survey data
#   2. Process raw TRY database files (.txt exports)
#   3. Download ecological indicators from floraveg.eu
#   4. Download dispersal traits from floraveg.eu
#   5. Output individual-level TRY observations (for merging with field data)
#
#   Species-level summaries and imputation are done in 03_merge_traits.R
#   after combining TRY data with field-collected traits.
#
# Input files:
#   - data/vegetation/2024_CAPHE_SpeDis_CleanData_20240214.csv (species list)
#   - data/try/try_*.txt (raw TRY database exports)
#
# Output files:
#   - output/try_traits_individual.csv (individual observations for ITV)
#   - output/try_trait_units.csv (TRY standardized unit per trait)
#   - output/indicators.csv (ecological indicators from floraveg.eu)
#   - output/dispersal.csv (dispersal traits from floraveg.eu)
#
# Requires:
#   - tidyTRY package (devtools::install_github("billurbektas/tidyTRY"))
# ==============================================================================

library(tidyverse)
library(here)
library(devtools)

#devtools::install_github("billurbektas/tidyTRY", force = TRUE, upgrade = "always")

library(tidyTRY)

# ==============================================================================
# 1. EXTRACT SPECIES LIST FROM VEGETATION DATA
# ==============================================================================
cat("\n=== Extracting species list from vegetation data ===\n")

# Read raw vegetation data to get species list (avoids dependency on Y matrix)
veg_raw = read_csv(here("Calanda_JSDM", "data", "vegetation",
                        "2024_CAPHE_SpeDis_CleanData_20240214.csv"),
                   show_col_types = FALSE)

# Extract unique species names and clean them
species_list = veg_raw %>%
  pull(taxon_global) %>%
  unique() %>%
  stringi::stri_trans_general("Latin-ASCII") %>%
  sort()

cat("Found", length(species_list), "unique species in vegetation data\n")

# ==============================================================================
# 2. INSPECT AVAILABLE TRAITS IN TRY DATA
# ==============================================================================
cat("\n=== Inspecting TRY database files ===\n")

# Path to raw TRY data files
try_path = here("Calanda_JSDM", "data", "try")

# List available TRY files
try_files = list_try_files(try_path)
print(try_files)

# # Read a sample to inspect available traits
# cat("\nReading TRY data to inspect available traits...\n")
# try_sample = read_try(try_path, species = species_list[1:10], progress = FALSE)
# 
# # Get trait information
# trait_info = get_trait_info(try_sample)
# cat("\nAvailable traits in TRY data:\n")
# print(trait_info, n = 50)
# 
# # Show trait availability matrix for sample species
# cat("\nTrait availability for sample species:\n")
# trait_avail = trait_availability(try_sample)
# print(trait_avail)
# 
# # Suggest trait mapping template (useful for adding new traits)
# cat("\nSuggested trait map template:\n")
# suggest_trait_map(try_sample) %>% print(n = 30)
# 
# cat("\n--- Review traits above, then proceed with processing ---\n\n")

# ==============================================================================
# 3. PROCESS RAW TRY DATA
# ==============================================================================
cat("\n=== Processing TRY database files ===\n")

# Define traits of interest
target_traits = c(
  "N_percent",           # Leaf nitrogen content
  "C_percent",           # Leaf carbon content
  "LDMC",                # Leaf dry matter content
  "SLA",                 # Specific leaf area
  "LA",                  # Leaf area
  "vegetative_height",   # Plant vegetative height
  "reproductive_height", # Plant generative/reproductive height
  "seed_mass",          # Seed dry mass
  "flowering_phenology"
)

# Create trait mapping for tidyTRY
trait_mapping = trait_map(

  N_percent = c(
    "Leaf nitrogen (N) content per leaf dry mass",
    "Leaf nitrogen content per leaf dry mass"
  ),
  C_percent = c(
    "Leaf carbon (C) content per leaf dry mass",
    "Leaf carbon content per leaf dry mass"
  ),
  LDMC = c(
    "Leaf dry mass per leaf fresh mass (leaf dry matter content, LDMC)",
    "Leaf dry mass per leaf fresh mass (LDMC)"
  ),
  SLA = c(
    "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded",
    "Leaf area per leaf dry mass (SLA or 1/LMA): petiole excluded",
    "Leaf area per leaf dry mass (SLA or 1/LMA): petiole included",
    "Leaf area per leaf dry mass (SLA)",
    "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded",
    "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole included",
    "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA) petiole, rhachis and midrib excluded"
  ),
  LA = c(
    "Leaf area (in case of compound leaves: leaflet, undefined if petiole is in- or excluded)",
    "Leaf area",
    "Leaf area (in case of compound leaves: leaf)",
    "Leaf area (in case of compound leaves: leaflet, petiole excluded)",
    "Leaf area (in case of compound leaves: leaf, petiole excluded)",
    "Leaf area (in case of compound leaves: leaflet, petiole included)",
    "Leaf area (in case of compound leaves: leaf, petiole included)",
    "Leaf area (in case of compound leaves undefined if leaf or leaflet, undefined if petiole is in- or excluded)",
    "Leaf area (in case of compound leaves: leaf, undefined if petiole in- or excluded)"
  ),
  vegetative_height = c(
    "Plant height vegetative",
    "Plant vegetative height"
  ),
  reproductive_height = c(
    "Plant height generative",
    "Plant generative height"
  ),
  seed_mass = c(
    "Seed dry mass",
    "Seed mass"
  ),
  flowering_phenology = c(
    "Plant reproductive phenology timing (flowering time)"
  )
)

# Process TRY data with tidyTRY
cat("\nProcessing TRY files for", length(species_list), "species...\n")

try_result = process_try(
  qualitative_ids = c(335),
  files = try_path,
  species = species_list,
  trait_map = trait_mapping,
  resolve_taxonomy = TRUE,
  resolve_method = "tnrs",
  max_error_risk = 4,
  extract_location = TRUE,
  chunk_size = 100000L
)

# Extract quantitative traits
try_quanti = try_result$quantitative
cat("Extracted", nrow(try_quanti), "trait observations\n")
cat("Species with TRY data:", n_distinct(try_quanti$AccSpeciesName), "\n")

# ==============================================================================
# 4. EXTRACT CLIMATE ZONES (optional)
# ==============================================================================
cat("\n=== Extracting climate zones ===\n")

coords_with_climate = extract_climate_zones(try_result$coordinates)

# Merge climate zones back to trait data
try_quanti = try_quanti %>%
  left_join(
    coords_with_climate %>% select(ObservationID, Climate_code = climate_code),
    by = "ObservationID"
  )

cat("Climate zones assigned\n")


# ==============================================================================
# 5. DOWNLOAD ECOLOGICAL INDICATORS FROM FLORAVEG.EU
# ==============================================================================
cat("\n=== Downloading ecological indicators from floraveg.eu ===\n")

indicators = read_indicators(
  species = try_result$taxonomy,
  indicators = c("Light", "Moisture", "Nutrients", "Disturbance.Severity")
)

cat("Indicators retrieved for", nrow(indicators), "species\n")

# Rename for consistency
indicators = indicators %>%
  rename(
    species = species_original,
    species_TNRS = species_resolved,
    disturbance = Disturbance.Severity
  ) %>%
  select(species, species_TNRS, Light, Moisture, Nutrients, disturbance)

# Save indicators
write_csv(indicators, here("Calanda_JSDM", "output", "indicators.csv"))
cat("Saved indicators to output/indicators.csv\n")

# ==============================================================================
# 6. DOWNLOAD DISPERSAL TRAITS FROM FLORAVEG.EU
# ==============================================================================
cat("\n=== Downloading dispersal traits from floraveg.eu ===\n")

dispersal = read_dispersal(
  species = try_result$taxonomy
)

cat("Dispersal data retrieved for", nrow(dispersal), "species\n")

# Rename for consistency
dispersal = dispersal %>%
  rename(
    species = species_original,
    species_TNRS = species_resolved,
    dispersal = dispersal_distance_class
  ) %>%
  select(species, species_TNRS, dispersal)

# Save dispersal
write_csv(dispersal, here("Calanda_JSDM", "output", "dispersal.csv"))
cat("Saved dispersal to output/dispersal.csv\n")

# ==============================================================================
# 7. FILTER BY CLIMATE ZONE
# ==============================================================================
cat("\n=== Filtering by climate zone ===\n")

# Exclude tropical and arid climates (keep temperate, continental, polar)
excluded_climates = c("Af", "Am", "Aw", "BWh", "BWk", "BSh", "BSk")

try_filtered = try_quanti %>%
  filter(!(Climate_code %in% excluded_climates | is.na(Climate_code)))

cat("After climate filtering:", nrow(try_filtered), "observations\n")

# ==============================================================================
# 8. CREATE INDIVIDUAL-LEVEL OUTPUT (for ITV analysis)
# ==============================================================================
cat("\n=== Creating individual-level trait dataset ===\n")

# Keep individual observations with relevant columns
try_individual = try_filtered %>%
  select(
    species = SpeciesName,
    species_TNRS = AccSpeciesName,
    TraitName,
    Value = StdValue,
    ObservationID,
    Climate_code,
    Dataset = DatasetID,
    Unit = UnitName
  ) %>%
  filter(!is.na(Value))

write_csv(try_individual, here("Calanda_JSDM", "output", "try_traits_individual.csv"))
cat("Saved", nrow(try_individual), "individual observations to output/try_traits_individual.csv\n")

# Save TRY unit reference (one row per trait, majority unit from StdValue)
write_csv(try_result$units, here("Calanda_JSDM", "output", "try_trait_units.csv"))
cat("Saved TRY trait unit reference to output/try_trait_units.csv\n")

# ==============================================================================
# 9. SUMMARY
# ==============================================================================
cat("\n=== Summary ===\n")
cat("Species in vegetation data:", length(species_list), "\n")
cat("TRY individual observations:", nrow(try_individual), "\n")
cat("Species with TRY data:", n_distinct(try_individual$species_TNRS), "\n")
cat("Species with indicators:", n_distinct(indicators$species_TNRS), "\n")
cat("Species with dispersal:", n_distinct(dispersal$species_TNRS), "\n")

cat("\nOutput files:\n")
cat("  - output/try_traits_individual.csv (for merging with field data)\n")
cat("  - output/try_trait_units.csv (TRY standardized units per trait)\n")
cat("  - output/indicators.csv\n")
cat("  - output/dispersal.csv\n")

cat("\nNext step: Run 02_prepare_trait_data.R for field traits,")
cat("\n           then 03_merge_traits.R to combine and calculate species summaries.\n")

cat("\n=== Script completed successfully ===\n")
