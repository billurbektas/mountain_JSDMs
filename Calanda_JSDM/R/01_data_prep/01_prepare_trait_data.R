# ==============================================================================
# Script: 01_prepare_trait_data.R
# Purpose: Comprehensive trait data preparation and analysis pipeline
#
# Description:
#   This script performs complete data processing from raw Excel files to
#   analysis-ready datasets including PCA and correlation analyses.
#
#   Processing steps:
#   1. Extract N and C isotope data from Excel files (N+C_Samples*.xls, Trial*.xls)
#      - Handles multi-section sheets with repeating headers
#      - Filters numeric sample IDs only (removes calibration samples)
#   2. Map CN analysis IDs to plant IDs using sample ID mapping file
#   3. Merge trait datasets: traits_clean, leaf_area, and nutrients
#   4. Calculate derived traits: LMA, LDMC, height ratios, LNC/LCC ratio
#   5. Validate data quality:
#      - Check dry < fresh weight measurements
#      - Flag outliers in nutrient data
#   6. Handle duplicate measurements by averaging
#   7. Create two output versions:
#      - Full dataset with error flags
#      - Clean dataset (only problematic values set to NA, not entire rows)
#   8. Calculate species-level summaries with mean and kCV (robust CV)
#   9. Generate correlation plot (all trait means and kCV values combined)
#   10. Perform PCA on trait means and kCV values separately
#       - Biplots with abbreviated species names (e.g., "GeumMont")
#       - Species colored by sample size
#
# Input files:
#   - data/traits/N+C_Samples*.xls (isotope data)
#   - data/traits/Trial*.xls (isotope data)
#   - data/traits/2025_CAPHE_traits_sample_id.csv (ID mapping)
#   - data/traits/2025_TRAITS_clean_20251023.csv (cleaned trait data)
#   - data/traits/leaf_area.csv (leaf area measurements)
#
# Output files:
#   - output/all_nc_samples_cleaned.csv
#   - output/leaf_nutrients_cleaned.csv
#   - output/all_traits_merged.csv
#   - output/final_traits_for_analysis.csv (with error flags)
#   - output/final_traits_clean.csv (clean version)
#   - output/species_trait_summary.csv
#   - output/pca_*_eigenvalues.csv and pca_*_var_contrib.csv
#   - output/correlation_all_traits.csv
#   - plot/corrplot_all_traits.pdf
#   - plot/pca_trait_means_biplot.pdf
#   - plot/pca_kcv_biplot.pdf
#
# ==============================================================================

# Load required libraries ----
library(readxl)
library(tidyverse)
library(corrplot)
library(FactoMineR)
library(factoextra)
library(ggrepel)
library(here)

# Function to clean N and C data from Excel files ----
clean_nc_data = function(file_path) {

  # Extract file name for identification
  file_name = basename(file_path)

  # Read entire sheets without any skipping to identify all sections
  d15n_raw = read_excel(file_path, sheet = "d15N", col_names = FALSE)
  d13c_raw = read_excel(file_path, sheet = "d13C", col_names = FALSE)

  # Function to find rows that contain "Sample #" (header rows)
  find_header_rows = function(df) {
    which(apply(df, 1, function(row) any(grepl("Sample #", row, fixed = TRUE))))
  }

  # Find all header row positions
  header_rows_n = find_header_rows(d15n_raw)
  header_rows_c = find_header_rows(d13c_raw)

  # Function to extract data from a section
  extract_section = function(df, header_row) {
    # The actual column names are in the header_row
    # The sub-headers (like "%N corr.") are in header_row + 1
    col_names_row1 = as.character(df[header_row, ])
    col_names_row2 = as.character(df[header_row + 1, ])

    # Find column indices
    sample_col = which(col_names_row1 == "Sample #")[1]

    # Try to find %N/C corr in row2 (sub-header) first, then in row1 (main header)
    corr_col = which(grepl("%N corr\\.|%C corr\\.", col_names_row2))[1]
    if (is.na(corr_col)) {
      corr_col = which(grepl("%N corr\\.|%C corr\\.", col_names_row1))[1]
    }

    # If %N/C corr column not found, skip this section
    if (is.na(corr_col)) {
      return(NULL)
    }

    # Determine the end of this section (next header or end of data)
    next_header_idx = which(header_rows_n > header_row)[1]
    if (is.na(next_header_idx)) {
      end_row = nrow(df)
    } else {
      end_row = header_rows_n[next_header_idx] - 1
    }

    # Extract data rows (skip 2 rows after header for sub-headers)
    data_start = header_row + 3
    if (data_start > end_row) {
      return(NULL)
    }

    # Extract the relevant columns
    section_data = df[data_start:end_row, c(sample_col, corr_col)]
    colnames(section_data) = c("sample_number", "corr_value")

    # Clean the data
    section_data = section_data %>%
      filter(!is.na(sample_number)) %>%
      filter(sample_number != "NA") %>%
      filter(!is.na(corr_value)) %>%
      mutate(corr_value = as.numeric(corr_value)) %>%
      filter(!is.na(corr_value))

    return(section_data)
  }

  # Extract d15N data from all sections
  d15n_sections = map(header_rows_n, ~ extract_section(d15n_raw, .))
  d15n_clean = bind_rows(d15n_sections) %>%
    rename(n_corr = corr_value) %>%
    mutate(source_file = file_name)

  # Extract d13C data from all sections
  d13c_sections = map(header_rows_c, ~ extract_section(d13c_raw, .))
  d13c_clean = bind_rows(d13c_sections) %>%
    rename(c_corr = corr_value) %>%
    mutate(source_file = file_name)

  # Merge d15N and d13C data by sample number
  # Add row numbers to handle duplicates properly
  d15n_clean = d15n_clean %>% mutate(row_id = row_number())
  d13c_clean = d13c_clean %>% mutate(row_id = row_number())

  nc_combined = full_join(d15n_clean, d13c_clean,
                           by = c("sample_number", "source_file", "row_id")) %>%
    select(-row_id)

  # Keep only numeric sample numbers (removes calibration samples like Ali, Caf, Tyr, etc.)
  nc_combined = nc_combined %>%
    filter(grepl("^[0-9]+$", sample_number))

  return(nc_combined)
}

# Process N+C_Samples13.xls ----
nc_data_13 = clean_nc_data(here("Calanda_JSDM", "data", "traits", "N+C_Samples13.xls"))

# Preview the cleaned data
print(nc_data_13)

# Save cleaned data
write_csv(nc_data_13, here("Calanda_JSDM", "output", "nc_samples_13_cleaned.csv"))

cat("Data cleaning completed for N+C_Samples13.xls\n")
cat("Cleaned data saved to: output/nc_samples_13_cleaned.csv\n")
cat("Number of samples:", nrow(nc_data_13), "\n")

# Process all N+C_Samples files ----
cat("\n--- Processing all N+C_Samples files ---\n")

# Find all N+C_Samples*.xls files
nc_files = c(list.files(here("Calanda_JSDM", "data", "traits"),
                       pattern = "N\\+C_Samples.*\\.xls$",
                       full.names = TRUE),
              list.files(here("Calanda_JSDM", "data", "traits"),
                         pattern = "Trial",
                         full.names = TRUE)
             )

cat("Found", length(nc_files), "files:\n")
print(basename(nc_files))

# Process all files
all_nc_data = map_df(nc_files, function(file) {
  cat("\nProcessing:", basename(file), "...")
  result = clean_nc_data(file)
  cat(" Done! (", nrow(result), "samples )\n")
  return(result)
})

# Read the sample ID mapping file to match cn_ID to plant_ID
sample_id_mapping = read_csv(here("Calanda_JSDM", "data", "traits", "2025_CAPHE_traits_sample_id.csv"))

cat("\n--- Matching CN IDs to Plant IDs ---\n")
cat("Sample ID mapping file:", nrow(sample_id_mapping), "rows\n")
cat("CN data samples:", nrow(all_nc_data), "samples\n")

# Rename and prepare CN data
all_nc_data =
  all_nc_data %>%
  mutate(sample_number = as.integer(sample_number))%>%
  rename(cn_ID = sample_number, LNC = n_corr, LCC = c_corr)

# Match cn_ID to plant_ID using the mapping file
# Convert cn_ID to integer in mapping file to match
all_nc_data = all_nc_data %>%
  left_join(sample_id_mapping %>%
              mutate(cn_ID = as.integer(cn_ID)) %>%
              filter(!is.na(cn_ID)) %>%  # Remove rows where cn_ID cannot be converted
              select(cn_ID, plant_ID, plant_species, site),
            by = "cn_ID")

cat("Samples matched to plant_ID:", sum(!is.na(all_nc_data$plant_ID)), "\n")
cat("Samples without plant_ID match:", sum(is.na(all_nc_data$plant_ID)), "\n")

# See if there is anything wrong with the data
all_nc_data %>%
  select(LCC, LNC, cn_ID)%>%
  pivot_longer(cols = !cn_ID)%>%
  ggplot(aes(value))+
    facet_grid(.~name, scales = "free_x")+
    geom_histogram()

outlier_c =
  all_nc_data %>%
  filter(LCC > 70 | LCC < 20) %>% # Above/below these numbers probably measurement error.
  pull(cn_ID)

#LNC seems to be OK but would be good to check against TRYDB

nutrients =
  all_nc_data %>%
  mutate(keep = ifelse(cn_ID %in% outlier_c, FALSE, TRUE))%>%
  rename(plant_id = plant_ID)%>%
  filter(!is.na(plant_id))%>%  # Keep only rows with valid plant_id
  mutate(plant_id = as.numeric(plant_id))  # Convert plant_id to numeric to match traits_clean

# Check for duplicate measurements (same plant_id measured multiple times)
cat("\n--- Handling duplicate nutrient measurements ---\n")
duplicates = nutrients %>%
  group_by(plant_id) %>%
  filter(n() > 1) %>%
  ungroup()

if(nrow(duplicates) > 0) {
  cat("Found", n_distinct(duplicates$plant_id), "plant_ids with duplicate measurements\n")
  cat("Total duplicate measurements:", nrow(duplicates), "\n")

  # Average nutrient values for duplicates
  nutrients = nutrients %>%
    group_by(plant_id) %>%
    summarize(
      LNC = mean(LNC, na.rm = TRUE),
      LCC = mean(LCC, na.rm = TRUE),
      keep = all(keep),  # Keep FALSE if any measurement is flagged as outlier
      source_file = paste(unique(source_file), collapse = "; "),  # Combine source files
      plant_species = first(plant_species),
      site = first(site),
      .groups = "drop"
    )

  cat("After averaging duplicates:", nrow(nutrients), "unique plant_ids\n")
} else {
  cat("No duplicate measurements found\n")
}

# Save combined data
write_csv(nutrients, here("Calanda_JSDM", "output", "leaf_nutrients_cleaned.csv"))

cat("\n=== Summary ===\n")
cat("Total samples across all files:", nrow(nutrients), "\n")
cat("Unique sample numbers:", n_distinct(nutrients$plant_id), "\n")
cat("Files processed:", length(nc_files), "\n")
cat("\nCombined data saved to: output/all_nc_samples_cleaned.csv\n")

# Show summary statistics
cat("\n--- Summary statistics ---\n")
cat("LNC range:", min(nutrients$LNC[nutrients$keep], na.rm = TRUE), "to",
    max(nutrients$LNC[nutrients$keep], na.rm = TRUE), "\n")
cat("LCC range:", min(nutrients$LCC[nutrients$keep], na.rm = TRUE), "to",
    max(nutrients$LCC[nutrients$keep], na.rm = TRUE), "\n")

### NOW let's merge all the traits data ----

# Read traits_clean data
traits_clean = read_csv(here("Calanda_JSDM", "data", "traits", "2025_TRAITS_clean_20251023.csv"))

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

# Remove temporary diff columns
traits_clean = traits_clean %>%
  select(-biomass_weight_diff, -flowering_unit_weight_diff, -leaves_weight_diff)

# Read leaf_area data
leaf_area = read_csv(here("Calanda_JSDM", "data", "traits", "leaf_area.csv"))

# Check dimensions
cat("\n--- Merging datasets ---\n")
cat("traits_clean dimensions:", nrow(traits_clean), "rows x", ncol(traits_clean), "columns\n")
cat("leaf_area dimensions:", nrow(leaf_area), "rows x", ncol(leaf_area), "columns\n")
cat("nutrients dimensions:", nrow(nutrients), "rows x", ncol(nutrients), "columns\n")

# Merge all three datasets by plant_id
all_traits = traits_clean %>%
  left_join(leaf_area %>%
              select(plant_id, n_scans, n_leaves, sum_leaf_cm2, mean_leaf_cm2, LMA),
            by = "plant_id") %>%
  left_join(nutrients %>%
              select(plant_id, LNC, LCC, keep, source_file),
            by = "plant_id")

# Show merge results
cat("\n--- Merge results ---\n")
cat("Final merged dataset:", nrow(all_traits), "rows x", ncol(all_traits), "columns\n")
cat("Samples with leaf_area data:", sum(!is.na(all_traits$LMA)), "\n")
cat("Samples with nutrient data:", sum(!is.na(all_traits$LNC)), "\n")
cat("Samples with both leaf_area and nutrients:",
    sum(!is.na(all_traits$LMA) & !is.na(all_traits$LNC)), "\n")

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

cat("LDMC calculated for", sum(!is.na(all_traits$LDMC)), "samples\n")
cat("LDMC range:", round(min(all_traits$LDMC, na.rm=TRUE), 2), "to",
    round(max(all_traits$LDMC, na.rm=TRUE), 2), "mg/g\n")

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
    LNC,              # Leaf Nitrogen Content (corrected %)
    LCC,              # Leaf Carbon Content (corrected %)
    n_leaves,         # Number of leaves scanned
    sum_leaf_cm2,     # Total leaf area

    # Quality flags
    keep,                        # Nutrient outlier flag (FALSE = outlier)
    weight_measurement_error,    # Weight measurement error flag
    source_file                  # CN analysis source file
  )

# Summary of final dataset
cat("\n=== Final Trait Dataset Summary ===\n")
cat("Total samples:", nrow(final_traits), "\n")
cat("Total variables:", ncol(final_traits), "\n\n")

cat("Data completeness:\n")
cat("  - Plant height (vegetative):", sum(!is.na(final_traits$vegetative_height)),
    "/", nrow(final_traits),
    " (", round(100*sum(!is.na(final_traits$vegetative_height))/nrow(final_traits), 1), "%)\n", sep="")
cat("  - Leaf fresh weight:", sum(!is.na(final_traits$leaves_fresh_weight)),
    "/", nrow(final_traits),
    " (", round(100*sum(!is.na(final_traits$leaves_fresh_weight))/nrow(final_traits), 1), "%)\n", sep="")
cat("  - Leaf dry weight:", sum(!is.na(final_traits$leaves_dry_weight)),
    "/", nrow(final_traits),
    " (", round(100*sum(!is.na(final_traits$leaves_dry_weight))/nrow(final_traits), 1), "%)\n", sep="")
cat("  - LMA:", sum(!is.na(final_traits$LMA)),
    "/", nrow(final_traits),
    " (", round(100*sum(!is.na(final_traits$LMA))/nrow(final_traits), 1), "%)\n", sep="")
cat("  - LNC:", sum(!is.na(final_traits$LNC)),
    "/", nrow(final_traits),
    " (", round(100*sum(!is.na(final_traits$LNC))/nrow(final_traits), 1), "%)\n", sep="")
cat("  - LCC:", sum(!is.na(final_traits$LCC)),
    "/", nrow(final_traits),
    " (", round(100*sum(!is.na(final_traits$LCC))/nrow(final_traits), 1), "%)\n", sep="")
cat("  - LDMC:", sum(!is.na(final_traits$LDMC)),
    "/", nrow(final_traits),
    " (", round(100*sum(!is.na(final_traits$LDMC))/nrow(final_traits), 1), "%)\n", sep="")

cat("\nSamples with complete trait data (all 7 key traits):",
    sum(complete.cases(final_traits[, c("vegetative_height", "leaves_fresh_weight", "leaves_dry_weight", "LMA", "LDMC", "LNC", "LCC")])),
    "\n")

cat("\nQuality flags:\n")
cat("  - Weight measurement errors:", sum(final_traits$weight_measurement_error, na.rm=TRUE), "\n")
cat("  - Nutrient outliers (keep=FALSE):", sum(!final_traits$keep, na.rm=TRUE), "\n")

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
    LNC = ifelse(!is.na(keep) & !keep, NA_real_, LNC),
    LCC = ifelse(!is.na(keep) & !keep, NA_real_, LCC)
  )

cat("Removed problematic values:\n")
cat("  - LMA/LDMC values removed due to weight errors:",
    sum(final_traits$weight_measurement_error & !is.na(final_traits$LMA), na.rm=TRUE), "\n")
cat("  - LNC/LCC values removed due to nutrient outliers:",
    sum(!is.na(final_traits$keep) & !final_traits$keep & !is.na(final_traits$LNC), na.rm=TRUE), "\n")

cat("\nCleaned dataset - Data completeness:\n")
cat("  - LMA:", sum(!is.na(final_traits_clean$LMA)),
    "/", nrow(final_traits_clean),
    " (", round(100*sum(!is.na(final_traits_clean$LMA))/nrow(final_traits_clean), 1), "%)\n", sep="")
cat("  - LDMC:", sum(!is.na(final_traits_clean$LDMC)),
    "/", nrow(final_traits_clean),
    " (", round(100*sum(!is.na(final_traits_clean$LDMC))/nrow(final_traits_clean), 1), "%)\n", sep="")
cat("  - LNC:", sum(!is.na(final_traits_clean$LNC)),
    "/", nrow(final_traits_clean),
    " (", round(100*sum(!is.na(final_traits_clean$LNC))/nrow(final_traits_clean), 1), "%)\n", sep="")
cat("  - LCC:", sum(!is.na(final_traits_clean$LCC)),
    "/", nrow(final_traits_clean),
    " (", round(100*sum(!is.na(final_traits_clean$LCC))/nrow(final_traits_clean), 1), "%)\n", sep="")

cat("\nSamples with complete trait data in cleaned dataset:",
    sum(complete.cases(final_traits_clean[, c("vegetative_height", "LMA", "LDMC", "LNC", "LCC")])),
    "\n")

# Save all datasets
write_csv(all_traits, here("Calanda_JSDM", "output", "all_traits_merged.csv"))
write_csv(final_traits, here("Calanda_JSDM", "output", "final_traits_for_analysis.csv"))
write_csv(final_traits_clean, here("Calanda_JSDM", "output", "final_traits_clean.csv"))

cat("\n=== Output files ===\n")
cat("Full merged data: output/all_traits_merged.csv\n")
cat("Final trait dataset (with error flags): output/final_traits_for_analysis.csv\n")
cat("Clean trait dataset (errors removed): output/final_traits_clean.csv\n")

# Calculate species-level trait summaries with kCV ----
cat("\n--- Calculating species-level trait summaries ---\n")

# Define kCV function
kcv = function(x, na.rm = TRUE) {
  # remove NAs if requested
  if (na.rm) x = x[!is.na(x)]

  # need at least 2 values
  if (length(x) < 2L) return(NA_real_)

  m  = mean(x)
  s  = stats::sd(x)

  # avoid division by zero
  if (m == 0) return(NA_real_)

  cv = s / m
  kcv = sqrt(cv^2 / (1 + cv^2))
  return(kcv)
}

# Calculate ratios and summarize by species
species_traits = final_traits_clean %>%
  mutate(
    # Calculate ratios
    repro_veg_ratio = reproductive_height / vegetative_height,
    LNC_LCC_ratio = LNC / LCC
  ) %>%
  group_by(plant_species) %>%
  summarize(
    n_samples = n(),

    # Vegetative height
    vegetative_height_mean = mean(vegetative_height, na.rm = TRUE),
    vegetative_height_kcv = kcv(vegetative_height, na.rm = TRUE),
    vegetative_height_n = sum(!is.na(vegetative_height)),

    # Reproductive height
    reproductive_height_mean = mean(reproductive_height, na.rm = TRUE),
    reproductive_height_kcv = kcv(reproductive_height, na.rm = TRUE),
    reproductive_height_n = sum(!is.na(reproductive_height)),

    # Reproductive/Vegetative ratio
    repro_veg_ratio_mean = mean(repro_veg_ratio, na.rm = TRUE),
    repro_veg_ratio_kcv = kcv(repro_veg_ratio, na.rm = TRUE),
    repro_veg_ratio_n = sum(!is.na(repro_veg_ratio)),

    # LNC
    LNC_mean = mean(LNC, na.rm = TRUE),
    LNC_kcv = kcv(LNC, na.rm = TRUE),
    LNC_n = sum(!is.na(LNC)),

    # LCC
    LCC_mean = mean(LCC, na.rm = TRUE),
    LCC_kcv = kcv(LCC, na.rm = TRUE),
    LCC_n = sum(!is.na(LCC)),

    # LMA
    LMA_mean = mean(LMA, na.rm = TRUE),
    LMA_kcv = kcv(LMA, na.rm = TRUE),
    LMA_n = sum(!is.na(LMA)),

    # LDMC
    LDMC_mean = mean(LDMC, na.rm = TRUE),
    LDMC_kcv = kcv(LDMC, na.rm = TRUE),
    LDMC_n = sum(!is.na(LDMC)),

    # LNC/LCC ratio
    LNC_LCC_ratio_mean = mean(LNC_LCC_ratio, na.rm = TRUE),
    LNC_LCC_ratio_kcv = kcv(LNC_LCC_ratio, na.rm = TRUE),
    LNC_LCC_ratio_n = sum(!is.na(LNC_LCC_ratio)),

    .groups = "drop"
  )

cat("Species-level summary calculated for", nrow(species_traits), "species\n")
cat("Traits included: vegetative_height, reproductive_height, repro/veg ratio,\n")
cat("                 LNC, LCC, LMA, LDMC, LNC/LCC ratio\n")

# Save species-level summary
write_csv(species_traits, here("Calanda_JSDM", "output", "species_trait_summary.csv"))

cat("\nSpecies trait summary saved to: output/species_trait_summary.csv\n")

# Show preview of species summary
cat("\nPreview of species summary (first 5 species):\n")
print(species_traits %>%
        select(plant_species, n_samples,
               vegetative_height_mean, vegetative_height_kcv,
               LMA_mean, LMA_kcv, LNC_mean, LNC_kcv) %>%
        head(5))

# Correlation plots ----
cat("\n--- Creating correlation plots ---\n")

# Combine means and kCV data for comprehensive correlation analysis
all_traits_data = species_traits %>%
  select(plant_species,
         # Means
         vegetative_height_mean, reproductive_height_mean, repro_veg_ratio_mean,
         LNC_mean, LCC_mean, LMA_mean, LDMC_mean, LNC_LCC_ratio_mean,
         # kCV values
         vegetative_height_kcv, reproductive_height_kcv, repro_veg_ratio_kcv,
         LNC_kcv, LCC_kcv, LMA_kcv, LDMC_kcv, LNC_LCC_ratio_kcv) %>%
  filter(complete.cases(.)) %>%  # Keep only complete cases
  column_to_rownames("plant_species")

cat("Correlation matrix calculated for", nrow(all_traits_data), "species with complete data\n")
cat("Including", ncol(all_traits_data), "variables (8 means + 8 kCV values)\n")

# Calculate correlation matrix
cor_all = cor(all_traits_data, use = "complete.obs")

# Create comprehensive correlation plot
pdf(here("Calanda_JSDM", "plot", "corrplot_all_traits.pdf"), width = 14, height = 14)
corrplot(cor_all,
         method = "color",
         type = "upper",
         order = "hclust",
         addCoef.col = "black",
         number.cex = 0.5,
         tl.col = "black",
         tl.srt = 45,
         tl.cex = 0.8,
         diag = FALSE,
         title = "Correlation - All Trait Means and kCV Values",
         mar = c(0, 0, 2, 0))
dev.off()

# Save correlation matrix
write_csv(as.data.frame(cor_all) %>% rownames_to_column("trait"),
          here("Calanda_JSDM", "output", "correlation_all_traits.csv"))

cat("\nCorrelation plot saved to plot/corrplot_all_traits.pdf\n")
cat("Correlation matrix saved to output/correlation_all_traits.csv\n")

# PCA Analysis ----
cat("\n--- Running PCA on trait means and kCV values ---\n")

# Function to abbreviate species names (4 letters genus + 4 letters species)
# e.g., "Geum montanum" -> "GeumMont"
abbreviate_species = function(species_name) {
  parts = strsplit(species_name, " ")[[1]]
  if (length(parts) >= 2) {
    genus = substr(parts[1], 1, 4)
    species = substr(parts[2], 1, 4)
    # Capitalize first letter of species epithet
    species = paste0(toupper(substr(species, 1, 1)), substr(species, 2, 4))
    return(paste0(genus, species))
  } else {
    return(substr(species_name, 1, 8))
  }
}

# Prepare data for PCA (remove species with too many missing values)
# PCA 1: Mean trait values
pca_means_data = species_traits %>%
  select(plant_species,
         vegetative_height_mean,
         repro_veg_ratio_mean,
         LNC_mean, LCC_mean, LMA_mean, LDMC_mean) %>%
  filter(complete.cases(.)) %>%  # Keep only complete cases
  column_to_rownames("plant_species")

cat("PCA on MEANS: Using", nrow(pca_means_data), "species with complete trait data\n")

# Run PCA on means
pca_means = PCA(pca_means_data, scale.unit = TRUE, graph = FALSE)

# PCA 2: kCV values
pca_kcv_data = species_traits %>%
  select(plant_species,
         vegetative_height_kcv,
         repro_veg_ratio_kcv,
         LNC_kcv, LCC_kcv, LMA_kcv, LDMC_kcv) %>%
  filter(complete.cases(.)) %>%  # Keep only complete cases
  column_to_rownames("plant_species")

cat("PCA on kCV: Using", nrow(pca_kcv_data), "species with complete kCV data\n")

# Run PCA on kCV
pca_kcv = PCA(pca_kcv_data, scale.unit = TRUE, graph = FALSE)

# Visualize PCAs - Create biplots only
cat("\nCreating PCA biplots...\n")

# Create abbreviated species names for plotting
abbreviated_means = sapply(rownames(pca_means_data), abbreviate_species)
abbreviated_kcv = sapply(rownames(pca_kcv_data), abbreviate_species)

# Get sample sizes for coloring
sample_sizes_means = species_traits %>%
  filter(plant_species %in% rownames(pca_means_data)) %>%
  select(plant_species, n_samples) %>%
  arrange(match(plant_species, rownames(pca_means_data)))

sample_sizes_kcv = species_traits %>%
  filter(plant_species %in% rownames(pca_kcv_data)) %>%
  select(plant_species, n_samples) %>%
  arrange(match(plant_species, rownames(pca_kcv_data)))

# Replace rownames with abbreviated names BEFORE plotting
rownames(pca_means$ind$coord) = abbreviated_means
pca_means$call$X = pca_means_data
rownames(pca_means$call$X) = abbreviated_means

# PCA of mean trait values - Biplot
# Create data frame for plotting with sample sizes
pca_means_coords = as.data.frame(pca_means$ind$coord[, 1:2])
pca_means_coords$species = rownames(pca_means_coords)
pca_means_coords$n_samples = sample_sizes_means$n_samples

pca_means_var = as.data.frame(pca_means$var$coord[, 1:2])
pca_means_var$variable = rownames(pca_means_var)

pdf(here("Calanda_JSDM", "plot", "pca_trait_means_biplot.pdf"), width = 12, height = 10)
p1 = ggplot() +
  # Add variable arrows
  geom_segment(data = pca_means_var,
               aes(x = 0, y = 0, xend = Dim.1 * 4, yend = Dim.2 * 4),
               arrow = arrow(length = unit(0.3, "cm")),
               color = "black", linewidth = 0.8) +
  geom_text(data = pca_means_var,
            aes(x = Dim.1 * 4.5, y = Dim.2 * 4.5, label = variable),
            color = "black", size = 4, fontface = "bold") +
  # Add species points and labels colored by n_samples
  geom_point(data = pca_means_coords,
             aes(x = Dim.1, y = Dim.2, color = n_samples),
             size = 2) +
  geom_text_repel(data = pca_means_coords,
                  aes(x = Dim.1, y = Dim.2, label = species, color = n_samples),
                  size = 3, max.overlaps = 15) +
  scale_color_viridis_c(name = "Sample size") +
  labs(title = "PCA - Mean Trait Values",
       subtitle = paste("Species:", nrow(pca_means_data)),
       x = paste0("PC1 (", round(pca_means$eig[1,2], 1), "%)"),
       y = paste0("PC2 (", round(pca_means$eig[2,2], 1), "%)")) +
  theme_bw() +
  theme(legend.position = "right")
print(p1)
dev.off()

# Replace rownames with abbreviated names BEFORE plotting
rownames(pca_kcv$ind$coord) = abbreviated_kcv
pca_kcv$call$X = pca_kcv_data
rownames(pca_kcv$call$X) = abbreviated_kcv

# PCA of kCV values - Biplot
# Create data frame for plotting with sample sizes
pca_kcv_coords = as.data.frame(pca_kcv$ind$coord[, 1:2])
pca_kcv_coords$species = rownames(pca_kcv_coords)
pca_kcv_coords$n_samples = sample_sizes_kcv$n_samples

pca_kcv_var = as.data.frame(pca_kcv$var$coord[, 1:2])
pca_kcv_var$variable = rownames(pca_kcv_var)

pdf(here("Calanda_JSDM", "plot", "pca_kcv_biplot.pdf"), width = 12, height = 10)
p2 = ggplot() +
  # Add variable arrows
  geom_segment(data = pca_kcv_var,
               aes(x = 0, y = 0, xend = Dim.1 * 4, yend = Dim.2 * 4),
               arrow = arrow(length = unit(0.3, "cm")),
               color = "black", linewidth = 0.8) +
  geom_text(data = pca_kcv_var,
            aes(x = Dim.1 * 4.5, y = Dim.2 * 4.5, label = variable),
            color = "black", size = 4, fontface = "bold") +
  # Add species points and labels colored by n_samples
  geom_point(data = pca_kcv_coords,
             aes(x = Dim.1, y = Dim.2, color = n_samples),
             size = 2) +
  geom_text_repel(data = pca_kcv_coords,
                  aes(x = Dim.1, y = Dim.2, label = species, color = n_samples),
                  size = 3, max.overlaps = 15) +
  scale_color_viridis_c(name = "Sample size") +
  labs(title = "PCA - kCV Values",
       subtitle = paste("Species:", nrow(pca_kcv_data)),
       x = paste0("PC1 (", round(pca_kcv$eig[1,2], 1), "%)"),
       y = paste0("PC2 (", round(pca_kcv$eig[2,2], 1), "%)")) +
  theme_bw() +
  theme(legend.position = "right")
print(p2)
dev.off()

# Save PCA summaries
cat("\nSaving PCA summaries...\n")

# Eigenvalues
write_csv(as.data.frame(pca_means$eig), here("Calanda_JSDM", "output", "pca_means_eigenvalues.csv"))
write_csv(as.data.frame(pca_kcv$eig), here("Calanda_JSDM", "output", "pca_kcv_eigenvalues.csv"))

# Variable contributions
write_csv(as.data.frame(pca_means$var$contrib), here("Calanda_JSDM", "output", "pca_means_var_contrib.csv"))
write_csv(as.data.frame(pca_kcv$var$contrib), here("Calanda_JSDM", "output", "pca_kcv_var_contrib.csv"))

cat("\nPCA biplots saved to plot/ directory\n")
cat("PCA summaries saved to output/ directory\n")

# Print variance explained
cat("\n=== PCA on MEAN trait values ===\n")
cat("Variance explained by PC1:", round(pca_means$eig[1,2], 2), "%\n")
cat("Variance explained by PC2:", round(pca_means$eig[2,2], 2), "%\n")
cat("Cumulative variance (PC1+PC2):", round(sum(pca_means$eig[1:2,2]), 2), "%\n")

cat("\n=== PCA on kCV values ===\n")
cat("Variance explained by PC1:", round(pca_kcv$eig[1,2], 2), "%\n")
cat("Variance explained by PC2:", round(pca_kcv$eig[2,2], 2), "%\n")
cat("Cumulative variance (PC1+PC2):", round(sum(pca_kcv$eig[1:2,2]), 2), "%\n")

cat("\n=== Script completed successfully! ===\n")
