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

