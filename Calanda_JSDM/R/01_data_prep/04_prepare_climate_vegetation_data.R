# ==============================================================================
# Script: 01_prepare_data.R
# Purpose: Prepare environmental + vegetation data for JSDM analysis
#
# Inputs:
#   - data/vegetation/2024_CAPHE_SpeDis_CleanData_20240214.csv
#   - data/ecostress/*.csv (LST data)
#   - data/modis/*.csv (ET data)
#   - data/wekeo/ (Copernicus snow data)
#   - data/vegetation/veg.coord.csv
#   - data/mask/*.shp (study region shapefiles)
#   - output/traits.csv (created by 01b_fetch_try_traits.R)
#
# Outputs:
#   - output/veg_clim.csv
#   - output/veg_tree.csv (tree/shrub cover + woody species abundances)
#   - output/starter_data_25.04.25.RData
#   - plot/land_use_variable.pdf
#   - plot/climate_calanda_surveys.pdf
#
# Requires:
#   - R/00_setup/functions_calanda.R
#   - Must run 01b_fetch_try_traits.R first to create output/traits.csv
# ==============================================================================

library(tidyverse)
library(conflicted)
library(stringi)
library(sf)
library(terra)
library(corrplot)
library(FactoMineR)
library(factoextra)
library(patchwork)
library(gt)
library(here)

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("extract", "tidyr")
conflict_prefer("intersect", "base")

source(here("Calanda_JSDM", "R", "00_setup", "functions_calanda.R"))

# ==============================================================================
# VEGETATION DATA
# ==============================================================================
cat("\n=== Loading vegetation data ===\n")

veg = read_csv(here("Calanda_JSDM", "data", "vegetation",
                     "2024_CAPHE_SpeDis_CleanData_20240214.csv")) %>%
  mutate(plot_id_releve = paste0(plot_id, releve_id)) %>%
  mutate(across(starts_with("soil_depth"), ~ as.numeric(.))) %>%
  rowwise() %>%
  mutate(soil_depth_mean = mean(c_across(starts_with("soil_depth")), na.rm = TRUE)) %>%
  mutate(soil_depth_var = var(c_across(starts_with("soil_depth")), na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(southness = abs(aspect - 180)) %>%
  mutate(altitude = rowMeans(select(., altitude_max, altitude_min), na.rm = TRUE)) %>%
  select(-altitude_min, -altitude_max) %>%
  select(plot_id_releve, x, y, slope, southness, trees_cover, shrubs_cover,
                rocks_cover, soil_depth_mean, soil_depth_var, open_close,
                taxon_global, species_cover, altitude) %>%
  distinct() %>%
  rename(Latitude = y, Longitude = x) %>%
  mutate(taxon_global = stri_trans_general(taxon_global, "Latin-ASCII"))%>%
  mutate(soil_depth_var = ifelse(is.nan(soil_depth_var), NA, soil_depth_var)) # happens because there are NAs for soil depth.

cat(length(unique(veg$taxon_global)), "species and",
    length(unique(veg$plot_id_releve)), "plots.\n")

# ==============================================================================
# ENVIRONMENTAL DATA
# ==============================================================================

if (file.exists(here("Calanda_JSDM", "output", "veg_clim.csv"))) {
  veg_clim = read_csv(here("Calanda_JSDM", "output", "veg_clim.csv")) %>%
    filter(fdd > -1000) %>%
    distinct()
} else {

  # Calanda study region mask ----
  if (file.exists(here("Calanda_JSDM", "output", "calanda_mask.shp"))) {
    calanda_mask = st_read(here("Calanda_JSDM", "output", "calanda_mask.shp"))
  } else {
    calanda_1 = st_read(here("Calanda_JSDM", "data", "mask", "study_region_2024.shp"))
    calanda_1 = st_transform(calanda_1, crs = "EPSG:4326")

    calanda_2 = st_read(here("Calanda_JSDM", "data", "mask", "study_region.shp"))
    calanda_2 = st_transform(calanda_2, crs = "EPSG:4326")
    common_cols = intersect(names(calanda_1), names(calanda_2))
    calanda_1 = calanda_1[, common_cols]
    calanda_2 = calanda_2[, common_cols]
    calanda_combined = rbind(calanda_1, calanda_2)
    calanda_dissolved = st_union(calanda_combined)
    calanda_mask = st_cast(st_make_valid(calanda_dissolved), "POLYGON")

    st_write(
      obj = calanda_mask,
      dsn = here("Calanda_JSDM", "output", "calanda_mask.shp"),
      driver = "ESRI Shapefile",
      append = FALSE
    )
  }

  # Topography ----
  cat("Processing topography...\n")
  if (dir.exists(here("Calanda_JSDM", "output", "metrics"))) {
    topo = extract_topography(
      topo_raster_dir = here("Calanda_JSDM", "output", "metrics"),
      veg_coords_path = here("Calanda_JSDM", "data", "vegetation", "veg.coord.csv")
    )
  } else {
    topo = process_swissalti3d_microtopo(
      csv_file = here("Calanda_JSDM", "data", "ch.swisstopo.swissalti3d-HJOyAIjH.csv"),
      mask_shapefile = NULL,
      output_dir = here("Calanda_JSDM", "output"),
      cut_left = 5,
      cut_right = 18,
      cut_top = 2,
      cut_bottom = 16
    )

    topo = extract_topography(
      topo_raster_dir = here("Calanda_JSDM", "output", "metrics"),
      veg_coords_path = here("Calanda_JSDM", "data", "vegetation", "veg.coord.csv")
    )
  }

  # Assess topography variables via PCA
  print(corrplot(
    cor(as.matrix(topo %>% select(-c(plot_id, releve_id, x, y, dem)) %>% na.omit())),
    method = "number", type = "lower", number.cex = 0.6
  ))
  res_pca = PCA(topo %>% select(c(flowdir, tpi, roughness)) %>% na.omit(), scale.unit = TRUE)

  p1 = fviz_pca_biplot(res_pca,
    repel = TRUE, label = "var",
    col.var = "black",
    col.ind = topo %>% na.omit() %>% pull(dem)
  ) +
    labs(color = "Altitude") +
    scale_color_viridis_c(option = "turbo") +
    theme(legend.position = "bottom", legend.key.size = unit(1, "cm"))
  print(p1)

  # Select TPI, flowdir and roughness
  topo = topo %>%
    mutate(plot_id_releve = paste0(plot_id, releve_id)) %>%
    select(plot_id_releve, flowdir, tpi, roughness)

  # Snow data from Copernicus ----
  cat("Processing snow data...\n")
  if (file.exists(here("Calanda_JSDM", "output", "snow_metrics.csv"))) {
    snow_metrics = read_csv(here("Calanda_JSDM", "output", "snow_metrics.csv"))
  } else {
    snow = process_gfsc_data(
      copernicus_dir = here("Calanda_JSDM", "data", "wekeo"),
      output_dir = here("Calanda_JSDM", "output"),
      veg_coords_path = here("Calanda_JSDM", "data", "vegetation", "veg.coord.csv"),
      calanda_mask_path = here("Calanda_JSDM", "output", "calanda_mask.shp"),
      force_reprocess = FALSE
    )

    snow = read_csv(here("Calanda_JSDM", "output", "gf_vegetation_data.csv")) %>%
      filter(!is.na(raster_value)) %>%
      filter(!(raster_value > 100))

    processed_snow = process_snow_data(
      raw_snow_df = snow,
      expand_time_series = TRUE
    )

    snow_metrics = calculate_snow_metrics(processed_snow)
    write_csv(snow_metrics, here("Calanda_JSDM", "output", "snow_metrics.csv"))
  }

  snow_metrics = snow_metrics %>%
    mutate(plot_id_releve = paste0(plot_id, releve_id)) %>%
    select(plot_id_releve, year, snow_sum, snow_disappearance_doy)

  # Land surface temperature from ECOSTRESS ----
  cat("Processing temperature data...\n")
  lst = bind_rows(
    read_csv(here("Calanda_JSDM", "data", "ecostress",
                  "Calanda-LST-02-2019-ECO-L2T-LSTE-002-results.csv")),
    read_csv(here("Calanda_JSDM", "data", "ecostress",
                  "Calanda-LST-02-2023-ECO-L2T-LSTE-002-results.csv"))
  ) %>%
    filter(ECO_L2T_LSTE_002_QC_Data_quality_flag_Description == "Good quality L1B data") %>%
    select(Category, Latitude, Longitude, Date, ECO_L2T_LSTE_002_LST) %>%
    rename(lst = ECO_L2T_LSTE_002_LST) %>%
    mutate(lst = lst - 273.15) %>%
    distinct() %>%
    rename(plot_id_releve = Category) %>%
    mutate(year = year(Date)) %>%
    filter(!year %in% c(2019, 2025)) %>%
    select(-year)

  processed_temp = process_temperature_data(
    raw_temp_df = lst,
    temp_col = "lst",
    expand_time_series = TRUE
  )

  temp_metrics = calculate_temp_metrics(processed_temp)

  # Freezing degree days between snow disappearance and spring warming ----
  fdd = left_join(temp_metrics, snow_metrics) %>%
    filter(!year %in% c(2019, 2025)) %>%
    mutate(snow_disappearance_doy = ifelse(is.na(snow_disappearance_doy), 1,
                                           snow_disappearance_doy)) %>%
    left_join(processed_temp %>% select(plot_id_releve, year, doy, temp_interpolated)) %>%
    filter(snow_disappearance_doy < spring_warming_doy) %>%
    group_by(plot_id_releve, Latitude, Longitude, year) %>%
    filter(doy <= spring_warming_doy & doy >= snow_disappearance_doy) %>%
    filter(temp_interpolated < 0) %>%
    summarize(fdd = sum(temp_interpolated), .groups = "drop") %>%
    mutate(fdd = ifelse(is.na(fdd), 0, fdd)) %>%
    group_by(plot_id_releve, Latitude, Longitude) %>%
    summarize(fdd = mean(fdd, na.rm = TRUE)) %>%
    ungroup() %>%
    select(plot_id_releve, fdd) %>%
    distinct()

  # Mean across years
  snow_metrics = snow_metrics %>%
    group_by(plot_id_releve) %>%
    summarize(snow_sum = mean(snow_sum, na.rm = TRUE))

  temp_metrics = temp_metrics %>%
    group_by(plot_id_releve) %>%
    summarize(summer_temp = mean(summer_temp, na.rm = TRUE))

  # Evapotranspiration from MODIS ----
  cat("Processing evapotranspiration data...\n")
  ety = bind_rows(
    read_csv(here("Calanda_JSDM", "data", "modis",
                  "Calanda-ET-MODIS-Yearly-2-MOD16A3GF-061-results.csv")),
    read_csv(here("Calanda_JSDM", "data", "modis",
                  "Calanda-MODIS-ET-Yearly-MOD16A3GF-061-results.csv"))
  ) %>%
    mutate(year = year(Date)) %>%
    select(Category, ID, Latitude, Longitude, year, MOD16A3GF_061_ET_500m) %>%
    rename(et_annual = MOD16A3GF_061_ET_500m) %>%
    mutate(et_annual = ifelse(et_annual > 6000, NA, et_annual)) %>%
    rename(plot_id_releve = Category) %>%
    select(-ID) %>%
    group_by(plot_id_releve, Latitude, Longitude) %>%
    summarize(et_annual = mean(et_annual, na.rm = TRUE)) %>%
    mutate(Latitude = round(Latitude, 5),
           Longitude = round(Longitude, 5)) %>%
    ungroup() %>%
    select(plot_id_releve, et_annual) %>%
    distinct() %>%
    filter(!is.nan(et_annual))

  # Combine all environmental data ----
  cat("Combining environmental data...\n")
  veg_clim = veg %>%
    filter(plot_id_releve != "dupl_108.OID2979333") %>%
    select(plot_id_releve, Latitude, Longitude, altitude, slope, southness,
           trees_cover, shrubs_cover, rocks_cover, soil_depth_mean, soil_depth_var) %>%
    mutate(Latitude = round(Latitude, 3),
           Longitude = round(Longitude, 3)) %>%
    distinct() %>%
    left_join(ety) %>%
    left_join(snow_metrics) %>%
    left_join(fdd) %>%
    left_join(temp_metrics) %>%
    left_join(topo) %>%
    mutate(across(everything(), ~if_else(is.nan(.), NA, .))) %>%
    mutate(fdd = ifelse((is.na(fdd) & !is.na(summer_temp)), 0, fdd))

  # Impute missing environmental data ----
  cat("Imputing missing environmental data...\n")
  imp_clim = impute_environmental_data(
    veg_clim,
    variables_to_impute = c("et_annual", "soil_depth_mean", "soil_depth_var"),
    m = 20,
    maxiter = 50,
    num_trees = 500,
    seed = 123,
    validation_fraction = 0.3
  )

  veg_clim %>%
    summarise(across(everything(), ~mean(is.na(.)) * 100)) %>%
    pivot_longer(cols = everything(),
                 names_to = "variable",
                 values_to = "percent_missing") %>%
    arrange(desc(percent_missing)) %>%
    right_join(imp_clim$performance) %>%
    mutate(percent_missing = round(percent_missing),
           r_squared = round(r_squared, 2)) %>%
    select(variable, percent_missing, r_squared) %>%
    gt()

  veg_clim = imp_clim$imputed_data %>%
    select(plot_id_releve, Latitude, Longitude, altitude, slope, summer_temp, fdd,
           et_annual_final, soil_depth_mean_final, soil_depth_var_final,
           trees_cover, shrubs_cover, rocks_cover, flowdir, tpi, roughness, snow_sum) %>%
    rename_with(~str_remove(., "_final$"), ends_with("_final")) %>%
    filter(!is.na(slope)) %>%
    filter(!is.na(altitude))

  write_csv(veg_clim, here("Calanda_JSDM", "output", "veg_clim.csv"))
}
veg_clim = veg_clim %>% select(-c(`...1`))

# ==============================================================================
# SPECIES DATA
# ==============================================================================
cat("\n=== Processing species data ===\n")

veg = veg %>%
  rename(species = taxon_global) %>%
  filter(!is.na(word(species, 2)))

cat("Genus taken out:", length(unique(veg$species)), "species and",
    length(unique(veg$plot_id_releve)), "plots.\n")

veg_comm = veg %>%
  select(plot_id_releve, species_cover, species) %>%
  rename(cover = species_cover) %>%
  distinct() %>%
  filter(!is.na(species)) %>%
  mutate(cover = ifelse(is.na(cover), 0.001, cover)) %>%
  group_by(plot_id_releve, species) %>%
  summarize(cover = sum(cover)) %>%
  group_by(plot_id_releve) %>%
  mutate(total_cover = sum(cover)) %>%
  group_by(plot_id_releve, species) %>%
  mutate(rel_cover = cover / total_cover) %>%
  group_by(plot_id_releve) %>%
  mutate(total_cover = sum(rel_cover)) %>%
  select(-total_cover, -cover) %>%
  ungroup()

veg_pa = veg_comm %>%
  mutate(rel_cover = ifelse(rel_cover > 0, 1, 0)) %>%
  pivot_wider(names_from = species, values_from = rel_cover, values_fill = 0) %>%
  column_to_rownames(var = "plot_id_releve")

veg_rare = veg_pa %>%
  summarise(across(everything(), ~sum(., na.rm = TRUE))) %>%
  pivot_longer(cols = everything()) %>%
  filter(value < nrow(veg_pa) * 0.01) %>%
  pull(name)

veg_pa = veg_comm %>%
  filter(!species %in% veg_rare) %>%
  mutate(rel_cover = ifelse(rel_cover > 0, 1, 0)) %>%
  pivot_wider(names_from = species, values_from = rel_cover, values_fill = 0) %>%
  column_to_rownames(var = "plot_id_releve")

cat("Rare species taken out:", ncol(veg_pa), "species and", nrow(veg_pa), "plots.\n")

veg_abund = veg_comm %>%
  filter(!species %in% veg_rare) %>%
  group_by(plot_id_releve) %>%
  mutate(tot_rel_cover = sum(rel_cover)) %>%
  filter(tot_rel_cover >= 0.8) %>%
  select(-tot_rel_cover) %>%
  pivot_wider(names_from = species, values_from = rel_cover, values_fill = 0)

cat("Plots with less than 80% cover taken out:", ncol(veg_abund) - 1, "species and",
    length(unique(veg_abund$plot_id_releve)), "plots.\n")

veg_abund = veg_abund %>%
  filter(plot_id_releve %in% veg_clim$plot_id_releve) %>%
  column_to_rownames(var = "plot_id_releve")

cat("Plots with missing climate info taken out:", ncol(veg_abund), "species and",
    nrow(veg_abund), "plots.\n")

Y = as.matrix(veg_pa[rownames(veg_abund), ])
cat("Species occurrence counts:\n")
print(as.data.frame(colSums(Y)))

# ==============================================================================
# VEG_TREE: TREE/SHRUB COVER + WOODY SPECIES ABUNDANCES
# ==============================================================================
cat("\n=== Building veg_tree dataset ===\n")

tree_species = c(
  "Abies alba", "Acer campestre", "Acer platanoides", "Acer pseudoplatanus",
  "Fagus sylvatica", "Fraxinus excelsior", "Juglans regia", "Larix decidua",
  "Picea abies", "Pinus sylvestris", "Prunus avium", "Quercus petraea",
  "Salix caprea", "Sorbus aria", "Sorbus aucuparia", "Ulmus glabra"
)

shrub_species = c(
  "Berberis vulgaris", "Clematis vitalba", "Cornus sanguinea", "Corylus avellana",
  "Crataegus monogyna", "Daphne mezereum", "Daphne striata", "Euonymus europaeus",
  "Frangula alnus", "Juniperus communis", "Ligustrum vulgare", "Lonicera xylosteum",
  "Rosa pendulina", "Rubus caesius", "Rubus fruticosus", "Rubus idaeus",
  "Viburnum lantana"
)

# Keep only tree/shrub species present in the cleaned community matrix
cleaned_species = colnames(Y)
tree_species_present = intersect(tree_species, cleaned_species)
shrub_species_present = intersect(shrub_species, cleaned_species)

cat("Tree species in cleaned data:", length(tree_species_present), "\n")
cat(" ", paste(tree_species_present, collapse = ", "), "\n")
cat("Shrub species in cleaned data:", length(shrub_species_present), "\n")
cat(" ", paste(shrub_species_present, collapse = ", "), "\n")

# Tree + shrub abundances from the cleaned community matrix (rel_cover)
woody_species = c(tree_species_present, shrub_species_present)
woody_abund = veg_abund[, intersect(woody_species, colnames(veg_abund)), drop = FALSE]

# Plot-level tree and shrub cover from the original vegetation data
plot_cover = veg %>%
  filter(plot_id_releve %in% rownames(veg_abund)) %>%
  select(plot_id_releve, trees_cover, shrubs_cover) %>%
  distinct() %>%
  group_by(plot_id_releve) %>%
  summarize(trees_cover = mean(trees_cover, na.rm = TRUE),
            shrubs_cover = mean(shrubs_cover, na.rm = TRUE),
            .groups = "drop")

# Combine: plot-level cover + species abundances
veg_tree = woody_abund %>%
  rownames_to_column("plot_id_releve") %>%
  left_join(plot_cover, by = "plot_id_releve") %>%
  relocate(plot_id_releve, trees_cover, shrubs_cover)

# Add growth form label columns for convenience
tree_col_flag = colnames(veg_tree) %in% tree_species_present
shrub_col_flag = colnames(veg_tree) %in% shrub_species_present
growth_form = case_when(
  tree_col_flag ~ "tree",
  shrub_col_flag ~ "shrub",
  TRUE ~ NA_character_
)

cat("veg_tree:", nrow(veg_tree), "plots,",
    length(tree_species_present), "tree species,",
    length(shrub_species_present), "shrub species\n")

write_csv(veg_tree, here("Calanda_JSDM", "output", "veg_tree.csv"))
cat("Saved output/veg_tree.csv\n")

# ==============================================================================
# TRAIT DATA
# ==============================================================================
cat("\n=== Loading trait data ===\n")
traits = read_csv(here("Calanda_JSDM", "output", "traits_medians_imputed.csv"), show_col_types = FALSE) %>%
  rename(species = species_TNRS)
cat("Loaded traits for", nrow(traits), "species\n")

# Load imputation flags
means_flags = read_csv(here("Calanda_JSDM", "output", "traits_medians_imputation_flags.csv"),
                        show_col_types = FALSE) %>%
  rename(species = species_TNRS)

kcv_flags = read_csv(here("Calanda_JSDM", "output", "traits_kcv_imputation_flags.csv"),
                      show_col_types = FALSE) %>%
  rename(species = species_TNRS)

# --- Species imputation summary table ---
# Number of sites each species occurs in (from Y matrix)
species_sites = tibble(
  species = colnames(Y),
  n_sites = colSums(Y)
)

# Pivot means flags to long format
means_flags_long = means_flags %>%
  pivot_longer(-species, names_to = "trait", values_to = "status") %>%
  mutate(type = "Median")

# Pivot kCV flags to long format
kcv_flags_long = kcv_flags %>%
  pivot_longer(-species, names_to = "trait", values_to = "status") %>%
  mutate(
    trait = str_remove(trait, "^kCV_"),
    type = "kCV"
  )

# For median flags, clean trait names (remove Median_ prefix for display)
means_flags_long = means_flags_long %>%
  mutate(trait = str_remove(trait, "^Median_"))

# Combine and pivot wide: one column per trait showing Original/Imputed
imputation_summary = bind_rows(means_flags_long, kcv_flags_long) %>%
  filter(species %in% species_sites$species) %>%
  mutate(label = paste0(type, "_", trait)) %>%
  select(species, label, status) %>%
  pivot_wider(names_from = label, values_from = status) %>%
  left_join(species_sites, by = "species") %>%
  arrange(desc(n_sites)) %>%
  select(species, n_sites, everything())

# Count imputed traits per species, separately for medians and kCV
median_cols = grep("^Median_", colnames(imputation_summary), value = TRUE)
kcv_cols = grep("^kCV_", colnames(imputation_summary), value = TRUE)

imputation_summary = imputation_summary %>%
  mutate(
    n_imputed_median = rowSums(across(all_of(median_cols), ~ . == "Imputed"), na.rm = TRUE),
    n_imputed_kcv = rowSums(across(all_of(kcv_cols), ~ . == "Imputed"), na.rm = TRUE)
  )

cat("\n=== Species imputation summary ===\n")
cat(sprintf("Species in community matrix: %d\n", nrow(species_sites)))
cat(sprintf("Species with imputation flags: %d\n",
            sum(species_sites$species %in% means_flags$species)))

# Build gt table with per-cell conditional coloring
status_cols = setdiff(colnames(imputation_summary),
                      c("species", "n_sites", "n_imputed_median", "n_imputed_kcv"))

imputation_gt = imputation_summary %>%
  gt() %>%
  tab_header(
    title = "Species trait imputation summary",
    subtitle = "Original = observed, Imputed = filled by random forest, NA = excluded"
  ) %>%
  cols_label(species = "Species", n_sites = "N sites",
             n_imputed_median = "N imputed (median)",
             n_imputed_kcv = "N imputed (kCV)")

# Apply per-cell coloring: green for Original, orange for Imputed
for (col in status_cols) {
  imputation_gt = imputation_gt %>%
    tab_style(
      style = cell_fill(color = "#D4EDDA"),
      locations = cells_body(columns = !!sym(col),
                             rows = .data[[col]] == "Original")
    ) %>%
    tab_style(
      style = cell_fill(color = "#FDDBC7"),
      locations = cells_body(columns = !!sym(col),
                             rows = .data[[col]] == "Imputed")
    )
}

print(imputation_gt)

# Save as CSV for reference
write_csv(imputation_summary,
          here("Calanda_JSDM", "output", "species_imputation_summary.csv"))
cat("Saved output/species_imputation_summary.csv\n")

# ==============================================================================
# COMMUNITY TRAIT CALCULATIONS (with and without imputed values)
# ==============================================================================
cat("\n=== Calculating community traits ===\n")

cwm_trait_cols = c("Nutrients", "disturbance", "dispersal", "Median_LDMC", "Median_LMA", "Median_LNC", "Median_LCC", "Median_vegetative_height",
                   "Median_seed_mass", "Median_LA")
cwm_log_traits = c("Median_LDMC", "Median_LMA", "Median_LNC", "Median_LCC", "Median_vegetative_height","Median_seed_mass", "Median_LA")

# Version 1: All traits (including imputed)
community_traits = calculate_community_traits(
  community_data = veg_abund,
  traits_data = traits,
  species_col = "species",
  abundance_col = "rel_cover",
  trait_cols = cwm_trait_cols,
  log_traits = cwm_log_traits
)

# Version 2: Only original (non-imputed) trait values
# Set imputed trait values back to NA using the imputation flags
traits_original = traits
for (col in cwm_trait_cols) {
  flag_col = col
  if (flag_col %in% names(means_flags) && col %in% names(traits_original)) {
    imputed_species = means_flags$species[means_flags[[flag_col]] == "Imputed"]
    traits_original[[col]][traits_original$species %in% imputed_species] = NA
  }
}

community_traits_original = calculate_community_traits(
  community_data = veg_abund,
  traits_data = traits_original,
  species_col = "species",
  abundance_col = "rel_cover",
  trait_cols = cwm_trait_cols,
  log_traits = cwm_log_traits
)

# Compare the two versions
cat("\nCommunity traits with all (imputed) values:\n")
cat("  Plots:", nrow(community_traits), "\n")
cat("  NaN in CWMs:", sum(is.nan(as.matrix(community_traits %>% select(ends_with("_cwm")))), na.rm = TRUE), "\n")
cat("  Mean trait coverage (% of abundance with data):\n")
community_traits %>%
  summarise(across(ends_with("_coverage"), ~ round(mean(.x, na.rm = TRUE) * 100, 1))) %>%
  pivot_longer(everything(), names_to = "trait", values_to = "mean_coverage_pct") %>%
  mutate(trait = str_remove(trait, "_coverage$")) %>%
  print(n = Inf)

cat("\nCommunity traits with original values only:\n")
cat("  Plots:", nrow(community_traits_original), "\n")
cat("  NaN in CWMs:", sum(is.nan(as.matrix(community_traits_original %>% select(ends_with("_cwm")))), na.rm = TRUE), "\n")
cat("  Mean trait coverage (% of abundance with data):\n")

community_traits_original %>%
  summarise(across(ends_with("_coverage"), ~ round(mean(.x, na.rm = TRUE) * 100, 1))) %>%
  pivot_longer(everything(), names_to = "trait", values_to = "mean_coverage_pct") %>%
  mutate(trait = str_remove(trait, "_coverage$")) %>%
  print(n = Inf)

# Compare imputed vs original CWMs: scatter plots colored by coverage ----
cat("\n=== Plotting CWM comparison (imputed vs original) ===\n")

cwm_imputed_long = community_traits %>%
  select(plot_id_releve, ends_with("_cwm")) %>%
  pivot_longer(-plot_id_releve, names_to = "trait", values_to = "cwm_imputed") %>%
  mutate(trait = str_remove(trait, "_cwm$"))

cwm_original_long = community_traits_original %>%
  select(plot_id_releve, ends_with("_cwm")) %>%
  pivot_longer(-plot_id_releve, names_to = "trait", values_to = "cwm_original") %>%
  mutate(trait = str_remove(trait, "_cwm$"))

# Coverage from the original-only version (shows how much data was truly observed)
coverage_long = community_traits_original %>%
  select(plot_id_releve, ends_with("_coverage")) %>%
  pivot_longer(-plot_id_releve, names_to = "trait", values_to = "coverage") %>%
  mutate(trait = str_remove(trait, "_coverage$"))

cwm_comparison = cwm_imputed_long %>%
  left_join(cwm_original_long, by = c("plot_id_releve", "trait")) %>%
  left_join(coverage_long, by = c("plot_id_releve", "trait")) %>%
  filter(!is.nan(cwm_imputed) & !is.nan(cwm_original))

# Compute per-trait correlation labels
cor_labels = cwm_comparison %>%
  group_by(trait) %>%
  summarise(
    r_pearson = cor(cwm_original, cwm_imputed, use = "complete.obs", method = "pearson"),
    r_spearman = cor(cwm_original, cwm_imputed, use = "complete.obs", method = "spearman"),
    label = sprintf("Pearson: %.3f\nSpearman: %.3f", r_pearson, r_spearman),
    .groups = "drop"
  )
# Vegetative height comparisons have r = 1 because there are only a couple of imputed species in the final selected communities.
p_cwm_compare = ggplot(cwm_comparison, aes(x = cwm_original, y = cwm_imputed, color = coverage)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text(data = cor_labels, aes(x = -Inf, y = Inf, label = label),
            hjust = -0.1, vjust = 1.3, size = 3.2, color = "black", inherit.aes = FALSE) +
  scale_color_viridis_c(option = "mako", direction = -1,
                         labels = scales::percent, name = "Trait coverage\n(original)") +
  facet_wrap(~ trait, scales = "free") +
  labs(x = "CWM (original only)",
       y = "CWM (with imputed)",
       title = "Community-weighted means: original vs imputed trait values",
       subtitle = "Points colored by proportion of community abundance with observed trait data") +
  theme_bw() +
  theme(text = element_text(size = 11),
        strip.text = element_text(size = 11, face = "bold"),
        legend.position = "right")

pdf(here("Calanda_JSDM", "plot", "cwm_imputed_vs_original.pdf"), height = 8, width = 10)
print(p_cwm_compare)
dev.off()
cat("Saved plot/cwm_imputed_vs_original.pdf\n")

# Save CWM comparison and community traits
write_csv(cwm_comparison, here("Calanda_JSDM", "output", "cwm_comparison.csv"))
write_csv(community_traits, here("Calanda_JSDM", "output", "community_traits_imputed.csv"))
write_csv(community_traits_original, here("Calanda_JSDM", "output", "community_traits_original.csv"))
cat("Saved output/cwm_comparison.csv\n")
cat("Saved output/community_traits_imputed.csv\n")
cat("Saved output/community_traits_original.csv\n")

# ==============================================================================
# ADD CWM NUTRIENTS AND DISTURBANCE TO ENVIRONMENTAL DATA
# ==============================================================================
cat("\n=== Adding CWM nutrients and disturbance to environmental data ===\n")

cwm_scores = community_traits %>%
  select(plot_id_releve, Nutrients_cwm, disturbance_cwm)%>%
  rename(nutrient = Nutrients_cwm, disturbance = disturbance_cwm)

veg_clim = veg_clim %>%
  left_join(cwm_scores, by = "plot_id_releve") %>%
  na.omit()

# ==============================================================================
# FINAL ENVIRONMENTAL MATRIX
# ==============================================================================
cat("\n=== Building final matrices ===\n")

veg_env = veg_clim %>%
  group_by(plot_id_releve, Longitude, Latitude, altitude) %>%
  summarize(across(slope:disturbance, ~mean(.)), .groups = "drop") %>%
  column_to_rownames(var = "plot_id_releve") %>%
  na.omit()%>%
  mutate(across(Longitude:disturbance, ~as.numeric(scale(.))))

# PCA of environmental variables
res_pca_env = PCA(veg_env %>% select(-c(Latitude, Longitude, altitude)))

# Environmental data diagnostic plots
pdf(here("Calanda_JSDM", "plot", "climate_calanda_surveys.pdf"), height = 10, width = 15)
print(
  veg_clim %>%
    filter(plot_id_releve %in% rownames(veg_env)) %>%
    select(-c(plot_id_releve)) %>%
    pivot_longer(cols = everything()) %>%
    group_by(name) %>%
    mutate(median_value = median(value, na.rm = TRUE)) %>%
    ggplot(aes(value)) +
    facet_wrap(. ~ name, scales = "free") +
    theme_bw() +
    geom_histogram(fill = "skyblue") +
    geom_vline(aes(xintercept = median_value), color = "red",
               linetype = "dashed", linewidth = 0.8) +
    labs(x = "Values across all plots")
)

altitude_color = veg_clim %>%
  filter(plot_id_releve %in% rownames(veg_env)) %>%
  distinct(plot_id_releve, .keep_all = TRUE) %>%
  arrange(match(plot_id_releve, rownames(veg_env))) %>%
  pull(altitude)

p1 = fviz_pca_biplot(res_pca_env,
  repel = TRUE, label = "var",
  col.var = "black",
  col.ind = altitude_color
) +
  labs(color = "Altitude") +
  scale_color_viridis_c(option = "turbo") +
  theme(legend.position = "bottom", legend.key.size = unit(1, "cm"))

p2 = fviz_pca_biplot(res_pca_env,
  axes = c(2, 3),
  repel = TRUE, label = "var",
  col.var = "black",
  col.ind = altitude_color
) +
  labs(color = "Altitude") +
  scale_color_viridis_c(option = "turbo") +
  theme(legend.position = "bottom", legend.key.size = unit(1, "cm"))

print(p1 + p2)

print(corrplot(
  cor(as.matrix(veg_clim %>%
    filter(plot_id_releve %in% rownames(veg_env)) %>%
    select(-c(plot_id_releve, Latitude, Longitude)))),
  method = "number", type = "lower", number.cex = 0.6
))

dev.off()

# Final X and Y matrices ----
X = veg_env
plots = intersect(rownames(Y), rownames(X))
Y = Y[plots, ]
X = X[plots, ]
cat("Final:", ncol(veg_abund), "species and", nrow(X), "plots.\n")

data_calanda_jsdm = list()
data_calanda_jsdm$X = X
data_calanda_jsdm$Y = Y

save.image(here("Calanda_JSDM", "output", paste0("starter_data_",Sys.Date(),".RData")))
saveRDS(data_calanda_jsdm, file = here("Calanda_JSDM", "output", paste0("data_calanda_jsdm_",Sys.Date(),".rds")))
cat("\n=== Data preparation complete ===\n")
