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
#   - data/traits/try.quantitative_traits_2025-03-17.csv
#   - data/traits/indicators_cleaned_calanda_2025-03-17.csv
#   - data/traits/dispersal_cleaned_calanda_2025-03-17.csv
#
# Outputs:
#   - output/veg_clim.csv
#   - output/traits.csv
#   - output/starter_data_25.04.25.RData
#   - plot/land_use_variable.pdf
#   - plot/climate_calanda_surveys.pdf
#
# Requires: R/00_setup/functions_calanda.R
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
  mutate(taxon_global = stri_trans_general(taxon_global, "Latin-ASCII"))

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
# TRAIT DATA
# ==============================================================================
cat("\n=== Processing trait data ===\n")

if (file.exists(here("Calanda_JSDM", "output", "traits.csv"))) {
  traits = read_csv(here("Calanda_JSDM", "output", "traits.csv"))
} else {
  try_data = read_csv(here("Calanda_JSDM", "data", "traits",
                           "try.quantitative_traits_2025-03-17.csv"))
  indicators = read_csv(here("Calanda_JSDM", "data", "traits",
                             "indicators_cleaned_calanda_2025-03-17.csv"))
  dispersal = read_csv(here("Calanda_JSDM", "data", "traits",
                            "dispersal_cleaned_calanda_2025-03-17.csv"))

  traits = try_data %>%
    filter(Trait %in% c("N_percent", "LDMC", "vegetative_height", "SLA", "LA",
                        "C_percent", "seed_mass")) %>%
    mutate(species = ifelse(species == "Vaccinium uliginosum subsp. uliginosum",
                            "Vaccinium uliginosum", species)) %>%
    filter(species %in% colnames(Y)) %>%
    filter(!is.na(Value), !is.na(Trait)) %>%
    filter(!Climate_code %in% c("Af", "Am", "Aw", "BWh", "BWk", "BSh", "BSk")) %>%
    group_by(species, species_TNRS, Trait) %>%
    summarize(Mean = mean(Value, na.rm = TRUE),
              Var = var(Value, na.rm = TRUE),
              .groups = "drop") %>%
    pivot_wider(names_from = Trait,
                values_from = c(Mean, Var),
                names_glue = "{.value}_{Trait}") %>%
    ungroup()

  traits = left_join(dispersal, traits) %>%
    left_join(indicators) %>%
    select(species, species_TNRS,
           Mean_seed_mass, Var_seed_mass, dispersal_distance_class,
           Mean_LA, Var_LA, Mean_SLA, Var_SLA, Mean_LDMC, Var_LDMC,
           Mean_vegetative_height, Var_vegetative_height,
           Mean_N_percent, Var_N_percent, Mean_C_percent, Var_C_percent,
           Light, Moisture, Nutrients, Disturbance.Severity) %>%
    rename(disturbance = Disturbance.Severity,
           dispersal = dispersal_distance_class)

  imp_traits = impute_functional_traits(
    traits,
    variables_to_impute = c("Mean_seed_mass", "dispersal", "Mean_LA", "Mean_SLA",
                            "Mean_LDMC", "Mean_vegetative_height", "Mean_N_percent",
                            "Mean_C_percent", "Light", "Moisture", "Nutrients",
                            "disturbance"),
    m = 30,
    maxiter = 100,
    num_trees = 500,
    seed = 123,
    validation_fraction = 0.3
  )

  traits %>%
    summarise(across(everything(), ~mean(is.na(.)) * 100)) %>%
    pivot_longer(cols = everything(),
                 names_to = "variable",
                 values_to = "percent_missing") %>%
    arrange(desc(percent_missing)) %>%
    right_join(imp_traits$performance) %>%
    mutate(percent_missing = round(percent_missing),
           r_squared = round(r_squared, 2)) %>%
    select(variable, percent_missing, r_squared) %>%
    gt()

  traits = imp_traits$imputed_data %>%
    select(species, species_TNRS, starts_with("Var_"), contains("_final")) %>%
    rename_with(~str_remove(., "_final$"), ends_with("_final"))
  write_csv(traits, here("Calanda_JSDM", "output", "traits.csv"))
}

# ==============================================================================
# LAND-USE PCA
# ==============================================================================
cat("\n=== Creating land-use proxy via PCA ===\n")

community_traits = calculate_community_traits(
  community_data = veg_abund,
  traits_data = traits,
  species_col = "species",
  abundance_col = "rel_cover",
  trait_cols = c("Mean_LDMC", "Nutrients", "Mean_SLA")
)

res_pca_lu = left_join(
  community_traits %>% select(plot_id_releve, Nutrients_cwm, Mean_LDMC_cwm, Mean_SLA_cwm),
  veg_abund %>% rownames_to_column(var = "plot_id_releve") %>%
    select(plot_id_releve, `Nardus stricta`)
) %>%
  column_to_rownames(var = "plot_id_releve")

res_pca_lu = PCA(res_pca_lu, scale.unit = TRUE)

pdf(here("Calanda_JSDM", "plot", "land_use_variable.pdf"), height = 8, width = 8)
fviz_pca_biplot(res_pca_lu,
  repel = TRUE, label = "var",
  col.var = "black", col.ind = "grey70"
) +
  theme(legend.position = "bottom", legend.key.size = unit(1, "cm"))
dev.off()

land_use_scores = data.frame(
  land_use = res_pca_lu$ind$coord[, "Dim.1", drop = TRUE]
) %>%
  rownames_to_column(var = "plot_id_releve")

veg_clim = veg_clim %>%
  left_join(land_use_scores) %>%
  na.omit()

# ==============================================================================
# FINAL ENVIRONMENTAL MATRIX
# ==============================================================================
cat("\n=== Building final matrices ===\n")

veg_env = veg_clim %>%
  group_by(plot_id_releve, Longitude, Latitude, altitude) %>%
  summarize(across(slope:land_use, ~mean(.))) %>%
  ungroup() %>%
  column_to_rownames(var = "plot_id_releve") %>%
  na.omit() %>%
  mutate(across(Longitude:land_use, ~as.numeric(scale(.))))

# PCA of environmental variables
res_pca_env = PCA(
  veg_clim %>%
    filter(plot_id_releve %in% rownames(veg_env)) %>%
    select(-c(plot_id_releve, Latitude, Longitude, altitude, trees_cover,
              snow_sum, shrubs_cover)),
  scale.unit = TRUE
)

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

p1 = fviz_pca_biplot(res_pca_env,
  repel = TRUE, label = "var",
  col.var = "black",
  col.ind = veg_clim %>%
    filter(plot_id_releve %in% rownames(veg_env)) %>%
    pull(altitude)
) +
  labs(color = "Altitude") +
  scale_color_viridis_c(option = "turbo") +
  theme(legend.position = "bottom", legend.key.size = unit(1, "cm"))

p2 = fviz_pca_biplot(res_pca_env,
  axes = c(2, 3),
  repel = TRUE, label = "var",
  col.var = "black",
  col.ind = veg_clim %>%
    filter(plot_id_releve %in% rownames(veg_env)) %>%
    pull(altitude)
) +
  labs(color = "Altitude") +
  scale_color_viridis_c(option = "turbo") +
  theme(legend.position = "bottom", legend.key.size = unit(1, "cm"))

print(p1 + p2)

print(corrplot(
  cor(as.matrix(veg_clim %>%
    filter(plot_id_releve %in% rownames(veg_env)) %>%
    select(-c(plot_id_releve, Latitude, Longitude, trees_cover)))),
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

# save.image(here("Calanda_JSDM", "output", "starter_data_25.04.25.RData"))

cat("\n=== Data preparation complete ===\n")
