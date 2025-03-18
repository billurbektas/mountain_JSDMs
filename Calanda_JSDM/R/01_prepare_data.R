# Get Calanda limits
if(file.exists("output/calanda_mask.shp")){
  calanda_mask = st_read("output/calanda_mask.shp")
}else{
calanda_1 = st_read("data/mask/study_region_2024.shp")
calanda_1 = st_transform(calanda_1, crs = "EPSG:4326")

calanda_2 = st_read("data/mask/study_region.shp")
calanda_2 = st_transform(calanda_2, crs = "EPSG:4326")
common_cols = intersect(names(calanda_1), names(calanda_2))
calanda_1 = calanda_1[, common_cols]
calanda_2 = calanda_2[, common_cols]
calanda_combined = rbind(calanda_1, calanda_2)
calanda_dissolved = st_union(calanda_combined)
calanda_mask = st_cast(st_make_valid(calanda_dissolved), "POLYGON")

st_write(
  obj = calanda_mask,                    
  dsn = "output/calanda_mask.shp",   
  driver = "ESRI Shapefile", 
  append = FALSE
)
}
#plot(st_geometry(calanda_mask), border = "black", lwd = 2, main = "Merged Polygon without Internal Lines")

# Process snow data from Copernicus ----

if(file.exists("output/snow_metrics.csv")){
  snow_metrics = read_csv("output/snow_metrics.csv")[,-1]
}else{
snow = 
  process_gfsc_data(
  copernicus_dir = "data/wekeo",
  output_dir = "output",
  veg_coords_path = "data/vegetation/veg.coord.csv",
  calanda_mask_path = "output/calanda_mask.shp",
  force_reprocess = FALSE
)

snow = read_csv("output/gf_vegetation_data.csv")%>%
  filter(!is.na(raster_value))%>%
  filter(!(raster_value >100))

processed_snow = process_snow_data(
  raw_snow_df = snow,
  expand_time_series = TRUE
)
  
snow_metrics = calculate_snow_metrics(processed_snow) 
write.csv(snow_metrics, file = "output/snow_metrics.csv")
}

snow_metrics =
  snow_metrics %>%
  mutate(plot_id_releve = paste0(plot_id, releve_id))%>%
  select(plot_id_releve, year, snow_sum, snow_disappearance_doy)

# Get land surface temperature from ECOSTRESS
## The slight problem with this data is that it is not "daily average" because they are measured only once.
## However, all the points in our site are measured at the same hour of the day so there should be no bias.
### We calculate: summer temperature + freezing degree days after snow melt until summer?
lst = bind_rows(read_csv("data/ecostress/Calanda-LST-02-2019-ECO-L2T-LSTE-002-results.csv"),
                read_csv("data/ecostress/Calanda-LST-02-2023-ECO-L2T-LSTE-002-results.csv")
                )%>%
  filter(ECO_L2T_LSTE_002_QC_Data_quality_flag_Description=="Good quality L1B data")%>%
  select(Category, Latitude, Longitude, Date, ECO_L2T_LSTE_002_LST)%>%
  rename(lst = ECO_L2T_LSTE_002_LST)%>%
  mutate(lst = lst -273.15)%>% # Kelvin to Celcius
  distinct()%>%
  rename(plot_id_releve = Category)%>%
  mutate(year = year(Date))%>%
  filter(!year %in% c(2019, 2025))%>%
  select(-year)

processed_temp =
  process_temperature_data(
  raw_temp_df = lst,
  temp_col = "lst",  # Column with temperature values
  expand_time_series = TRUE
)

temp_metrics = calculate_temp_metrics(processed_temp)

ggplot(processed_temp %>% filter(plot_id_releve == "23.OID2979726") , aes(doy, temp_interpolated))+
  geom_line()+
  geom_line(data = processed_temp %>% filter(plot_id_releve == "23.OID2979726"), aes(doy, temp_raw),
            color = "grey50")+
  geom_vline(xintercept = c(151, 95))

# Calculate fdd between snow disseappearance and the first day when the temperatures = 1°C
## if the first day temperature == 1°C is earlier than snow diseappearance than by default fdd = 0
fdd = 
  left_join(temp_metrics, snow_metrics)%>%
  filter(!year %in% c(2019, 2025))%>%
  mutate(snow_disappearance_doy = ifelse(is.na(snow_disappearance_doy), 1, snow_disappearance_doy))%>% # no snow
  left_join(processed_temp %>% select(plot_id_releve, year, doy, temp_interpolated))%>%
  filter(snow_disappearance_doy < spring_warming_doy) %>%
  group_by(plot_id_releve, Latitude, Longitude, year) %>%
  filter(doy <= spring_warming_doy & doy >= snow_disappearance_doy) %>%
  filter(temp_interpolated < 0) %>%
  summarize(fdd = sum(temp_interpolated), .groups = "drop") %>%
  mutate(fdd = ifelse(is.na(fdd), 0, fdd))%>%
  group_by(plot_id_releve, Latitude, Longitude)%>%
  summarize(fdd = mean(fdd, na.rm = TRUE))%>%
  ungroup()%>%
  select(plot_id_releve, fdd)%>%
  distinct()

# Get mean across years
snow_metrics = 
  snow_metrics %>%
  group_by(plot_id_releve)%>%
  summarize(snow_sum = mean(snow_sum, na.rm = TRUE))

temp_metrics = 
  temp_metrics %>%
  group_by(plot_id_releve)%>%
  summarize(summer_temp = mean(summer_temp, na.rm = TRUE))

# Get the evapotranspiration from the MODIS data 
### Summer evapotranspiration ~ summer drought
### the summer ET and annual ET are very correlated. I am just taking the annual ET.
ety =  bind_rows(read_csv("data/modis/Calanda-ET-MODIS-Yearly-2-MOD16A3GF-061-results.csv"),
                read_csv("data/modis/Calanda-MODIS-ET-Yearly-MOD16A3GF-061-results.csv"))%>%
  mutate(year = year(Date))%>%
  select(Category, ID, Latitude, Longitude, year, MOD16A3GF_061_ET_500m)%>%
  rename(et.annual = MOD16A3GF_061_ET_500m)%>%
  mutate(et.annual = ifelse(et.annual>6000, NA, et.annual))%>%
  rename(plot_id_releve = Category)%>%
  select(-ID)%>%
  #replace_with_closest_site_data(.)%>%
  group_by(plot_id_releve, Latitude, Longitude)%>%
  summarize(et.annual = mean(et.annual, na.rm = TRUE))%>%
  mutate(Latitude = round(Latitude, 5),
         Longitude = round(Longitude, 5)) %>%
  ungroup()%>%
  select(plot_id_releve, et.annual)%>%
  distinct()%>%
  filter(!is.nan(et.annual))

# Get land-use intensity 
# First run - processes and saves the raster
land_use = extract_mowing_events(
  raster_pattern = "grassland-use_intensity_.*\\.tif$",
  coords_path = NULL,
  raster_dir = "data/land_use",
  processed_raster_path = "output/land_use.tif"
)

# Use the cached raster to extract points
land_use = extract_mowing_events(
  raster_pattern = "grassland-use_intensity_.*\\.tif$",
  coords_path = "data/vegetation/veg.coord.csv",
  output_path = "output/mowing_events_at_veg_points.csv",
  raster_dir = "data/land_use",
  processed_raster_path = "output/land_use.tif"
)

land_use = 
  land_use %>%
  mutate(plot_id_releve = paste0(plot_id, releve_id))%>%
  rename(Latitude = y, Longitude = x)%>%
  pivot_longer(cols = mowingEvents_2018:mowingEvents_2020)%>%
  group_by(plot_id_releve, Latitude, Longitude)%>%
  summarize(land_use = round(mean(value, na.rm = TRUE)))%>%
  mutate(land_use = ifelse(is.nan(land_use), NA, land_use))%>%
  ungroup()%>%
  select(plot_id_releve, land_use)%>%
  distinct()
 
# Get climate data and reproject
# clim = rast(list.files("data/climate", full.names = TRUE))
# crs(clim) = "EPSG:2056"
# clim = project(clim, "EPSG:4326")
# cropped_clim = crop(clim, ext(vect(calanda_mask)))
# cropped_clim_df = 
#   as.data.frame(cropped_clim, xy = TRUE, na.rm = TRUE)%>%
#   pivot_longer(cols = AI_8110_LV95:gdd5Y_8110_LV95)%>%
#   mutate(name = str_extract(name, "^[^_]+"))

#Get sp list for transplant Calanda
tp = CH_Calanda$community

# Step 1: Create the full list of species names
tp.species = tp %>% distinct(SpeciesName) %>% pull(SpeciesName)

if(file.exists("data/vegetation/tp.species.csv")){
  tp.species = read_csv(file = "data/vegetation/tp.species.csv")
}else{
  tp.species = TNRS(taxonomic_names = tp.species)
  write.csv(tp.species, file = "data/vegetation/tp.species.csv")
}  

tp.species =
  tp.species %>%
  dplyr::select(Name_submitted, Accepted_name) %>%
  rename(SpeciesName = Name_submitted, species = Accepted_name)

tp = left_join(tp, tp.species)

tp.species = tp %>% distinct(SpeciesName)

# Get vegetation surveys in Calanda

veg = read_csv("data/vegetation/2024_CAPHE_SpeDis_CleanData_20240214.csv") %>%
  mutate(plot_id_releve = paste0(plot_id, releve_id)) %>%
  mutate(across(starts_with("soil_depth"), ~ as.numeric(.))) %>% # Ensure all soil_depth columns are numeric
  rowwise() %>% 
  mutate(soil_depth_mean = mean(c_across(starts_with("soil_depth")), na.rm = TRUE)) %>%
  mutate(soil_depth_var = var(c_across(starts_with("soil_depth")), na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(southness = abs(aspect - 180)) %>%
  mutate(altitude = rowMeans(select(., altitude_max, altitude_min), na.rm = TRUE)) %>%
  dplyr::select(-altitude_min, -altitude_max) %>%
  dplyr::select(plot_id_releve, x, y, slope, southness, trees_cover, shrubs_cover, rocks_cover, soil_depth_mean, soil_depth_var, open_close, taxon_global, species_cover, altitude) %>%
  distinct() %>%
  rename(Latitude = y, Longitude = x)%>%
  mutate(taxon_global = stri_trans_general(taxon_global, "Latin-ASCII"))

# veg.clim = 
#   veg %>%
#   dplyr::select(plot_releve_id, x, y)%>%
#   vect(., geom = c("x", "y"), crs = "EPSG:4326")%>%
#   terra::extract(clim, ., xy = FALSE, ID = FALSE)%>%
#   rename_with(~ str_extract(., "^[^_]+"))%>%
#   bind_cols(veg %>% dplyr::select(plot_releve_id, altitude, slope, southness, trees_cover, shrubs_cover, rocks_cover, soil_depth_mean, soil_depth_var))
# 
# veg.climx =
#   veg.clim %>% 
#   distinct()%>%
#   na.omit()

veg.clim = 
  veg %>% 
  filter(plot_id_releve != "dupl_108.OID2979333")%>%
  select(plot_id_releve, Latitude, Longitude, altitude, slope, southness, trees_cover, shrubs_cover, rocks_cover, soil_depth_mean, soil_depth_var)%>%
  mutate(Latitude = round(Latitude, 3),
         Longitude = round(Longitude, 3))%>%
  distinct()%>%
  left_join(ety)%>%
  left_join(snow_metrics)%>%
  left_join(fdd)%>%
  left_join(temp_metrics)%>%
  left_join(land_use)%>%
  mutate(across(everything(), ~if_else(is.nan(.), NA, .)))%>%
  mutate(fdd = ifelse((is.na(fdd) & !is.na(summer_temp)), 0, fdd ))

imp.clim = 
impute_environmental_data(
  veg.clim,
  variables_to_impute = c("et.annual", "soil_depth_mean", "soil_depth_var", "land_use"),
  m = 20,
  maxiter = 50,
  num.trees = 500,
  seed = 123,
  validation_fraction = 0.3
)

veg.clim %>%
  summarise(across(everything(), ~mean(is.na(.)) * 100)) %>%
  pivot_longer(cols = everything(), 
               names_to = "variable", 
               values_to = "percent_missing") %>%
  arrange(desc(percent_missing))%>%
  right_join(imp.clim$performance) %>%
  mutate(percent_missing = round(percent_missing),
         r_squared = round(r_squared, 2))%>%
  select(variable, percent_missing, r_squared)%>%
  gt()

veg.clim = 
  imp.clim$imputed_data %>%
  select(plot_id_releve, Latitude, Longitude, altitude, slope, summer_temp, fdd, et.annual_final,
         soil_depth_mean_final, soil_depth_var_final, land_use_final, trees_cover, shrubs_cover, rocks_cover,
         snow_sum
         )%>%
  rename_with(~str_remove(., "_final$"), ends_with("_final"))%>%
  filter(!is.na(slope))%>%
  filter(!is.na(altitude))
  
write.csv(veg.clim, file = "output/veg.clim.csv")

res.pca = PCA(veg.clim %>% dplyr::select(-c(plot_id_releve, Latitude, Longitude, altitude)), scale.unit = TRUE)

# Plot using ggplot2
pdf("plot/climate_calanda_surveys.pdf", height = 10, width = 15)
print(veg.clim %>% dplyr::select(-c(plot_id_releve, Latitude, Longitude))%>%
  pivot_longer(cols = everything()) %>%
  group_by(name) %>%
  mutate(median_value = median(value, na.rm = TRUE)) %>%
  ggplot(aes(value)) +
  facet_wrap(. ~ name, scales = "free") +
  theme_minimal() +
  geom_histogram(fill = "skyblue") +
  geom_vline(aes(xintercept = median_value), color = "red", linetype = "dashed", size = 0.8) +
  labs(x = "Values across all plots (611 plots - when NAs omitted)"))

p1 = fviz_pca_biplot(res.pca, 
                repel = TRUE,
                label = "var",
                col.var = "black", # Variables color
                col.ind = veg.clim$altitude)+
  labs(color = "Altitude")+
  scale_color_viridis_c(option = "turbo")+
  theme(legend.position = "bottom",
        legend.key.size = unit(1, "cm"))

p2 = fviz_pca_biplot(res.pca, 
                axes = c(2,3),
                repel = TRUE,
                label = "var",
                col.var = "black", # Variables color
                col.ind = veg.clim$altitude)+
  labs(color = "Altitude")+
  scale_color_viridis_c(option = "turbo")+
  theme(legend.position = "bottom",
        legend.key.size = unit(1, "cm"))

print(p1+p2)

print(corrplot(cor(as.matrix(veg.clim %>% dplyr::select(-c(plot_id_releve, Latitude, Longitude)))), method = "number", type = "lower", number.cex = 0.6))

dev.off()

# sp.list = veg %>% pull(taxon_global) %>% unique()

# if(file.exists("data/vegetation/sp.list.csv")){
#   sp.list = read_csv(file = "data/vegetation/sp.list.csv")
# }else{
#   sp.list = TNRS(taxonomic_names = sp.list)
#   write.csv(sp.list, file = "data/vegetation/sp.list.csv")
# }  

# sp.list =
#   sp.list %>%
#   dplyr::select(Name_submitted, Accepted_name) %>%
#   rename(taxon_name = Name_submitted, species = Accepted_name)
# 
# veg = left_join(veg, sp.list)
# 
veg = 
  veg %>%
  rename(species = taxon_global)%>%
  filter(!is.na(word(species, 2))) #Take out genus

veg.env = 
  veg.clim %>% 
  column_to_rownames(var = "plot_id_releve")%>%
  na.omit(.)%>%
  mutate(across(Latitude:snow_sum, ~as.numeric(scale(.))))
  
veg.comm = veg %>% 
  filter(plot_id_releve %in% veg.clim$plot_id_releve)%>%
  select(plot_id_releve, species_cover, species)%>%
  rename(cover = species_cover)%>%
  distinct()%>%
  filter(!is.na(species))%>%
  mutate(cover = ifelse(is.na(cover), 0.001, cover))%>% #there are species that are super low abundance like one leaf.
  group_by(plot_id_releve, species)%>%
  summarize(cover = sum(cover))%>%
  group_by(plot_id_releve)%>%
  mutate(total_cover = sum(cover))%>%
  group_by(plot_id_releve, species)%>%
  mutate(rel_cover = cover/total_cover)%>%
  group_by(plot_id_releve)%>%
  mutate(total_cover = sum(rel_cover))%>% # Check if all sums up to 0
  select(-total_cover, -cover)%>%
  ungroup()

veg.PA = 
  veg.comm %>%
  mutate(rel_cover = ifelse(rel_cover>0, 1, 0))%>%
  pivot_wider(names_from = species, values_from = rel_cover, values_fill = 0)%>%
  column_to_rownames(var = "plot_id_releve")

tp.species =
  tp %>%
  dplyr::select(UniqueID, Rel_Cover, species)%>%
  distinct()%>%
  filter(!is.na(species))%>%
  ungroup()%>%
  mutate(Rel_Cover = ifelse(Rel_Cover>0, 1, 0))%>%
  pivot_wider(names_from = species, values_from = Rel_Cover, values_fill = 0)%>%
  column_to_rownames(var = "UniqueID")%>%
  summarise(across(everything(), ~sum(., na.rm = TRUE)))%>%
  pivot_longer(cols = everything())%>%
  filter(value<length(unique(tp$UniqueID))*0.01)%>%
  pull(name)

veg.rare = 
  veg.PA%>%
  summarise(across(everything(), ~sum(., na.rm = TRUE)))%>%
  pivot_longer(cols = everything())%>%
  filter(value<nrow(veg.PA)*0.01)%>%
  pull(name)

# 1) Take out species if rare species are common to both datasets
takeout.sp = intersect(veg.rare, tp.species)  
# 2) Take out species if species are rare in Calanda and it does not exist in TP
takeout.sp = c(takeout.sp, veg.rare[!veg.rare %in% unique(tp$species)])

veg.PA = 
  veg.comm %>%
  filter(!species %in% takeout.sp)%>%
  mutate(rel_cover = ifelse(rel_cover>0, 1, 0))%>%
  pivot_wider(names_from = species, values_from = rel_cover, values_fill = 0)%>%
  column_to_rownames(var = "plot_id_releve")

veg.abund = 
  veg.comm %>%
  filter(!species %in% takeout.sp)%>%
  mutate(rel_cover = rel_cover)%>%
  pivot_wider(names_from = species, values_from = rel_cover, values_fill = 0)%>%
  column_to_rownames(var = "plot_id_releve")

Y = as.matrix(veg.PA)
colSums(Y)
X = veg.env[rownames(veg.PA), ]
# lambda.env = 0.001
# alpha.env = 1.0
# lambda.sp = 0.002
# alpha.sp = 0.2
# lambda.bio = 0.001
# alpha.bio = 1.0

save.image("output/starter_data.RData")
calanda_sp_list = colnames(veg.PA)
saveRDS(calanda_sp_list, file = "output/calanda_sp_list.RDS")
