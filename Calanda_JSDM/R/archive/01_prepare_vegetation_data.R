load("data/vegetation/TransPlantNetwork_101024.RData")

# Traits ----
if(file.exists("output/traits.csv")){
  traits = read_csv("output/traits.csv")
}else{
  try =  read_csv("data/traits/try.quantitative_traits_2025-03-17.csv")
  indicators = read_csv("data/traits/indicators_cleaned_calanda_2025-03-17.csv")
  dispersal = read_csv("data/traits/dispersal_cleaned_calanda_2025-03-17.csv")
  
  traits =
    try %>%
    filter(Trait %in% c("N_percent", "LDMC", "vegetative_height", "SLA", "LA"))%>%
    filter(species %in% colnames(Y))%>%
    filter(!is.na(Value), !is.na(Trait))%>%
    filter(!Climate_code %in% c("Af", "Am", "Aw", "BWh", "BWk", "BSh", "BSk"))%>%
    group_by(species, species_TNRS, Trait)%>%
    summarize(Value = mean(Value, na.rm = TRUE))%>%
    pivot_wider(names_from = Trait, values_from = Value)%>%
    ungroup()
  
  traits = 
    left_join(dispersal %>% select(-`...1`), traits)%>%
    left_join(indicators %>% select(-`...1`))%>%
    select(species, seed_mass, dispersal_distance_class, LA, SLA, LDMC, vegetative_height, N_percent,
           Light, Moisture, Nutrients, Disturbance.Severity
    )%>%
    rename(disturbance = Disturbance.Severity,
           dispersal = dispersal_distance_class)
  
  imp.traits = impute_functional_traits(
    traits,
    variables_to_impute = c("seed_mass", "dispersal", "LA", "SLA", "LDMC", "vegetative_height", 
                            "N_percent", "Light", "Moisture", "Nutrients", "disturbance"),
    m = 30,
    maxiter = 100,
    num.trees = 500,
    seed = 123,
    validation_fraction = 0.3
  )
  traits %>%
    summarise(across(everything(), ~mean(is.na(.)) * 100)) %>%
    pivot_longer(cols = everything(), 
                 names_to = "variable", 
                 values_to = "percent_missing") %>%
    arrange(desc(percent_missing)) %>%
    right_join(imp.traits$performance) %>%
    mutate(percent_missing = round(percent_missing),
           r_squared = round(r_squared, 2)) %>%
    select(variable, percent_missing, r_squared) %>%
    gt()
  
  traits = 
    imp.traits$imputed_data %>%
    select(species, contains("final"))%>%
    rename_with(~str_remove(., "_final$"), ends_with("_final"))
  write.csv(traits, file = "output/traits.csv")
}

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

veg = 
  veg %>%
  rename(species = taxon_global)%>%
  filter(!is.na(word(species, 2))) #Take out genus

veg.comm = veg %>% 
  #filter(plot_id_releve %in% veg.clim$plot_id_releve)%>%
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
  ungroup()%>%
  mutate(species = ifelse(species == "Silene vulgaris subsp. glareosa", "Silene vulgaris", species))

veg.abund = 
  veg.comm %>%
  mutate(rel_cover = rel_cover)%>%
  pivot_wider(names_from = species, values_from = rel_cover, values_fill = 0)%>%
  column_to_rownames(var = "plot_id_releve")

community_traits =
  calculate_community_traits_extended(
  community_data = veg.abund,  
  traits_data = traits[,-1],           
  #species_col = "species",         
  abundance_col = "rel_cover"                
)

# Create a PCA for land-use ----
res.pca = left_join(community_traits %>% select(plot_id_releve, Nutrients, disturbance),
                    veg.abund %>% rownames_to_column(var = "plot_id_releve")%>% select(plot_id_releve, `Nardus stricta`))%>%
  #rename(nardus_stricta_abundance = `Nardus stricta`)%>%
  na.omit()%>%
  column_to_rownames(var = "plot_id_releve")

res.pca = PCA(res.pca, scale.unit = TRUE)

