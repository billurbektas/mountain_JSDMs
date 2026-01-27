library(readxl)
load("~/Desktop/Common_JSDM/Calanda_JSDM/output/starter_data.RData")
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

focal_species =
traits %>%
  filter(is.na(N_percent))%>%
  pull(species)

focal_species %in% colnames(Y) # all in species x sites matrix.

caphe =
bind_rows(
read_csv("~/Downloads/2024_CalTraits_Data_1-500.csv"),
read_csv("~/Downloads/2024_CalTraits_Data_1501-2000+.csv"),
read_csv("~/Downloads/2024_CalTraits_Data_1001-1500.csv"),
read_csv("~/Downloads/2024_CalTraits_Data_501-1000.csv"))%>%
rename(state_phenology = `state phenology`)%>%
  select(plant_ID, plant_species, individual_number, site, state_phenology, biomass)%>%
  distinct()

# Define phenology state categories for grouping
categorize_phenology = function(pheno) {
  if(is.na(pheno)) return("other")
  
  has_B = str_detect(pheno, "B")
  has_F = str_detect(pheno, "F") 
  has_D = str_detect(pheno, "D")
  has_V = str_detect(pheno, "V")
  
  # Categorize phenology states according to priority rules
  if(has_B || has_F) return("BF")
  if(has_V) return("V") 
  if(has_D) return("D")
  return("other")
}

# First calculate the number of sites per species AND phenology category
caphe_with_pheno_category = caphe %>%
  mutate(pheno_category = sapply(state_phenology, categorize_phenology))

sites_per_species_pheno = caphe_with_pheno_category %>%
  group_by(plant_species, pheno_category) %>%
  summarize(num_sites_pheno = n_distinct(site), .groups = "drop")

# Also calculate global site count per species (for reference)
sites_per_species = caphe %>%
  group_by(plant_species) %>%
  summarize(total_sites = n_distinct(site))

# Join this information back to the original data
caphe_with_sites = caphe_with_pheno_category %>%
  left_join(sites_per_species_pheno, by = c("plant_species", "pheno_category")) %>%
  left_join(sites_per_species, by = "plant_species")


# Create a prioritization function based on your criteria
prioritize_samples = function(data) {
  data = data %>%
    mutate(
      is_focal = plant_species %in% focal_species,
      has_BFD = if_else(!is.na(state_phenology), 
                        str_detect(state_phenology, "B|F|D"), 
                        FALSE),
      has_BF = if_else(!is.na(state_phenology), 
                       str_detect(state_phenology, "B|F"), 
                       FALSE),
      has_V = if_else(!is.na(state_phenology), 
                      str_detect(state_phenology, "V"), 
                      FALSE),
      has_D = if_else(!is.na(state_phenology), 
                      str_detect(state_phenology, "D"), 
                      FALSE),
      first_3 = individual_number <= 3
    )
  
  # Apply the priority rules
  data = data %>%
    mutate(priority = case_when(
      # Priority 1
      is_focal & has_BFD & first_3 ~ 1,
      
      # Priority 2-6: B or F phenology with different site counts
      !is_focal & has_BF & num_sites_pheno == 5 & first_3 ~ 2,
      !is_focal & has_BF & num_sites_pheno == 4 & first_3 ~ 3,
      !is_focal & has_BF & num_sites_pheno == 3 & first_3 ~ 4,
      !is_focal & has_BF & num_sites_pheno == 2 & first_3 ~ 5,
      !is_focal & has_BF & num_sites_pheno == 1 & first_3 ~ 6,
      
      # Priority 7-11: V phenology (but not B or F)
      !is_focal & has_V & !has_BF & num_sites_pheno == 5 & first_3 ~ 7,
      !is_focal & has_V & !has_BF & num_sites_pheno == 4 & first_3 ~ 8,
      !is_focal & has_V & !has_BF & num_sites_pheno == 3 & first_3 ~ 9,
      !is_focal & has_V & !has_BF & num_sites_pheno == 2 & first_3 ~ 10,
      !is_focal & has_V & !has_BF & num_sites_pheno == 1 & first_3 ~ 11,
      
      # Priority 12-16: D phenology (but not B, F, or V)
      !is_focal & has_D & !has_BF & !has_V & num_sites_pheno == 5 & first_3 ~ 12,
      !is_focal & has_D & !has_BF & !has_V & num_sites_pheno == 4 & first_3 ~ 13,
      !is_focal & has_D & !has_BF & !has_V & num_sites_pheno == 3 & first_3 ~ 14,
      !is_focal & has_D & !has_BF & !has_V & num_sites_pheno == 2 & first_3 ~ 15,
      !is_focal & has_D & !has_BF & !has_V & num_sites_pheno == 1 & first_3 ~ 16,
      
      # Priority 17-21: All remaining phenology
      !is_focal & !has_BF & !has_V & !has_D & num_sites_pheno == 5 & first_3 ~ 17,
      !is_focal & !has_BF & !has_V & !has_D & num_sites_pheno == 4 & first_3 ~ 18,
      !is_focal & !has_BF & !has_V & !has_D & num_sites_pheno == 3 & first_3 ~ 19,
      !is_focal & !has_BF & !has_V & !has_D & num_sites_pheno == 2 & first_3 ~ 20,
      !is_focal & !has_BF & !has_V & !has_D & num_sites_pheno == 1 & first_3 ~ 21,
      
      # Default: lowest priority
      TRUE ~ 22
    ))
  
  return(data)
}

# Apply the prioritization
caphe_prioritized = prioritize_samples(caphe_with_sites)

# Display result with prioritization
caphe_prioritized = caphe_prioritized %>%
  arrange(priority) %>%
  select(plant_ID, plant_species, individual_number, site, state_phenology, priority, everything())

# To save the prioritized data to a CSV file
write_csv(caphe_prioritized, "prioritized_plant_samples.csv")

# View the top records ordered by priority
head(caphe_prioritized)
# Calculate coverage statistics by priority
calculate_coverage_statistics = function(data) {
  # Get total counts for denominator
  total_species = data %>% 
    select(plant_species) %>% 
    n_distinct()
  
  total_sites = data %>%
    select(site) %>%
    n_distinct()
  
  # Calculate total possible site-species combinations
  total_site_species = data %>%
    select(site, plant_species) %>%
    distinct() %>%
    nrow()
  
  # Calculate cumulative coverage by priority
  coverage_stats = data %>%
    arrange(priority) %>%
    mutate(priority_group = if_else(priority <= 21, paste0("Priority ", priority), "Remaining")) %>%
    group_by(priority_group) %>%
    summarise(
      sample_count = n(),
      plant_ID_count = n_distinct(plant_ID),
      species_count = n_distinct(plant_species),
      site_count = n_distinct(site),
      site_species_count = n_distinct(paste0(site, "_", plant_species))
    ) %>%
    # Calculate percentages
    mutate(
      sample_percent = round(sample_count / nrow(data) * 100, 1),
      plant_ID_percent = round(plant_ID_count / n_distinct(data$plant_ID) * 100, 1),
      species_percent = round(species_count / total_species * 100, 1),
      site_percent = round(site_count / total_sites * 100, 1),
      site_species_percent = round(site_species_count / total_site_species * 100, 1)
    )
  
  # Calculate cumulative statistics
  # Create a summary of site coverage at different phenology states
  site_pheno_summary = data %>%
    group_by(plant_species, pheno_category, site) %>%
    summarize(has_samples = n() > 0, .groups = "drop") %>%
    summarize(
      species_pheno_site_combinations = n(),
      .by = c(pheno_category)
    )
  
  # Calculate cumulative statistics
  cumulative_stats = data %>%
    arrange(priority) %>%
    mutate(priority_threshold = priority) %>%
    group_by(priority_threshold) %>%
    mutate(rank = row_number()) %>%
    filter(rank == 1) %>%
    ungroup() %>%
    arrange(priority_threshold) %>%
    mutate(
      cumulative_samples = map_dbl(priority_threshold, function(p) {
        sum(data$priority <= p)
      }),
      cumulative_plant_IDs = map_dbl(priority_threshold, function(p) {
        data %>% filter(priority <= p) %>% pull(plant_ID) %>% n_distinct()
      }),
      cumulative_species = map_dbl(priority_threshold, function(p) {
        data %>% filter(priority <= p) %>% pull(plant_species) %>% n_distinct()
      }),
      cumulative_sites = map_dbl(priority_threshold, function(p) {
        data %>% filter(priority <= p) %>% pull(site) %>% n_distinct()
      }),
      cumulative_site_species = map_dbl(priority_threshold, function(p) {
        data %>% 
          filter(priority <= p) %>% 
          select(site, plant_species) %>%
          distinct() %>%
          nrow()
      }),
      cumulative_species_pheno_sites = map_dbl(priority_threshold, function(p) {
        data %>% 
          filter(priority <= p) %>% 
          group_by(plant_species, pheno_category, site) %>%
          summarize(has_samples = n() > 0, .groups = "drop") %>%
          nrow()
      })
    ) %>%
    # Calculate percentages
    mutate(
      cumulative_sample_percent = round(cumulative_samples / nrow(data) * 100, 1),
      cumulative_plant_IDs_percent = round(cumulative_plant_IDs / n_distinct(data$plant_ID) * 100, 1),
      cumulative_species_percent = round(cumulative_species / total_species * 100, 1),
      cumulative_site_percent = round(cumulative_sites / total_sites * 100, 1),
      cumulative_site_species_percent = round(cumulative_site_species / total_site_species * 100, 1),
      
      # Calculate the total possible species-phenology-site combinations
      total_species_pheno_sites = sum(site_pheno_summary$species_pheno_site_combinations),
      cumulative_species_pheno_site_percent = round(cumulative_species_pheno_sites / 
                                                      total_species_pheno_sites * 100, 1)
    ) %>%
    select(
      priority_threshold,
      cumulative_samples, cumulative_sample_percent,
      cumulative_plant_IDs, cumulative_plant_IDs_percent,
      cumulative_species, cumulative_species_percent,
      cumulative_sites, cumulative_site_percent,
      cumulative_site_species, cumulative_site_species_percent,
      cumulative_species_pheno_sites, cumulative_species_pheno_site_percent
    )
  
  return(list(
    by_priority = coverage_stats,
    cumulative = cumulative_stats
  ))
}

# Calculate coverage statistics
coverage_statistics = calculate_coverage_statistics(caphe_prioritized)

# Display the statistics
print("Coverage statistics by priority group:")
print(coverage_statistics$by_priority)

print("\nCumulative coverage by priority threshold:")
print(coverage_statistics$cumulative)

# Create visualization of the cumulative coverage
coverage_plot = ggplot(coverage_statistics$cumulative, aes(x = priority_threshold)) +
  geom_line(aes(y = cumulative_species_percent, color = "Species"), size = 1) +
  geom_line(aes(y = cumulative_sample_percent, color = "Samples"), size = 1) +
  geom_line(aes(y = cumulative_site_species_percent, color = "Site-Species Combinations"), size = 1) +
  geom_line(aes(y = cumulative_species_pheno_site_percent, color = "Species-Phenology-Site Combinations"), 
            size = 1, linetype = "dashed") +
  # Add labels for cumulative plant_IDs at each priority point
  geom_text(aes(y = cumulative_sample_percent, 
                label = paste0(cumulative_plant_IDs)), 
            vjust = -0.5, size = 3, color = "black") +
  scale_color_manual(values = c(
    "Species" = "blue", 
    "Samples" = "green",
    "Site-Species Combinations" = "orange",
    "Species-Phenology-Site Combinations" = "purple"
  )) +
  scale_x_continuous(breaks = seq(1, 22, by = 1)) +
  labs(
    title = "Cumulative Coverage by Priority",
    x = "Priority Threshold",
    y = "Percent Coverage (%)",
    color = "Category"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# To save the coverage statistics to CSV files
# write_csv(coverage_statistics$by_priority, "coverage_by_priority.csv")
# write_csv(coverage_statistics$cumulative, "cumulative_coverage.csv")

# Display the plot
print(coverage_plot)
ggplot2::ggsave("cumulative_coverage_plot.png", coverage_plot, width = 10, height = 6, dpi = 300)
