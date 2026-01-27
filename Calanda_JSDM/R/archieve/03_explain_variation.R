# Explain the variation in R2 with climate ----
an = readRDS("results_from_Max/an_sjsdm_calanda.rds")
model = readRDS("results_from_Max/model_sjsdm_calanda.rds")
res = internalStructure(an, fractions = "proportional")

## Community point of view ----
sites = res$internals$Sites
# Create RGB colors for plant communities
# Process sites data with RGB colors
sites = 
  sites %>%
  rownames_to_column("plot_id_releve")%>%
  left_join(veg.env %>% as.data.frame() %>% rownames_to_column("plot_id_releve"))%>%
  mutate(
    # Store RGB values for other plots
    R = spa,
    G = codist,
    B = env,
    color = rgb(R, G, B)
  )

# Link to the environment gradient 
res.pca = PCA(sites %>% select(slope, snow_sum, fdd, rocks_cover, trees_cover, summer_temp, land_use, soil_depth_var), scale.unit = TRUE)

# Extract PCA coordinates
ind_coords = as.data.frame(get_pca_ind(res.pca)$coord)
var_coords = as.data.frame(get_pca_var(res.pca)$coord)

# Calculate scaling factor for arrows (make them longer)
scaling_factor = 5  # Adjust this value to make arrows longer/shorter

# Extract PCA coordinates
ind_coords = as.data.frame(get_pca_ind(res.pca)$coord)
var_coords = as.data.frame(get_pca_var(res.pca)$coord)

# Add colors to ind_coords
ind_coords$colors = gsub('"', '', sites$color)  # Remove quotes from hex codes

# Calculate scaling factor for arrows
scaling_factor = 5  

p1.1 = 
ggplot() +
  # Add arrows for variables
  geom_segment(data = var_coords,
               aes(x = 0, y = 0, 
                   xend = Dim.1 * scaling_factor, 
                   yend = Dim.2 * scaling_factor),
               arrow = arrow(length = unit(0.5, "cm")),
               color = "black",
               size = 0.5) +
  
  # Add variable labels
  geom_text_repel(data = var_coords,
                  aes(x = Dim.1 * scaling_factor,
                      y = Dim.2 * scaling_factor,
                      label = rownames(var_coords)),
                  color = "black",
                  max.overlaps = Inf) +
  
  # Add points with custom colors
  geom_point(data = ind_coords,
             aes(x = Dim.1, y = Dim.2, color = I(colors)),  # Note the I() function
             size = 3,
             alpha = 0.6) +
  
  # Add theme elements
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # Add labels
  labs(title = "PCA - Biplot",
       x = paste0("Dim1 (", round(res.pca$eig[1,2], 1), "%)"),
       y = paste0("Dim2 (", round(res.pca$eig[2,2], 1), "%)")) +
  
  # Make it more square-shaped
  coord_fixed(ratio = 1) +
  
  # Expand the plot margins
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  scale_y_continuous(expand = expansion(mult = 0.2))


p1.2 = 
  ggplot() +
  # Add arrows for variables
  geom_segment(data = var_coords,
               aes(x = 0, y = 0, 
                   xend = Dim.2 * scaling_factor, 
                   yend = Dim.3 * scaling_factor),
               arrow = arrow(length = unit(0.5, "cm")),
               color = "black",
               size = 0.5) +
  
  # Add variable labels
  geom_text_repel(data = var_coords,
                  aes(x = Dim.2 * scaling_factor,
                      y = Dim.3 * scaling_factor,
                      label = rownames(var_coords)),
                  color = "black",
                  max.overlaps = Inf) +
  
  # Add points with custom colors
  geom_point(data = ind_coords,
             aes(x = Dim.2, y = Dim.3, color = I(colors)),  # Note the I() function
             size = 3,
             alpha = 0.6) +
  
  # Add theme elements
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # Add labels
  labs(title = "PCA - Biplot",
       x = paste0("Dim2 (", round(res.pca$eig[2,2], 1), "%)"),
       y = paste0("Dim3 (", round(res.pca$eig[3,2], 1), "%)")) +
  
  # Make it more square-shaped
  coord_fixed(ratio = 1) +
  
  # Expand the plot margins
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  scale_y_continuous(expand = expansion(mult = 0.2))

# Call the function with your sites data
p2 = create_ternary(sites)

# Create plot with only significant relationships -scaled
p3 = sites %>% 
  select(env, spa, codist, r2, fdd, land_use, et.annual, summer_temp, soil_depth_var, slope) %>%
  pivot_longer(cols = env:r2, names_to = "process", values_to = "variance_values") %>%
  group_by(process) %>%
  nest() %>%
  mutate(
    model = map(data, ~lm(variance_values ~ fdd + land_use + et.annual + summer_temp + soil_depth_var + slope, 
                          data = .x)),
    coef_table = map(model, ~as.data.frame(summary(.x)$coefficients) %>%
                       rownames_to_column("var"))
  ) %>%
  unnest(coef_table) %>%
  filter(var != "(Intercept)") %>%  # Remove intercept
  mutate(
    # Create lower and upper confidence intervals
    lower_ci = Estimate - 1.96 * `Std. Error`,
    upper_ci = Estimate + 1.96 * `Std. Error`,
    # Create a significance flag for transparency
    significant = `Pr(>|t|)` < 0.05) %>%
  mutate(process = factor(process, levels = c("env", "spa", "codist", "r2")))%>%
  ggplot(aes(x = var, y = Estimate, color = process)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 3, aes(alpha = significant)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci, alpha = significant), 
                width = 0.2, position = position_dodge(width = 0.5)) +
  facet_wrap(~process, ncol = 4) +
  coord_flip() +  # Horizontal layout for better readability
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3)) +
  scale_color_manual(values = c(
    "spa" = rgb(1, 0.3, 0.3),
    "env" = rgb(0.3, 0.3, 1),
    "codist" = rgb(0.3, 1, 0.3),
    "r2" = "grey50"
  )) +
  labs(
    y = "Coefficient Estimate",
    x = "Predictor Variable",
    title = "Effects of Environmental Variables on Different Processes",
    subtitle = "Transparent bars indicate non-significant effects (p > 0.05)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank()
  )

p3
  
# Create plot with only significant relationships -nonscaled
p4 = sites %>% 
  select(env, spa, codist, r2, slope:rocks_cover)%>%
  pivot_longer(cols = env:r2, names_to = "process", values_to = "variance_values") %>%
  pivot_longer(cols = slope:rocks_cover, names_to = "var", values_to = "var_values") %>%
  group_by(process, var) %>%
  nest() %>%
  mutate(
    model = map(data, ~lm(variance_values ~ var_values, data = .x)),
    p_value = map_dbl(model, ~summary(.x)$coefficients[2,4])) %>%
  filter(p_value < 0.05) %>% 
  unnest(data) %>%
  ggplot(aes(var_values, variance_values, color = process)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(var~., scales = "free") +
  theme_bw() +
  labs(
    x = "Variable Values",
    y = "Variance Component",
    color = "Variables"
  ) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "white")
  )+
  scale_color_manual(values = c("spa" = rgb(1,0.3,0.3),
                                "env" = rgb(0.3,0.3, 1),
                                "codist" = rgb(0.3,1,0.3),
                                "r2" = "grey50"))
p4

pdf("plot/calanda_jsdm_outputs.pdf", height = 6, width = 9)
plot(an)
plot(res)
plotAssemblyEffects(res)
dev.off()

pdf("plot/calanda_jsdm_var_communities.pdf", height = 10, width = 15)
print(p1.1+p2)
print(p1.2+p2)
print(p3)
print(p4)
dev.off()

## Species point of view ----
species = res$internals$Species

# Calculate species environmental means 
veg.env.raw = 
  veg.env = 
  veg.clim %>% 
  group_by(plot_id_releve, Longitude, Latitude, altitude)%>%
  summarize(across(slope:land_use, ~mean(.)))%>%
  ungroup()%>%
  column_to_rownames(var = "plot_id_releve")%>%
  na.omit(.)
  
sp_env = calculate_species_env_means(species_data = Y, env_data = veg.env.raw)

# Create RGB colors for plant communities
# Process sites data with RGB colors
species = species %>%
  rownames_to_column("species") %>%
  left_join(sp_env) %>%
  mutate(
    spa_scaled = (spa - min(spa, na.rm = TRUE)) / (max(spa, na.rm = TRUE) - min(spa, na.rm = TRUE)),
    codist_scaled = (codist - min(codist, na.rm = TRUE)) / (max(codist, na.rm = TRUE) - min(codist, na.rm = TRUE)),
    env_scaled = (env - min(env, na.rm = TRUE)) / (max(env, na.rm = TRUE) - min(env, na.rm = TRUE))
  ) %>%
  mutate(
    R = spa_scaled,
    G = codist_scaled, 
    B = env_scaled,
    color = rgb(R, G, B)
  )


# Link to the environment gradient
res.pca = PCA(species %>% select(slope:land_use), scale.unit = TRUE)

# Extract PCA coordinates
ind_coords = as.data.frame(get_pca_ind(res.pca)$coord)
var_coords = as.data.frame(get_pca_var(res.pca)$coord)

# Calculate scaling factor for arrows (make them longer)
scaling_factor = 5  # Adjust this value to make arrows longer/shorter

# Extract PCA coordinates
ind_coords = as.data.frame(get_pca_ind(res.pca)$coord)
var_coords = as.data.frame(get_pca_var(res.pca)$coord)

# Add colors to ind_coords
ind_coords$colors = gsub('"', '', species$color)  # Remove quotes from hex codes

# Calculate scaling factor for arrows
scaling_factor = 5  

p1.1  = 
  ggplot() +
  # Add arrows for variables
  geom_segment(data = var_coords,
               aes(x = 0, y = 0, 
                   xend = Dim.1 * scaling_factor, 
                   yend = Dim.2 * scaling_factor),
               arrow = arrow(length = unit(0.5, "cm")),
               color = "black",
               size = 0.5) +
  
  # Add variable labels
  geom_text_repel(data = var_coords,
                  aes(x = Dim.1 * scaling_factor,
                      y = Dim.2 * scaling_factor,
                      label = rownames(var_coords)),
                  color = "black",
                  max.overlaps = Inf) +
  
  # Add points with custom colors
  geom_point(data = ind_coords,
             aes(x = Dim.1, y = Dim.2, color = I(colors)),  # Note the I() function
             size = 3,
             alpha = 0.6) +
  
  # Add theme elements
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # Add labels
  labs(title = "PCA - Biplot",
       x = paste0("Dim1 (", round(res.pca$eig[1,2], 1), "%)"),
       y = paste0("Dim2 (", round(res.pca$eig[2,2], 1), "%)")) +
  
  # Make it more square-shaped
  coord_fixed(ratio = 1) +
  
  # Expand the plot margins
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  scale_y_continuous(expand = expansion(mult = 0.2))

p1.1

p1.2 = 
  ggplot() +
  # Add arrows for variables
  geom_segment(data = var_coords,
               aes(x = 0, y = 0, 
                   xend = Dim.2 * scaling_factor, 
                   yend = Dim.3 * scaling_factor),
               arrow = arrow(length = unit(0.5, "cm")),
               color = "black",
               size = 0.5) +
  
  # Add variable labels
  geom_text_repel(data = var_coords,
                  aes(x = Dim.2 * scaling_factor,
                      y = Dim.3 * scaling_factor,
                      label = rownames(var_coords)),
                  color = "black",
                  max.overlaps = Inf) +
  
  # Add points with custom colors
  geom_point(data = ind_coords,
             aes(x = Dim.2, y = Dim.3, color = I(colors)),  # Note the I() function
             size = 3,
             alpha = 0.6) +
  
  # Add theme elements
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # Add labels
  labs(title = "PCA - Biplot",
       x = paste0("Dim2 (", round(res.pca$eig[2,2], 1), "%)"),
       y = paste0("Dim3 (", round(res.pca$eig[3,2], 1), "%)")) +
  
  # Make it more square-shaped
  coord_fixed(ratio = 1) +
  
  # Expand the plot margins
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  scale_y_continuous(expand = expansion(mult = 0.2))
p1.2

# Call the function with your sites data
p2 = create_ternary(species)

# Create plot with only significant relationships - ratios
p3 = species %>% 
  select(env, spa, codist, r2, fdd, land_use, et.annual, summer_temp, soil_depth_var, slope) %>%
  mutate(across(fdd:slope, ~as.numeric(scale(.))))%>%
  pivot_longer(cols = env:r2, names_to = "process", values_to = "variance_values") %>%
  group_by(process) %>%
  nest() %>%
  mutate(
    model = map(data, ~lm(variance_values ~ fdd + land_use + et.annual + summer_temp + soil_depth_var + slope, 
                          data = .x)),
    coef_table = map(model, ~as.data.frame(summary(.x)$coefficients) %>%
                       rownames_to_column("var"))
  ) %>%
  unnest(coef_table) %>%
  filter(var != "(Intercept)") %>%  # Remove intercept
  mutate(
    # Create lower and upper confidence intervals
    lower_ci = Estimate - 1.96 * `Std. Error`,
    upper_ci = Estimate + 1.96 * `Std. Error`,
    # Create a significance flag for transparency
    significant = `Pr(>|t|)` < 0.05
  ) %>%
  mutate(process = factor(process, levels = c("env", "spa", "codist", "r2")))%>%
  ggplot(aes(x = var, y = Estimate, color = process)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 3, aes(alpha = significant)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci, alpha = significant), 
                width = 0.2, position = position_dodge(width = 0.5)) +
  facet_wrap(~process, ncol = 4) +
  coord_flip() +  # Horizontal layout for better readability
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3)) +
  scale_color_manual(values = c(
    "spa" = rgb(1, 0.3, 0.3),
    "env" = rgb(0.3, 0.3, 1),
    "codist" = rgb(0.3, 1, 0.3),
    "r2" = "grey50"
  )) +
  labs(
    y = "Coefficient Estimate",
    x = "Predictor Variable",
    title = "Effects of Environmental Variables on Different Processes",
    subtitle = "Transparent bars indicate non-significant effects (p > 0.05)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank()
  )

p3

# Create plot with only significant relationships - not ratio
p4 = species %>% 
  select(env, spa, codist, r2, slope:land_use)%>%
  pivot_longer(cols = env:r2, names_to = "process", values_to = "variance_values") %>%
  pivot_longer(cols = slope:land_use, names_to = "var", values_to = "var_values") %>%
  group_by(process, var) %>%
  nest() %>%
  mutate(
    model = map(data, ~lm(variance_values ~ var_values, data = .x)),
    p_value = map_dbl(model, ~summary(.x)$coefficients[2,4])) %>%
  filter(p_value < 0.05) %>% 
  unnest(data) %>%
  ggplot(aes(var_values, variance_values, color = process)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(var~., scales = "free") +
  theme_bw() +
  ylim(0, 1)+
  labs(
    x = "Variable Values",
    y = "Variance Component",
    color = "Variables"
  ) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "white")
  )+
  scale_color_manual(values = c("spa" = rgb(1,0.3,0.3),
                                "env" = rgb(0.3,0.3, 1),
                                "codist" = rgb(0.3,1,0.3),
                                "r2" = "black"))

p4

pdf("plot/calanda_jsdm_var_species.pdf", height = 10, width = 15)
print(p1.1+p2)
print(p1.2+p2)
print(p3)
print(p4)
dev.off()

# Traits ----
## Clean traits and impute ----
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

# Check the percentage of missing values for each trait and imputation performance
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

## Community point of view ----
sites = res$internals$Sites
# Get CM traits

community_traits <- calculate_community_traits(
  community_data = Y,  
  traits_data = traits,           
  species_col = "species",         
  abundance_col = NULL                
)

# Create RGB colors for plant communities
# Process sites data with RGB colors
sites = 
  sites %>%
  rownames_to_column("plot_id_releve")%>%
  left_join(community_traits)%>%
  mutate(
    # Scale the variables to ensure they sum to 1 for ternary plotting
    total = spa + codist + env,
    spa_scaled = spa/total,
    codist_scaled = codist/total,
    env_scaled = env/total,
    
    # Store RGB values for other plots
    R = spa_scaled,
    G = codist_scaled, 
    B = env_scaled,
    color = rgb(R, G, B)
  )

p1.1 = create_pca_biplot(
  data = sites,
  variables = c("LDMC_cwm", "SLA_cwm", "N_percent_cwm", "vegetative_height_cwm", "seed_mass_cwm"),
  point_size = 2, point_alpha = 0.4
)
p1.1

p1.2 = create_pca_biplot(
  data = sites,
  variables = c("Moisture_cwm", "Nutrients_cwm", "Light_cwm", "disturbance_cwm"),
  point_size = 2, point_alpha = 0.6
)
p1.2

p1.3 = create_pca_biplot(
  data = sites,
  variables = c("seed_mass_cwm", "dispersal_cwm", "vegetative_height_cwm"),
  point_size = 2, point_alpha = 0.6
)
p1.3

# Call the function with your sites data
p2 = create_ternary(sites)

# Create plot with only significant relationships -scaled
p3 = sites %>% 
  select(env_scaled, spa_scaled, codist_scaled, contains("cwm"))%>%
  pivot_longer(cols = env_scaled:codist_scaled, names_to = "process", values_to = "variance_values") %>%
  pivot_longer(cols = contains("cwm"), names_to = "var", values_to = "var_values") %>%
  group_by(process, var) %>%
  nest() %>%
  mutate(
    model = map(data, ~lm(variance_values ~ var_values, data = .x)),
    p_value = map_dbl(model, ~summary(.x)$coefficients[2,4])) %>%
  filter(p_value < 0.05) %>% 
  unnest(data) %>%
  ggplot(aes(var_values, variance_values, color = process)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(var~., scales = "free") +
  theme_bw() +
  labs(
    x = "Variable Values",
    y = "Variance Component",
    color = "Variables"
  ) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "white")
  )+
  scale_color_manual(values = c("spa_scaled" = rgb(1,0.3,0.3),
                                "env_scaled" = rgb(0.3,0.3, 1),
                                "codist_scaled" = rgb(0.3,1,0.3)))
p3

# Create plot with only significant relationships -nonscaled
p4 = sites %>% 
  select(env, spa, codist, contains("cwm"))%>%
  pivot_longer(cols = env:codist, names_to = "process", values_to = "variance_values") %>%
  pivot_longer(cols = contains("cwm"), names_to = "var", values_to = "var_values") %>%
  group_by(process, var) %>%
  nest() %>%
  mutate(
    model = map(data, ~lm(variance_values ~ var_values, data = .x)),
    p_value = map_dbl(model, ~summary(.x)$coefficients[2,4])) %>%
  filter(p_value < 0.05) %>% 
  unnest(data) %>%
  ggplot(aes(var_values, variance_values, color = process)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(var~., scales = "free") +
  theme_bw() +
  labs(
    x = "Variable Values",
    y = "Variance Component",
    color = "Variables"
  ) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "white")
  )+
  scale_color_manual(values = c("spa" = rgb(1,0.3,0.3),
                                "env" = rgb(0.3,0.3, 1),
                                "codist" = rgb(0.3,1,0.3)))
p4


pdf("plot/calanda_jsdm_var_communities_traits.pdf", height = 10, width = 15)
print(p1.1$plot+p2)
print(p1.2$plot+p2)
print(p3)
print(p4)
dev.off()





## Species point of view ----
species = res$internals$Species

# Create RGB colors for plant communities
# Process sites data with RGB colors
species = 
  species %>%
  rownames_to_column("species")%>%
  left_join(traits)%>%
  mutate(
    # Scale the variables to ensure they sum to 1 for ternary plotting
    total = spa + codist + env,
    spa_scaled = spa/total,
    codist_scaled = codist/total,
    env_scaled = env/total) %>%
  filter(total != 0)%>%
  mutate(
    # Store RGB values for other plots
    R = spa_scaled,
    G = codist_scaled, 
    B = env_scaled,
    color = rgb(R, G, B)
  )

# Link to the environment gradient
res.pca = PCA(species %>% select(SLA, vegetative_height, N_percent, LA, LDMC), scale.unit = TRUE)

# Extract PCA coordinates
ind_coords = as.data.frame(get_pca_ind(res.pca)$coord)
var_coords = as.data.frame(get_pca_var(res.pca)$coord)

# Calculate scaling factor for arrows (make them longer)
scaling_factor = 5  # Adjust this value to make arrows longer/shorter

# Extract PCA coordinates
ind_coords = as.data.frame(get_pca_ind(res.pca)$coord)
var_coords = as.data.frame(get_pca_var(res.pca)$coord)

# Add colors to ind_coords
ind_coords$colors = gsub('"', '', species$color)  # Remove quotes from hex codes

# Calculate scaling factor for arrows
scaling_factor = 5  

p1.1  = 
  ggplot() +
  # Add arrows for variables
  geom_segment(data = var_coords,
               aes(x = 0, y = 0, 
                   xend = Dim.1 * scaling_factor, 
                   yend = Dim.2 * scaling_factor),
               arrow = arrow(length = unit(0.5, "cm")),
               color = "black",
               size = 0.5) +
  
  # Add variable labels
  geom_text_repel(data = var_coords,
                  aes(x = Dim.1 * scaling_factor,
                      y = Dim.2 * scaling_factor,
                      label = rownames(var_coords)),
                  color = "black",
                  max.overlaps = Inf) +
  
  # Add points with custom colors
  geom_point(data = ind_coords,
             aes(x = Dim.1, y = Dim.2, color = I(colors)),  # Note the I() function
             size = 3,
             alpha = 0.6) +
  
  # Add theme elements
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # Add labels
  labs(title = "PCA - Biplot",
       x = paste0("Dim1 (", round(res.pca$eig[1,2], 1), "%)"),
       y = paste0("Dim2 (", round(res.pca$eig[2,2], 1), "%)")) +
  
  # Make it more square-shaped
  coord_fixed(ratio = 1) +
  
  # Expand the plot margins
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  scale_y_continuous(expand = expansion(mult = 0.2))

p1.1

# Call the function with your sites data
p2 = create_ternary(species)

# Create plot with only significant relationships - ratios
p3 = species %>% 
  select(env, spa, codist, seed_mass:disturbance)%>%
  pivot_longer(cols = env:codist, names_to = "process", values_to = "variance_values") %>%
  pivot_longer(cols = seed_mass:disturbance, names_to = "var", values_to = "var_values") %>%
  group_by(process, var) %>%
  nest() %>%
  mutate(
    model = map(data, ~lm(variance_values ~ var_values, data = .x)),
    p_value = map_dbl(model, ~summary(.x)$coefficients[2,4])) %>%
  filter(p_value < 0.05) %>% 
  unnest(data) %>%
  ggplot(aes(var_values, variance_values, color = process)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(var~., scales = "free") +
  theme_bw() +
  ylim(0, 1)+
  labs(
    x = "Variable Values",
    y = "Variance Component",
    color = "Variables"
  ) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "white")
  )+
  scale_color_manual(values = c("spa" = rgb(1,0.3,0.3),
                                "env" = rgb(0.3,0.3, 1),
                                "codist" = rgb(0.3,1,0.3)))

p3

# Create plot with only significant relationships - not ratio
p4 = species %>% 
  select(env_scaled, spa_scaled, codist_scaled,seed_mass:disturbance)%>%
  pivot_longer(cols = env_scaled:codist_scaled, names_to = "process", values_to = "variance_values") %>%
  pivot_longer(cols = seed_mass:disturbance, names_to = "var", values_to = "var_values") %>%
  group_by(process, var) %>%
  nest() %>%
  mutate(
    model = map(data, ~lm(variance_values ~ var_values, data = .x)),
    p_value = map_dbl(model, ~summary(.x)$coefficients[2,4])) %>%
  filter(p_value < 0.05) %>% 
  unnest(data) %>%
  ggplot(aes(var_values, variance_values, color = process)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(var~., scales = "free") +
  theme_bw() +
  ylim(0, 1)+
  labs(
    x = "Variable Values",
    y = "Variance Component",
    color = "Variables"
  ) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "white")
  )+
  scale_color_manual(values = c("spa_scaled" = rgb(1,0.3,0.3),
                                "env_scaled" = rgb(0.3,0.3, 1),
                                "codist_scaled" = rgb(0.3,1,0.3)))

p4

pdf("plot/calanda_jsdm_var_species_traits.pdf", height = 10, width = 15)
print(p1.1+p2)
print(p3)
print(p4)
dev.off()
