# ==============================================================================
# TERNARY PLOT OF VARIANCE COMPONENTS
# ==============================================================================
cat("\n=== Creating ternary plot ===\n")

# Extract sites with variance components
sites =
  res$internals$Sites %>%
  rownames_to_column("plot_id_releve") %>%
  left_join(veg.env %>% as.data.frame() %>% rownames_to_column("plot_id_releve")) %>%
  left_join(veg %>% select(plot_id_releve, open_close) %>% distinct()) %>%
  select(plot_id_releve, env, spa, codist, Longitude, Latitude, open_close)

# Normalize by row sums
sites = sites %>%
  mutate(
    total = env + spa + codist,
    env_prop = env / total,
    spa_prop = spa / total,
    codist_prop = codist / total
  ) %>%
  filter(!is.nan(codist_prop))

# Create ternary plot with shapes for open/close
p_ternary = ggtern(sites, aes(x = env_prop, y = codist_prop, z = spa_prop)) +
  geom_point(aes(color = rgb(spa_prop, codist_prop, env_prop), shape = open_close),
             size = 3, alpha = 0.7) +
  scale_color_identity() +
  scale_shape_manual(
    values = c("open" = 16, "close" = 17),
    labels = c("Open", "Close"),
    name = "Habitat"
  ) +
  labs(
    x = "Environment",
    y = "Species associations",
    z = "Space",
    title = "Ternary plot of variance components",
    subtitle = "Colors represent RGB mixing | Shapes indicate habitat type"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    tern.axis.title = element_text(size = 14, face = "bold"),
    legend.position = "right"
  )

cat("Ternary plot created\n")

# ==============================================================================
# 2D CONTOUR PLOT WITH DEM
# ==============================================================================
cat("\n=== Creating 2D contour plot with DEM ===\n")

# Convert sites to sf object
sites_sf = sites %>%
  select(-Latitude, -Longitude) %>%
  left_join(veg %>% select(plot_id_releve, Longitude, Latitude)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

# Load Calanda mask and DEM
calanda_mask = st_read("output/calanda_mask.shp", quiet = TRUE)
dem_original = rast("output/metrics/dem.tif")

# Crop to Calanda mask extent
dem_cropped = crop(dem_original, calanda_mask)

# Aggregate to reduce resolution for faster plotting
cat("Aggregating DEM for faster plotting...\n")
dem_agg = aggregate(dem_cropped, fact = 3, fun = mean)

# Create RGB colors for visualization
sites_sf = sites_sf %>%
  mutate(rgb_color = rgb(spa_prop, codist_prop, env_prop))

# Create contour plot
p_contour = ggplot() +
  tidyterra::geom_spatraster(data = dem_agg, alpha = 0.7) +
  scale_fill_gradient(low = "white", high = "black", name = "Elevation") +
  tidyterra::geom_spatraster_contour(data = dem_agg,
                                     color = "white", alpha = 0.5, linewidth = 0.3) +
  geom_sf(data = calanda_mask, fill = NA, color = "black", linewidth = 1) +
  geom_sf(data = sites_sf, aes(color = rgb_color, shape = open_close),
          size = 2, alpha = 0.8) +
  scale_color_identity() +
  scale_shape_manual(
    values = c("open" = 16, "close" = 17),
    labels = c("Open", "Close"),
    name = "Habitat"
  ) +
  coord_sf() +
  labs(
    title = "Calanda: RGB variance components on DEM",
    subtitle = "Contour lines show elevation | Point colors show variance components",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    legend.position = "right"
  )

cat("DEM contour plot created\n")

# ==============================================================================
# VARIANCE PARTITIONING OVERVIEW
# ==============================================================================
cat("\n=== Creating variance partitioning plots ===\n")

# Prepare data for violin plot
violin_data = sites %>%
  select(env, codist, spa, open_close) %>%
  rename(Environment = env,
         `Species associations` = codist,
         Space = spa) %>%
  pivot_longer(cols = Environment:Space)

# Create summary for labels (one per violin)
label_data = violin_data %>%
  group_by(name, open_close) %>%
  summarize(
    median_val = median(value),
    max_val = max(value),
    .groups = "drop"
  ) %>%
  mutate(
    # Position labels at the top of each violin
    label_y = max_val * 1.05,
    # Create dodge position manually
    name_numeric = as.numeric(factor(name)),
    dodge_offset = ifelse(open_close == "open", -0.2, 0.2),
    x_pos = name_numeric + dodge_offset
  )

p_variance = ggplot(violin_data, aes(name, value, color = name, group = interaction(open_close, name))) +
  geom_violin(trim = FALSE, position = position_dodge(width = 0.8)) +
  geom_jitter(alpha = 0.1) +
  stat_summary(fun = median, geom = "crossbar",
               width = 0.6, fatten = 2, size = 0.7, position = position_dodge(width = 0.8))+
  stat_summary(fun = median, geom = "point",size = 3 , position = position_dodge(width = 0.8))+
  stat_summary(fun = median, geom = "text",
               aes(label = sprintf("%.2f", ..y..)),
               vjust = -1, size = 3.5 , position = position_dodge(width = 0.8)) +
  # Add single label per violin
  geom_text(data = label_data,
            aes(x = x_pos, y = 0.50, label = open_close, color = "grey50"),
            size = 5, fontface = "bold",
            inherit.aes = FALSE) +
  theme_minimal() +
  xlab("Component") +
  ylab("R² (Variance explained)")+
  scale_color_manual(values = c(
    "Environment" = rgb(red = 0,blue =1, green =0),
    "Space" = rgb(red = 1,blue =0, green =0),
    "Species associations" = rgb(red = 0,blue =0, green =1)
  ))+guides(color = "none")+
  theme(text = element_text(size = 15))

cat("Variance partitioning plot created\n")

p_histogram = sites %>%
  select(env, codist, spa) %>%
  rename(Environment = env,
         `Species associations` = codist,
         Space = spa) %>%
  pivot_longer(cols = Environment:Space) %>%
  ggplot(aes(log(value)))+
  geom_histogram()+
  facet_grid(~name)+
  theme_minimal()

cat("Variance partitioning plots created\n")
