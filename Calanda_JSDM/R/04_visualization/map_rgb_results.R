# ==============================================================================
# Script: map_rgb_results.R
# Purpose: RGB visualizations of JSDM variance components projected onto
#          Calanda mountain DEM using rayshader
#
# Inputs:
#   - output/starter_data_25.04.25.RData (veg.env, veg)
#   - results_from_Max/res_sjsdm_calanda.rds
#   - output/calanda_mask.shp
#   - output/metrics/dem.tif
#
# Outputs:
#   - plot/ternary_plot_variance_components.pdf
#   - plot/calanda_rgb_3d.png
#   - plot/calanda_rgb_2d_map.pdf
#   - plot/rgb_legend.pdf
#   - plot/calanda_dem_contour_rgb.pdf
# ==============================================================================

library(tidyverse)
library(conflicted)
library(sf)
library(terra)
library(raster)
library(rayshader)
library(ggtern)
library(elevatr)
library(tidyterra)
library(here)

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
# ==============================================================================
# LOAD DATA
# ==============================================================================
cat("\n=== Loading JSDM results ===\n")

# Load JSDM results with variance components
load(here("Calanda_JSDM", "output", "starter_data_25.04.25.RData"))
res = readRDS(here("Calanda_JSDM", "results_from_Max", "res_sjsdm_calanda.rds"))

# Rename dotted objects from RData to snake_case
veg_env = veg.env
rm(veg.env)

# Extract sites with variance components and coordinates
sites =
  res$internals$Sites %>%
  rownames_to_column("plot_id_releve") %>%
  left_join(veg_env %>% as.data.frame() %>% rownames_to_column("plot_id_releve")) %>%
  left_join(veg %>% select(plot_id_releve, open_close) %>% distinct()) %>%
  select(plot_id_releve, env, spa, codist, Longitude, Latitude, open_close)

cat(sprintf("Loaded %d sites with variance components\n", nrow(sites)))

# ==============================================================================
# CREATE RGB COLORS FROM VARIANCE COMPONENTS
# ==============================================================================
cat("\n=== Creating RGB colors ===\n")

# Normalize by row sums (exactly as sjSDM does)
sites = sites %>%
  mutate(
    total = env + spa + codist,
    env_prop = env / total,
    spa_prop = spa / total,
    codist_prop = codist / total
  ) %>%
  filter(!is.nan(codist_prop))

# Create RGB colors (matching your color scheme)
sites = sites %>%
  mutate(
    rgb_color = rgb(
      r = spa_prop,
      g = codist_prop,
      b = env_prop
    )
  )

cat("RGB colors created successfully\n")

# ==============================================================================
# CREATE TERNARY PLOT
# ==============================================================================
cat("\n=== Creating ternary plot ===\n")

# Create ternary plot with shapes for open/close
p_ternary = ggtern(sites, aes(x = env_prop, y = codist_prop, z = spa_prop)) +
  geom_point(aes(color = rgb_color, shape = open_close), size = 3, alpha = 0.7) +
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

# Save ternary plot
pdf(here("Calanda_JSDM", "plot", "ternary_plot_variance_components.pdf"), height = 8, width = 10)
print(p_ternary)
dev.off()

cat("Ternary plot saved to: plot/ternary_plot_variance_components.pdf\n")

# ==============================================================================
# PREPARE DATA FOR RAYSHADER (FOLLOWING WORKING EXAMPLE)
# ==============================================================================
cat("\n=== Preparing data for rayshader ===\n")

# Convert sites to sf object
sites_sf = sites %>%
  select(-Latitude, -Longitude)%>%
  left_join(veg %>% select(plot_id_releve, Longitude, Latitude))%>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

# Load Calanda mask
calanda_mask = st_read(here("Calanda_JSDM", "output", "calanda_mask.shp"))

# Get elevation raster using elevatr (like the working example)
# z parameter controls resolution: higher = better resolution but slower
cat("Downloading elevation data (this may take a moment)...\n")
calanda_elevation_raster = get_elev_raster(calanda_mask, z = 11)

# Convert to matrix using raster package (CRITICAL - this is what works!)
cat("Converting elevation to matrix...\n")
calanda_elevation_matrix = raster_to_matrix(calanda_elevation_raster)

# Calculate shadows
cat("Calculating shadows...\n")
calanda_shadow = ray_shade(calanda_elevation_matrix)
calanda_amb = ambient_shade(calanda_elevation_matrix)

# Get extent and CRS
calanda_extent = ext(calanda_elevation_raster)
calanda_crs = crs(calanda_elevation_raster)

cat("\n=== Diagnostic Information ===\n")
cat("Elevation raster extent:\n")
print(calanda_extent)
cat("\nElevation raster dimensions:", dim(calanda_elevation_raster), "\n")
cat("Elevation raster CRS:", as.character(calanda_crs), "\n")

# Transform coordinates to raster CRS
sites_coords_transformed = st_transform(sites_sf, calanda_crs)
sites_coords_xy = st_coordinates(sites_coords_transformed)

cat("\nOriginal site coordinates (WGS84) range:\n")
cat("  Longitude:", range(sites$Longitude), "\n")
cat("  Latitude:", range(sites$Latitude), "\n")

cat("\nTransformed site coordinates range:\n")
cat("  X:", range(sites_coords_xy[, 1]), "\n")
cat("  Y:", range(sites_coords_xy[, 2]), "\n")

# Extract terrain elevations at plot locations
terrain_elevations = raster::extract(calanda_elevation_raster, sites_coords_transformed)

cat("\nTerrain elevations at plot locations:\n")
cat("  Range:", range(terrain_elevations, na.rm = TRUE), "\n")
cat("  Number of NA values:", sum(is.na(terrain_elevations)), "out of", length(terrain_elevations), "\n")

cat("\nElevation matrix dimensions:", dim(calanda_elevation_matrix), "\n")
cat("Elevation matrix value range:", range(calanda_elevation_matrix, na.rm = TRUE), "\n")

cat("\n=== End Diagnostics ===\n")
cat("Data preparation complete\n")

# ==============================================================================
# CREATE RAYSHADER VISUALIZATION
# ==============================================================================
cat("\n=== Creating 3D visualization ===\n")

if (!file.exists(here("Calanda_JSDM", "plot", "calanda_rgb_3d.png"))) {

  # Create the 3D plot
  calanda_elevation_matrix %>%
    sphere_shade(texture = create_texture("#A7C7E7", "#89CFF0", "#98E4D8", "#AAF0D1", "grey")) %>%
    add_shadow(calanda_shadow, 0.5) %>%
    add_shadow(calanda_amb, 0) %>%
    plot_3d(calanda_elevation_matrix,
            zscale = 10,
            fov = 0,
            theta = 30,
            phi = 45,
            zoom = 0.3,
            windowsize = c(1200, 1000))

  Sys.sleep(0.7)

  # Calculate proper altitude offset (relative to terrain)
  # Elevate points 500m above their terrain elevation
  point_altitude_offset = 50
  point_altitudes = terrain_elevations + point_altitude_offset

  cat(sprintf("Adding %d RGB points to 3D map...\n", nrow(sites)))
  cat("Point altitude range:", range(point_altitudes, na.rm = TRUE), "\n")

  # Add vertical lines from terrain to points FIRST (before points)
  cat("Adding vertical lines...\n")
  for(i in 1:nrow(sites_coords_xy)) {
    if(!is.na(terrain_elevations[i])) {
      render_path(
        extent = calanda_extent,
        lat = c(sites_coords_xy[i, 2], sites_coords_xy[i, 2]),
        long = c(sites_coords_xy[i, 1], sites_coords_xy[i, 1]),
        altitude = c(terrain_elevations[i], point_altitudes[i]),
        zscale = 10,
        color = "black",
        linewidth = 3
      )
    }
  }

  # Now add RGB points at elevated positions
  cat("Adding RGB points...\n")
  render_points(
    extent = calanda_extent,
    lat = sites_coords_xy[, 2],
    long = sites_coords_xy[, 1],
    altitude = point_altitudes,
    zscale = 10,
    color = sites$rgb_color,
    size = 4
  )

  # Add compass
  render_compass(position = "E", compass_radius = 80)

  # Save snapshot
  render_snapshot(here("Calanda_JSDM", "plot", "calanda_rgb_3d.png"), clear = TRUE)

  cat("3D visualization saved to: plot/calanda_rgb_3d.png\n")
} else {
  cat("3D visualization already exists, skipping render\n")
}

# ==============================================================================
# CREATE 2D MAP WITH RGB POINTS
# ==============================================================================
cat("\n=== Creating 2D map ===\n")

# Create 2D map with actual plot points
p_map = ggplot() +
  geom_sf(data = calanda_mask, fill = "grey90", color = "black", linewidth = 1) +
  geom_sf(data = sites_sf, aes(color = rgb_color), size = 1, alpha = 0.8) +
  scale_color_identity() +
  coord_sf() +
  labs(
    title = "Calanda: RGB representation of variance components at survey plots",
    subtitle = "Point colors represent mixing of Environment, Space, and Species associations",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    legend.position = "none"
  )

# Save 2D map
pdf(here("Calanda_JSDM", "plot", "calanda_rgb_2d_map.pdf"), height = 10, width = 12)
print(p_map)
dev.off()

cat("2D map saved to: plot/calanda_rgb_2d_map.pdf\n")

# ==============================================================================
# CREATE LEGEND
# ==============================================================================
cat("\n=== Creating RGB legend ===\n")

# Create a grid showing RGB combinations
legend_grid = expand.grid(
  env = seq(0, 1, length.out = 10),
  spa = seq(0, 1, length.out = 10)
) %>%
  mutate(
    codist = 0.5,  # Fixed at middle value
    rgb_color = rgb(spa, codist, env)
  )

p_legend = ggplot(legend_grid, aes(x = env, y = spa)) +
  geom_tile(aes(fill = rgb_color)) +
  scale_fill_identity() +
  labs(
    x = "Environment contribution",
    y = "Space contribution",
    title = "RGB Color Legend",
    subtitle = "Species associations fixed at 0.5"
  ) +
  theme_bw() +
  theme(text = element_text(size = 12))

# Save legend
pdf(here("Calanda_JSDM", "plot", "rgb_legend.pdf"), height = 6, width = 7)
print(p_legend)
dev.off()

cat("Legend saved to: plot/rgb_legend.pdf\n")

# ==============================================================================
# CREATE 2D CONTOUR PLOT WITH ORIGINAL DEM
# ==============================================================================
cat("\n=== Creating 2D contour plot with DEM ===\n")

# Load original DEM
dem_original = rast(here("Calanda_JSDM", "output", "metrics", "dem.tif"))

# Crop to Calanda mask extent
dem_cropped = crop(dem_original, calanda_mask)

# Aggregate to reduce resolution for faster plotting (factor of 3)
cat("Aggregating DEM for faster plotting...\n")
dem_agg = aggregate(dem_cropped, fact = 3, fun = mean)

# Create contour plot using tidyterra (no data frame conversion needed!)
p_contour = ggplot() +
  tidyterra::geom_spatraster(data = dem_agg, alpha = 0.7) +
  scale_fill_gradient(low = "white", high = "black") +
  tidyterra::geom_spatraster_contour(data = dem_agg,
                                     color = "white", alpha = 0.5, linewidth = 0.3) +
  geom_sf(data = calanda_mask, fill = NA, color = "black", linewidth = 1) +
  geom_sf(data = sites_sf, aes(color = rgb_color, shape = open_close), size = 2, alpha = 0.1) +
  scale_color_identity() +
  coord_sf() +
  labs(
    title = "Calanda: RGB variance components on DEM",
    subtitle = "Contour lines show elevation | Point colors show variance components",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    legend.position = "right"
  )

# Save contour plot
pdf(here("Calanda_JSDM", "plot", "calanda_dem_contour_rgb.pdf"), height = 10, width = 14)
print(p_contour)
dev.off()

cat("DEM contour plot saved to: plot/calanda_dem_contour_rgb.pdf\n")

cat("\n=== RGB mapping completed successfully! ===\n")
cat("\nOutputs created:\n")
cat("  - plot/ternary_plot_variance_components.pdf\n")
cat("  - plot/calanda_rgb_3d.png\n")
cat("  - plot/calanda_rgb_2d_map.pdf\n")
cat("  - plot/rgb_legend.pdf\n")
cat("  - plot/calanda_dem_contour_rgb.pdf\n")
