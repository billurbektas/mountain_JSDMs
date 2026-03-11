# ==============================================================================
# Script: map.R
# Purpose: Create study area map of Calanda with vegetation survey points
#
# Inputs:
#   - output/starter_data_25.04.25.RData (veg.clim)
#   - output/calanda_mask.shp
#
# Outputs:
#   - plot/map_calanda.pdf
# ==============================================================================

library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(here)

# Load data ----
load(here("Calanda_JSDM", "output", "starter_data_25.04.25.RData"))
veg_clim = veg.clim
rm(veg.clim)

# Get Switzerland map data
switzerland = ne_countries(scale = "medium", country = "Switzerland", returnclass = "sf")

# Convert vegetation coordinates to sf object
veg_points = veg_clim %>%
  filter(trees_cover == 0) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

# Compute bounding box from vegetation points
veg_bbox = st_bbox(veg_points)

# Load study area shapefile
calanda_mask = st_read(here("Calanda_JSDM", "output", "calanda_mask.shp"))

# Overview map of Switzerland with study area
main_map = ggplot() +
  geom_sf(data = switzerland, fill = "white", color = "gray70") +
  geom_sf(data = calanda_mask, fill = "grey50", alpha = 0.7) +
  theme_bw() +
  labs(title = "Switzerland with Study Area")

# Zoomed map of the study area
zoomed_plot = ggplot() +
  geom_sf(data = switzerland, fill = "white", color = "gray70") +
  geom_sf(data = calanda_mask, fill = "grey50", alpha = 0.5, color = "black") +
  geom_sf(data = veg_points, size = 1, aes(color = altitude)) +
  coord_sf(xlim = c(veg_bbox["xmin"], veg_bbox["xmax"]),
           ylim = c(veg_bbox["ymin"], veg_bbox["ymax"])) +
  scale_color_viridis_c(option = "turbo") +
  theme_bw() +
  theme(plot.background = element_rect(fill = "white", color = "black"),
        panel.background = element_rect(fill = "white")) +
  labs(title = "Vegetation Surveys in Study Area",
       subtitle = "Zoomed view of study site")

pdf(here("Calanda_JSDM", "plot", "map_calanda.pdf"), height = 6, width = 6)
print(main_map)
print(zoomed_plot)
dev.off()

cat("\n=== Map created ===\n")
