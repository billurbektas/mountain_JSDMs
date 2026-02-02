# Load required libraries
library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)

# Get Switzerland map data
switzerland <- ne_countries(scale = "medium", country = "Switzerland", returnclass = "sf")

# Load vegetation survey coordinates from CSV
#veg_coord <- read_csv("data/vegetation/veg.coord.csv")

# Convert vegetation coordinates to sf object
# Note: The x and y coordinates appear to be in EPSG:4326 (WGS84) format
veg_points <- veg.clim %>%
  filter(trees_cover == 0)%>%
  st_as_sf(coords = c( "Longitude", "Latitude"), crs = 4326)

# Load your shapefile
my_shapefile <- st_read("output/calanda_mask.shp")

# Basic plot of Switzerland with your shapefile
main_map <- ggplot() +
  geom_sf(data = switzerland, fill = "white", color = "gray70") +
  geom_sf(data = my_shapefile, fill = "grey50", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Switzerland with Study Area")

# Create a zoomed plot for the inset
zoomed_plot <- ggplot() +
                            geom_sf(data = switzerland, fill = "white", color = "gray70") +
                            geom_sf(data = my_shapefile, fill = "grey50", alpha = 0.5, color = "black") +
                            geom_sf(data = veg_points, size = 1, aes(color = altitude)) +
                            coord_sf(xlim = c(veg_bbox["xmin"], veg_bbox["xmax"]),
                                     ylim = c(veg_bbox["ymin"], veg_bbox["ymax"])) +
                            theme_minimal() +
                            scale_color_viridis_c(option = "turbo")+
                            theme(plot.background = element_rect(fill = "white", color = "black"),
                                  panel.background = element_rect(fill = "white")) +
                            labs(title = "Vegetation Surveys in Study Area",
                                 subtitle = "Zoomed view of study site")


# Print the inset map
pdf("plot/map_calanda.pdf", height = 6, width = 6)
print(main_map)
print(zoomed_plot)
dev.off()
