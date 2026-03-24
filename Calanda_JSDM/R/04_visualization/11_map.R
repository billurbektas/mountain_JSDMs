# ==============================================================================
# Script: 11_map.R
# Author: Billur Bektas
# Claude (Anthropic) was used to assist with code refactoring, validation, and documentation.
#
# Purpose: Study area map with modeled communities colored by climate RGB
#          triangle (summer_temp = Red, et.annual = Green, |fdd| = Blue).
#          Points sized by log(species richness), shaped by open/closed habitat.
#          Sites not in functional analysis shown at alpha = 0.5.
#          Switzerland inset and RGB triangle legend embedded in figure.
#
# Inputs:
#   - output/veg_clim.csv (unscaled environmental data with coordinates)
#   - output/data_calanda_jsdm_2026-03-06.rds (X, Y matrices — defines modeled sites)
#   - output/community_traits_unweighted_imputed.csv (defines functional analysis sites)
#   - output/calanda_mask.shp (study area polygon)
#   - output/metrics/dem.tif (high-resolution DEM)
#   - data/vegetation/2024_CAPHE_SpeDis_CleanData_20240214.csv (open_close habitat)
#
# Outputs:
#   - plot/map_calanda.pdf
#
# Requires: R/00_setup/functions_calanda.R
# ==============================================================================

library(tidyverse)
library(conflicted)
library(sf)
library(terra)
library(tidyterra)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(ggnewscale)
library(here)

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(tidyr::extract)

source(here("Calanda_JSDM", "R", "00_setup", "functions_calanda.R"))

# ==============================================================================
# LOAD DATA
# ==============================================================================
cat("\n=== Building study area map ===\n")

# Environmental data (unscaled)
veg_clim = read_csv(here("Calanda_JSDM", "output", "veg_clim.csv"),
                     show_col_types = FALSE) %>%
  select(-any_of("...1"))

# JSDM data (to identify modeled sites + species richness)
data_calanda = readRDS(here("Calanda_JSDM", "output", "data_calanda_jsdm_2026-03-06.rds"))
Y = data_calanda$Y
modeled_sites = rownames(Y)
sp_richness = rowSums(Y)

# Community traits (defines which modeled sites enter functional analysis)
comm_traits = read_csv(here("Calanda_JSDM", "output",
                             "community_traits_unweighted_imputed.csv"),
                        show_col_types = FALSE)
functional_sites = comm_traits$plot_id_releve

# Open/closed habitat from raw vegetation data
raw_veg = read_csv(here("Calanda_JSDM", "data", "vegetation",
                         "2024_CAPHE_SpeDis_CleanData_20240214.csv"),
                    show_col_types = FALSE) %>%
  mutate(plot_id_releve = paste0(plot_id, releve_id)) %>%
  distinct(plot_id_releve, open_close)

# Study area mask
calanda_mask_original = st_read(here("Calanda_JSDM", "output", "calanda_mask.shp"), quiet = TRUE)
# Slightly buffered version (300m) for display — captures 14 plots just outside the boundary
# To revert: replace calanda_mask with calanda_mask_original below
# Buffer to capture all outlier points, smooth, then clip to DEM extent
dem_extent = st_as_sf(st_as_sfc(st_bbox(rast(here("Calanda_JSDM", "output", "metrics", "dem.tif")))))
calanda_mask = st_buffer(calanda_mask_original, dist = 300) %>%
  st_union() %>%
  smoothr::smooth(method = "ksmooth", smoothness = 2) %>%
  st_as_sf() %>%
  st_intersection(dem_extent)

# DEM
dem = rast(here("Calanda_JSDM", "output", "metrics", "dem.tif"))

# Switzerland boundary
switzerland = ne_countries(scale = "medium", country = "Switzerland", returnclass = "sf")

# ==============================================================================
# BUILD PLOT DATA — only modeled sites
# ==============================================================================
cat("  Building plot data...\n")

# Filter to modeled sites only, add metadata
plot_data = veg_clim %>%
  filter(plot_id_releve %in% modeled_sites) %>%
  left_join(raw_veg, by = "plot_id_releve") %>%
  mutate(
    sp_rich = sp_richness[plot_id_releve],
    in_functional = plot_id_releve %in% functional_sites,
    habitat = case_when(
      open_close == "open" ~ "Open",
      open_close == "close" ~ "Closed",
      TRUE ~ "Unknown"
    )
  )

# --- RGB from climate variables ---
# summer_temp -> Red, et.annual -> Green, |fdd| -> Blue
plot_data = plot_data %>%
  mutate(
    r_val = scales::rescale(summer_temp, to = c(0, 1)),
    g_val = scales::rescale(et.annual, to = c(0, 1)),
    b_val = scales::rescale(abs(fdd), to = c(0, 1)),
    rgb_color = rgb(r_val, g_val, b_val),
    log_sp_rich = log(sp_rich),
    # Bin species richness for discrete size legend
    sp_rich_bin = factor(cut(sp_rich,
      breaks = c(0, 7, 15, 30, Inf),
      labels = c("5", "10", "20", "40"),
      include.lowest = TRUE
    ), levels = c("5", "10", "20", "40")),
    point_alpha = ifelse(in_functional, 0.8, 0.5)
  )

# Convert to sf
plot_sf = plot_data %>%
  filter(!is.na(Latitude) & !is.na(Longitude)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

# Split for layered plotting (functional behind, non-functional in front with lower alpha)
sf_functional = plot_sf %>% filter(in_functional)
sf_nonfunctional = plot_sf %>% filter(!in_functional)

cat(sprintf("  %d modeled sites: %d in functional analysis, %d not (alpha=0.5)\n",
            nrow(plot_sf), nrow(sf_functional), nrow(sf_nonfunctional)))

# ==============================================================================
# DEM PREPARATION
# ==============================================================================
cat("  Preparing DEM...\n")

veg_bbox = st_bbox(plot_sf)
dem_crop = crop(dem, ext(
  veg_bbox["xmin"] - 0.005, veg_bbox["xmax"] + 0.005,
  veg_bbox["ymin"] - 0.005, veg_bbox["ymax"] + 0.005
))
dem_plot = aggregate(dem_crop, fact = 5, fun = "mean")

# Clip the study area boundary to the DEM footprint (non-NA cells)
# so the dashed line never extends beyond the visible relief
dem_footprint = as.polygons(dem_crop > -Inf) %>% st_as_sf() %>% st_union()
calanda_mask = st_intersection(calanda_mask, dem_footprint)

# ==============================================================================
# MAIN MAP: CALANDA DEM + RGB COMMUNITIES
# ==============================================================================
cat("  Creating Calanda DEM map...\n")

# Add an alpha_label column for the legend
sf_nonfunctional$alpha_label = "Not in functional analysis"
sf_functional$alpha_label = "In functional analysis"

p_main = ggplot() +
  # DEM greyscale relief — with altitude legend
  geom_spatraster(data = dem_plot) +
  scale_fill_gradientn(
    colors = grey.colors(10, start = 0.3, end = 1.0),
    na.value = "transparent",
    name = "Altitude (m)",
    guide = guide_colorbar(order = 4, barwidth = 1, barheight = 6)
  ) +
  new_scale_fill() +
  # Study area boundary removed — DEM + points define the area
  # Non-functional sites (alpha = 0.5, behind)
  geom_sf(data = sf_nonfunctional,
          aes(shape = habitat, size = sp_rich_bin, alpha = alpha_label),
          fill = sf_nonfunctional$rgb_color,
          color = "black", stroke = 0.3) +
  # Functional sites (alpha = 0.8, on top)
  geom_sf(data = sf_functional,
          aes(shape = habitat, size = sp_rich_bin, alpha = alpha_label),
          fill = sf_functional$rgb_color,
          color = "black", stroke = 0.3) +
  # Alpha scale: functional vs non-functional
  scale_alpha_manual(
    values = c("In functional analysis" = 0.8, "Not in functional analysis" = 0.5),
    name = "Analysis",
    guide = guide_legend(order = 3,
                         override.aes = list(fill = "grey50", shape = 21, size = 3))
  ) +
  # Shapes: fillable (21 = circle, 24 = triangle)
  scale_shape_manual(
    values = c("Open" = 21, "Closed" = 24, "Unknown" = 22),
    name = "Habitat",
    guide = guide_legend(order = 1,
                         override.aes = list(fill = "grey50", size = 3, alpha = 0.8))
  ) +
  # Species richness size with explicit breaks
  scale_size_manual(
    values = c("5" = 0.5, "10" = 1.2, "20" = 2.0, "40" = 3.5),
    name = "Species richness",
    guide = guide_legend(order = 2,
                         override.aes = list(fill = "grey50", shape = 21, alpha = 0.8))
  ) +
  # Coordinate limits
  coord_sf(
    xlim = c(veg_bbox["xmin"] - 0.003, veg_bbox["xmax"] + 0.003),
    ylim = c(veg_bbox["ymin"] - 0.003, veg_bbox["ymax"] + 0.003)
  ) +
  # North arrow
  annotation_north_arrow(
    location = "tr", which_north = "true",
    height = unit(1.2, "cm"), width = unit(1.0, "cm"),
    style = north_arrow_fancy_orienteering(text_size = 8)
  ) +
  # Scale bar
  annotation_scale(
    location = "bl", width_hint = 0.2,
    style = "ticks", text_cex = 0.8,
    line_width = 0.5, height = unit(0.2, "cm")
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 8),
    legend.position = "right",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold"),
    panel.grid = element_line(color = "grey90", linewidth = 0.2)
  )

# ==============================================================================
# SWITZERLAND INSET (drawn as grob, placed inside the main plot)
# ==============================================================================
cat("  Creating Switzerland inset...\n")

p_swiss = ggplot() +
  geom_sf(data = switzerland, fill = "grey95", color = "grey40", linewidth = 0.3) +
  geom_sf(data = calanda_mask, fill = "grey40", alpha = 0.8, color = NA) +
  annotate("text",
           x = st_coordinates(st_centroid(st_union(calanda_mask)))[1] + 0.35,
           y = st_coordinates(st_centroid(st_union(calanda_mask)))[2] + 0.18,
           label = "Calanda", size = 2.8, fontface = "bold", color = "grey20") +
  theme_void() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(0, 0, 0, 0)
  )

# ==============================================================================
# RGB TRIANGLE LEGEND
# ==============================================================================
cat("  Creating RGB triangle legend...\n")

n_grid = 50
tri_grid = expand_grid(
  i = seq(0, n_grid),
  j = seq(0, n_grid)
) %>%
  filter(i + j <= n_grid) %>%
  mutate(
    k = n_grid - i - j,
    r = i / n_grid,
    g = j / n_grid,
    b = k / n_grid,
    x = 0.5 * (2 * g + b) / (r + g + b + 1e-10),
    y = (sqrt(3) / 2) * b / (r + g + b + 1e-10),
    color = rgb(r, g, b)
  )

p_triangle = ggplot() +
  geom_point(data = tri_grid, aes(x = x, y = y), color = tri_grid$color,
             size = 2.5, shape = 15) +
  # Vertex labels — placed further out to avoid cutoff
  # Labels outside the triangle at each vertex
  annotate("label", x = 0, y = -0.08,
           label = "Summer Temp.",
           size = 2.3, fontface = "bold", color = "red3", hjust = 0.5,
           fill = "white", alpha = 0.8, label.size = 0, label.padding = unit(0.1, "lines")) +
  annotate("label", x = 1, y = -0.08,
           label = "Annual ET",
           size = 2.3, fontface = "bold", color = "green4", hjust = 0.5,
           fill = "white", alpha = 0.8, label.size = 0, label.padding = unit(0.1, "lines")) +
  annotate("label", x = 0.5, y = sqrt(3)/2 + 0.08,
           label = "|FDD|",
           size = 2.3, fontface = "bold", color = "blue3", hjust = 0.5,
           fill = "white", alpha = 0.8, label.size = 0, label.padding = unit(0.1, "lines")) +
  # Triangle outline
  annotate("segment", x = 0, xend = 1, y = 0, yend = 0, linewidth = 0.4) +
  annotate("segment", x = 1, xend = 0.5, y = 0, yend = sqrt(3)/2, linewidth = 0.4) +
  annotate("segment", x = 0.5, xend = 0, y = sqrt(3)/2, yend = 0, linewidth = 0.4) +
  coord_fixed(clip = "off") +
  theme_void() +
  theme(plot.margin = margin(10, 8, 10, 8))

# ==============================================================================
# COMPOSE: inset CH map and RGB triangle inside the main plot
# ==============================================================================
cat("  Composing final figure...\n")

# Convert insets to grobs
swiss_grob = ggplotGrob(p_swiss)
triangle_grob = ggplotGrob(p_triangle)

# Add insets to main plot using annotation_custom
p_final = p_main +
  # Switzerland inset — top-left corner, pushed slightly up
  annotation_custom(
    grob = swiss_grob,
    xmin = veg_bbox["xmin"] - 0.003,
    xmax = veg_bbox["xmin"] + 0.048,
    ymin = veg_bbox["ymax"] - 0.035,
    ymax = veg_bbox["ymax"] + 0.008
  ) +
  # RGB triangle — bottom-right corner, shifted slightly left
  annotation_custom(
    grob = triangle_grob,
    xmin = veg_bbox["xmax"] - 0.045,
    xmax = veg_bbox["xmax"] - 0.004,
    ymin = veg_bbox["ymin"] - 0.003,
    ymax = veg_bbox["ymin"] + 0.038
  )

pdf(here("Calanda_JSDM", "plot", "map_calanda.pdf"), width = 11, height = 9)
print(p_final)
dev.off()

cat("  Saved plot/map_calanda.pdf\n")
cat("\n=== Map complete ===\n")
