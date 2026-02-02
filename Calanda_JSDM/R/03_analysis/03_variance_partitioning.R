# ==============================================================================
# Script: 03_variance_partitioning.R
# Purpose: Variance partitioning visualization (Venn diagram + ternary plots)
#
# Inputs:
#   - output/starter_data_25.04.25.RData
#   - results_from_Max/an_sjsdm_calanda.rds
#   - results_from_Max/model_sjsdm_calanda.rds
#   - results_from_Max/res_sjsdm_calanda.rds
#
# Outputs:
#   - plot/venn.pdf
#   - plot/ternary_sites.pdf
#   - plot/ternary_species.pdf
#
# Requires: R/00_setup/functions_calanda.R (plot_tern_sites, plot_tern_species,
#           plot_anova_custom)
# ==============================================================================

library(tidyverse)
library(sjSDM)
library(here)

source(here("Calanda_JSDM", "R", "00_setup", "functions_calanda.R"))

# Load data ----
load(here("Calanda_JSDM", "output", "starter_data_25.04.25.RData"))
# Rename dotted objects from RData to snake_case
veg_clim = veg.clim
rm(veg.clim)
an = readRDS(here("Calanda_JSDM", "results_from_Max", "an_sjsdm_calanda.rds"))
model_jsdm = readRDS(here("Calanda_JSDM", "results_from_Max", "model_sjsdm_calanda.rds"))
environment(model_jsdm$get_model)$device = "cpu"
res = readRDS(here("Calanda_JSDM", "results_from_Max", "res_sjsdm_calanda.rds"))

# ==============================================================================
# VARIANCE PARTITIONING: VENN DIAGRAM + TERNARY PLOTS
# ==============================================================================
cat("\n=== Creating variance partitioning visualizations ===\n")

# Define custom colors
color_codist = "#00bd89"  # Green for species associations
color_spa = "#d00000"     # Red for space
color_env = "#81caf3"     # Blue for environment

# Venn diagram ----
pdf(file = "plot/venn.pdf", height = 8, width = 8)
p_venn = plot_anova_custom(an)
dev.off()

# Species altitude ranges ----
cat("Calculating species altitude ranges...\n")
species_altitude_ranges = veg_clim %>%
  select(plot_id_releve, altitude) %>%
  distinct() %>%
  left_join(
    veg %>% select(plot_id_releve, species, species_cover),
    by = "plot_id_releve"
  ) %>%
  filter(species_cover > 0) %>%
  group_by(species) %>%
  summarize(
    min_altitude = min(altitude, na.rm = TRUE),
    max_altitude = max(altitude, na.rm = TRUE),
    altitude_range = max_altitude - min_altitude,
    mean_altitude = mean(altitude, na.rm = TRUE),
    .groups = "drop"
  )

# Ternary plots ----
p1 = plot_tern_sites(
  res = res,
  veg_clim = veg_clim,
  color_env = color_env,
  color_codist = color_codist,
  color_spa = color_spa
)

p2 = plot_tern_species(
  res = res,
  veg = veg,
  veg_clim = veg_clim,
  color_env = color_env,
  color_codist = color_codist,
  color_spa = color_spa
)

cat("Combined variance partitioning plot created\n")

pdf(file = "plot/ternary_sites.pdf", height = 10, width = 12)
p1
dev.off()

pdf(file = "plot/ternary_species.pdf", height = 10, width = 12)
p2
dev.off()

cat("\n=== Variance partitioning complete ===\n")
