# ==============================================================================
# JSDM PREDICTIONS FOR TRANSPLANT EXPERIMENTS (CALANDA)
# ==============================================================================
# This script generates predictions from the fitted JSDM for transplant experiment
# communities (CH_Calanda, CH_Calanda2) under different conditioning scenarios
# to evaluate the role of biotic interactions in novel community assembly.

# Load required libraries ----
library(tidyverse)
library(sjSDM)
library(vegan)
library(boot)
library(patchwork)
library(here)
setwd(here("Calanda_JSDM"))

# ==============================================================================
# 1. SETUP & LOAD DATA
# ==============================================================================
cat("\n=== Loading data and models ===\n")

# Load starter data (contains veg, veg.env, etc.)
load("output/starter_data_25.04.25.RData")
load("data/vegetation/TransPlantNetwork_101024.RData")

# Load fitted JSDM model
model_sjsdm = readRDS("results_from_Max/model_sjsdm_calanda.rds")
environment(model_sjsdm$get_model)$device = "cpu"

#model_sjsdm = readRDS("output/model_sjsdm_calanda_260425.rds")
cat("JSDM model loaded successfully\n")

# Load transplant experiment data
# TransPlant data structure:
# - 'dat' is the global dataframe with ALL vegetation data across experiments
# - Each experiment object (CH_Calanda, etc.) contains only metadata in $meta
# - We filter 'dat' by Region to get experiment-specific data

# Check if global data exists
if(!exists("dat")) {
  stop("Global 'dat' dataframe not found. Please load TransPlantNetwork.RData first:\n  load('path/to/TransPlantNetwork.RData')")
}

cat("Transplant experiment data found\n")

# ==============================================================================
# 2. ORGANIZE TRANSPLANT COMMUNITIES
# ==============================================================================
cat("\n=== Organizing transplant communities ===\n")

# Following the TransPlant protocol from 02_organize.R lines 96-108
# Key concepts:
# - originSiteID = site where turf came from (HIGH elevation)
# - destSiteID = site where turf is now
# - Treatment == "LocalControl": not transplanted (stayed at origin/dest)
# - Treatment == "Warm": transplanted from high to low

# Filter to Calanda experiments only
calanda_regions = c("CH_Calanda")

calanda_dat = dat %>%
  filter(Region %in% calanda_regions) %>%
  # Re-calculate relative cover (following TransPlant protocol)
  group_by(Region, Year, originSiteID, originBlockID, destSiteID, destBlockID,
           destPlotID, Treatment, turfID, UniqueID) %>%
  mutate(Total_Cover = sum(Cover, na.rm = TRUE),
         Rel_Cover = Cover/Total_Cover) %>%
  ungroup() %>%
  filter(!is.na(Rel_Cover)) %>%
  filter(!is.na(SpeciesName))

cat(sprintf("Raw Calanda data: %d records\n", nrow(calanda_dat)))

# Following EXACT TransPlant protocol (lines 96-108 of 02_organize.R)
# Get distinct warmed treatment paths (Region + originSiteID + destSiteID)
warmed_paths = calanda_dat %>%
  filter(Treatment == "Warm") %>%
  distinct(Region, originSiteID, destSiteID)

cat("\nWarmed treatment paths found:\n")
print(warmed_paths)

# For each warmed path, extract corresponding control and warmed plots
# Using pmap to process each path separately (prevents Cartesian products)
calanda_nested = warmed_paths %>%
  mutate(comm = pmap(
    .l = list(R = Region, O = originSiteID, D = destSiteID),
    .f = function(R, O, D) {
      bind_rows(
        # originControls (NH): stayed at HIGH/origin site
        originControls = calanda_dat %>%
          filter(Region == R &
                   originSiteID == O &
                   destSiteID == O &
                   Treatment == "LocalControl"),

        # destControls (NL): stayed at LOW/destination site
        destControls = calanda_dat %>%
          filter(Region == R &
                   originSiteID == D &
                   destSiteID == D &
                   Treatment == "LocalControl"),

        # warmed (W): transplanted from HIGH to LOW
        warmed = calanda_dat %>%
          filter(Region == R &
                   originSiteID == O &
                   destSiteID == D &
                   Treatment == "Warm"),

        .id = "ODT"
      )
    }
  ))

# Unnest and rename plot types
calanda_data = calanda_nested %>%
  filter(originSiteID == "Cal")%>%
  select(comm) %>%

  unnest(comm) %>%
  mutate(
    plot_type = case_when(
      ODT == "originControls" ~ "high_elevation",
      ODT == "destControls" ~ "low_elevation",
      ODT == "warmed" ~ "warmed"
    ),
    experiment = Region  # Add experiment column for later use
  ) %>%
  select(-ODT)

cat(sprintf("\nOrganized Calanda data: %d records\n", nrow(calanda_data)))
cat(sprintf("  High elevation plots: %d\n", sum(calanda_data$plot_type == "high_elevation")))
cat(sprintf("  Low elevation plots: %d\n", sum(calanda_data$plot_type == "low_elevation")))
cat(sprintf("  Warmed plots: %d\n", sum(calanda_data$plot_type == "warmed")))

# Convert to presence/absence and create wide format
# Create unique plot identifier
calanda_data = calanda_data %>%
  mutate(plot_id = paste(experiment, destSiteID, destBlockID, destPlotID, Year, sep = "_"))

# Convert to wide format (plot × species)
comm_wide = calanda_data %>%
  select(plot_id, plot_type, SpeciesName, Cover, experiment,
         originSiteID, destSiteID, Year) %>%
  # Convert to presence/absence
  mutate(presence = ifelse(Cover > 0, 1, 0)) %>%
  select(-Cover) %>%
  distinct() %>%
  pivot_wider(names_from = SpeciesName,
              values_from = presence,
              values_fill = 0)

# Extract plot metadata
plot_meta = comm_wide %>%
  select(plot_id, plot_type, experiment, originSiteID, destSiteID, Year)

# Extract species matrix
Y_all = comm_wide %>%
  select(-plot_id, -plot_type, -experiment, -originSiteID, -destSiteID, -Year) %>%
  as.matrix()

rownames(Y_all) = comm_wide$plot_id

cat(sprintf("Community matrix: %d plots × %d species\n", nrow(Y_all), ncol(Y_all)))

# Filter species: keep only those occurring ≥5 times
species_counts = colSums(Y_all)
keep_species = names(species_counts)[species_counts >= 5]

Y_all = Y_all[, keep_species]
cat(sprintf("After filtering (≥5 occurrences): %d species retained\n", ncol(Y_all)))

# Split into treatment groups
Y_high = Y_all[plot_meta$plot_type == "high_elevation", ]
Y_low = Y_all[plot_meta$plot_type == "low_elevation", ]
Y_warm = Y_all[plot_meta$plot_type == "warmed", ]

cat(sprintf("\nTreatment groups:\n"))
cat(sprintf("  Y_NH: %d plots × %d species\n", nrow(Y_high), ncol(Y_high)))
cat(sprintf("  Y_NL: %d plots × %d species\n", nrow(Y_low), ncol(Y_low)))
cat(sprintf("  Y_W: %d plots × %d species\n", nrow(Y_warm), ncol(Y_warm)))

# ==============================================================================
# 3. CLASSIFY SPECIES POOLS
# ==============================================================================
cat("\n=== Classifying species pools ===\n")

# Get species lists from each treatment
sp_high = colnames(Y_high)[colSums(Y_high) > 0]
sp_low = colnames(Y_low)[colSums(Y_low) > 0]
sp_warm = colnames(Y_warm)[colSums(Y_warm) > 0]

# Classify following TransPlant logic
high_unique = setdiff(sp_high, sp_low)  # Strictly high elevation
overlap = intersect(sp_high, sp_low)     # Found in both
low_only = setdiff(sp_low, sp_high)     # Strictly low elevation

# For warmed plots
warmed_residents = intersect(sp_warm, c(high_unique, overlap))  # R_high
warmed_colonizers = setdiff(sp_warm, warmed_residents)          # C_low (includes outsiders)

# Create species pool classification dataframe
species_pools = data.frame(
  species = colnames(Y_all),
  stringsAsFactors = FALSE
) %>%
  mutate(
    pool = case_when(
      species %in% high_unique ~ "high_unique",
      species %in% overlap ~ "overlap",
      species %in% low_only ~ "low_only",
      TRUE ~ "other"
    ),
    is_resident = species %in% warmed_residents,
    is_colonizer = species %in% warmed_colonizers,
    in_warmed = species %in% sp_warm
  )

cat(sprintf("Species pool classification:\n"))
cat(sprintf("  High unique: %d\n", sum(species_pools$pool == "high_unique")))
cat(sprintf("  Overlap: %d\n", sum(species_pools$pool == "overlap")))
cat(sprintf("  Low only: %d\n", sum(species_pools$pool == "low_only")))
cat(sprintf("\nIn warmed plots:\n"))
cat(sprintf("  Residents (R_high): %d\n", sum(species_pools$is_resident)))
cat(sprintf("  Colonizers (C_low): %d\n", sum(species_pools$is_colonizer)))

# ==============================================================================
# 4. EXTRACT ENVIRONMENTAL DATA
# ==============================================================================
cat("\n=== Extracting environmental data ===\n")

# Get elevation information from veg data
# Match plot IDs to get elevations
# First, we need to map transplant plots to environmental data

# For this example, I'll use the site-level averages
# You may need to adjust this based on your actual data structure

high_site_elevation = 2000
low_site_elevation = 1400
longitude_high = CH_Calanda$meta$Longitude[CH_Calanda$meta$destSiteID == "Cal"]
longitude_low = CH_Calanda$meta$Longitude[CH_Calanda$meta$destSiteID == "Nes"]
latitude_high = CH_Calanda$meta$Latitude[CH_Calanda$meta$destSiteID == "Cal"]
latitude_low = CH_Calanda$meta$Latitude[CH_Calanda$meta$destSiteID == "Nes"]

# Get plots within ±100m elevation
high_buffer_plots = veg %>%
  filter(abs(altitude - high_site_elevation) <= 100) %>%
  pull(plot_id_releve)%>%
  unique()

low_buffer_plots = veg %>%
  filter(abs(altitude - low_site_elevation) <= 100) %>%
  pull(plot_id_releve)%>%
  unique()


# Extract environmental variables for high and low
# Make sure these plots exist in veg.env
high_buffer_plots = intersect(high_buffer_plots, rownames(veg.env))
low_buffer_plots = intersect(low_buffer_plots, rownames(veg.env))

if(length(high_buffer_plots) > 0 & length(low_buffer_plots) > 0) {
  X_high = colMeans(veg.env[high_buffer_plots, , drop=FALSE], na.rm = TRUE)
  X_low = colMeans(veg.env[low_buffer_plots, , drop=FALSE], na.rm = TRUE)

  cat(sprintf("X_high: averaged over %d plots\n", length(high_buffer_plots)))
  cat(sprintf("X_low: averaged over %d plots\n", length(low_buffer_plots)))
} else {
  stop("Could not find matching plots in veg.env for elevation buffer")
}

# Create environmental matrices for predictions
# Each plot gets the same environment within treatment
n_high = nrow(Y_high)
n_low = nrow(Y_low)
n_warm = nrow(Y_warm)

env_high = matrix(rep(X_high, n_high), nrow = n_high, byrow = TRUE)
env_low = matrix(rep(X_low, n_low), nrow = n_low, byrow = TRUE)
env_warm = matrix(rep(X_low, n_warm), nrow = n_warm, byrow = TRUE)

colnames(env_high) = colnames(env_low) = colnames(env_warm) = names(X_high)

cat("Environmental matrices created\n")

# ==============================================================================
# 5. ALIGN SPECIES BETWEEN JSDM AND TRANSPLANT DATA
# ==============================================================================
cat("\n=== Aligning species between JSDM and transplant data ===\n")

# Get species from JSDM model
# This depends on how sjSDM stores species names
# Typically in model$species or similar
jsdm_species = model_sjsdm$species  # Y is the species matrix used to fit the model

# Find common species
common_species = intersect(colnames(Y_all), jsdm_species)
cat(sprintf("Common species between JSDM and transplant: %d\n", length(common_species)))

# Subset to common species
Y_high = Y_high[, common_species]
Y_low = Y_low[, common_species]
Y_warm = Y_warm[, common_species]

# Update species pool classifications
species_pools = species_pools %>%
  filter(species %in% common_species)

warmed_residents = species_pools %>%
  filter(is_resident) %>%
  pull(species)

warmed_colonizers = species_pools %>%
  filter(is_colonizer) %>%
  pull(species)

cat(sprintf("After alignment:\n"))
cat(sprintf("  Total species: %d\n", length(common_species)))
cat(sprintf("  Residents: %d\n", length(warmed_residents)))
cat(sprintf("  Colonizers: %d\n", length(warmed_colonizers)))

# ==============================================================================
# 6. PREDICTION FUNCTIONS
# ==============================================================================
cat("\n=== Setting up prediction functions ===\n")

# Function for marginal prediction (environment only)
predict_marginal = function(model, env_data, Y_obs, latitude, longitude) {
  SP_data = data.frame(Latitude = rep(latitude, nrow(env_data)),
                       Longitude = rep(longitude, nrow(env_data)))
  # Use sjSDM predict function with environment only
  pred = predict(object = model, newdata = env_data, SP = SP_data, type = "link")
  colnames(pred) = model$species
  pred = pred[, colnames(Y_obs)]  # This is problematic maybe bc it should predict the absences too.
  return(pred)
}

# Function for conditional prediction (environment + codistribution)
predict_conditional = function(model, env_data, Y_obs, latitude, longitude) {
  SP_data = data.frame(Latitude = rep(latitude, nrow(env_data)),
                       Longitude = rep(longitude, nrow(env_data)))
  
  # Pad Y_obs to match training data dimensions
  missing_species = setdiff(model$species, colnames(Y_obs))
  if(length(missing_species) > 0) {
    zeros = matrix(0, nrow(Y_obs), length(missing_species))
    colnames(zeros) = missing_species
    Y_full = cbind(Y_obs, zeros)[, model$species]
  } else {
    Y_full = Y_obs[, model$species]
  }
  
  # Initialize result matrix
  result = matrix(NA, nrow(Y_obs), ncol(Y_obs))
  colnames(result) = colnames(Y_obs)
  result = as.data.frame(result)
  
  # Predict each observed species conditional on others
  for(sp in colnames(Y_obs)) {
    Y_cond = Y_full
    Y_cond[, sp] = NA  # Set target species to NA
    
    pred = predict(model, newdata = env_data, Y = Y_cond, SP = SP_data, type = "link")
    pred = as.data.frame(pred)
    colnames(pred) = colnames(Y_cond[,which(is.na(colSums(Y_cond))), drop = FALSE])
    result[,sp] = pred[,sp, drop = FALSE]  # Extract prediction for target species
  }
  
  return(result)
}
# Note: The above conditional prediction may need adjustment based on sjSDM's API
# sjSDM may not have built-in conditional prediction
# Alternative approach: manually implement using model coefficients

cat("Prediction functions defined\n")
cat("WARNING: Conditional prediction may need adjustment based on sjSDM API\n")

# ==============================================================================
# 7. RUN PREDICTION SCENARIOS ----
# ==============================================================================
cat("\n=== Running prediction scenarios ===\n")

# Initialize results list
predictions = list()

# First, get training data predictions (S0a, S0b)
cat("\n=== Training Data Predictions ===\n")

# Load training data
Y_train = Y[, model_sjsdm$species]

# Scenario 0a: Training Marginal
cat("\nScenario 0a: Training - Marginal\n")
pred_train_marginal = predict(model_sjsdm,
                               newdata = X[,model_sjsdm$names[-1]],
                               SP = X[, c("Latitude", "Longitude")],
                               type = "link")
colnames(pred_train_marginal) = model_sjsdm$species
predictions$s0a_train_marginal = pred_train_marginal

# Scenario 0b: Training Conditional (leave-one-out)
cat("\nScenario 0b: Training - Conditional (leave-one-out)\n")
# For conditional, we predict each species given others present
pred_train_conditional = matrix(NA, nrow(Y_train), ncol(Y_train))
colnames(pred_train_conditional) = colnames(Y_train)

for(j in 1:ncol(pred_train_conditional)) {
  sp = colnames(pred_train_conditional)[j]
  cat(sprintf("  Predicting %s (%d/%d)\r", sp, j, ncol(pred_train_conditional)))

  # Set target species to NA
  Y_cond = Y_train
  Y_cond[, sp] = NA

  # Predict
  pred = predict(model_sjsdm,
                 newdata = X[,model_sjsdm$names[-1]],
                 Y = Y_cond,
                 SP = X[, c("Latitude", "Longitude")],
                 type = "link")
  colnames(pred) =  sp
  # Extract prediction for target species
  pred_train_conditional[, j] = pred[, sp]
}
cat("\n")

predictions$s0b_train_conditional = pred_train_conditional

cat("Training predictions complete\n")

# Now run transplant experiment predictions

# Scenario 1: High elevation - Marginal
cat("\nScenario 1: High elevation - Marginal\n")
predictions$s1_high_elevation_marginal = predict_marginal(model = model_sjsdm, env_data = env_high, Y_obs = Y_high, latitude = latitude_high, longitude = longitude_high)

# Scenario 2: Native High - Conditional
cat("\nScenario 2: NH - Conditional (all others)\n")
predictions$s2_high_elevation_conditional = predict_conditional(model = model_sjsdm, env_data = env_high, Y_obs = Y_high, 
                                                                latitude = latitude_high, longitude = longitude_high)

# Scenario 3: Native Low - Marginal
cat("\nScenario 3: NL - Marginal\n")
predictions$s3_low_elevation_marginal = predict_marginal(model = model_sjsdm, env_data = env_low, Y_obs = Y_low, latitude = latitude_low, longitude = longitude_low)

# Scenario 4: Native Low - Conditional
cat("\nScenario 4: NL - Conditional (all others)\n")
predictions$s4_low_elevation_conditional = predict_conditional(model = model_sjsdm, env_data = env_low, Y_obs = Y_low, 
                                                               latitude = latitude_low, longitude = longitude_low)

# Scenario 5: Warmed - Marginal
cat("\nScenario 5: Warmed - Marginal\n")
predictions$s5_warmed_marginal = predict_marginal(model = model_sjsdm, env_data = env_warm, Y_obs = Y_warm, latitude = latitude_low, longitude = longitude_low)

# Scenario 6: Warmed - Residents Only
cat("\nScenario 6: Warmed - Residents only\n")
# Create Y with only residents observed, colonizers set to NA
Y_warm_residents = Y_warm
Y_warm_residents[, warmed_colonizers] = NA
predictions$s6_warmed_residents = predict_conditional(model = model_sjsdm, env_data = env_warm, Y_obs = Y_warm_residents, 
                                                      latitude = latitude_low, longitude = longitude_low)

# Scenario 7: Warmed - Colonizers Only
cat("\nScenario 7: Warmed - Colonizers only\n")
# Create Y with only colonizers observed, residents set to NA
Y_warm_colonizers = Y_warm
Y_warm_colonizers[,warmed_residents]= NA
predictions$s7_warmed_colonizers = predict_conditional(model = model_sjsdm, env_data = env_warm, Y_obs = Y_warm_colonizers, 
                                                        latitude = latitude_low, longitude = longitude_low)

# Scenario 8: Warmed - All Others
cat("\nScenario 8: Warmed - All others\n")
predictions$s8_warmed_all = predict_conditional(model = model_sjsdm, env_data = env_warm, Y_obs = Y_warm, 
                                                       latitude = latitude_low, longitude = longitude_low)

cat("\nAll scenarios completed\n")

# ==============================================================================
# 8. CALCULATE METRICS
# ==============================================================================
cat("\n=== Calculating evaluation metrics ===\n")

# Function to calculate log-loss
calc_logloss = function(y_obs, y_pred) {
  # Clip probabilities to avoid log(0)
  y_pred = pmax(pmin(y_pred, 1 - 1e-15), 1e-15)

  # Calculate log-loss
  logloss = -(y_obs * log(y_pred) + (1 - y_obs) * log(1 - y_pred))

  return(logloss)
}

# Function to calculate Bray-Curtis dissimilarity
calc_bray_curtis = function(y_obs, y_pred) {
  # For each plot (row), calculate BC dissimilarity
  bc = numeric(nrow(y_obs))

  for(i in 1:nrow(y_obs)) {
    # Using probability-based BC
    bc[i] = sum(abs(y_obs[i, ] - y_pred[i, ])) / sum(y_obs[i, ] + y_pred[i, ])
  }

  return(bc)
}

# Calculate metrics for each scenario
metrics = list()

# Helper function to compute metrics for a scenario
compute_metrics = function(y_obs, y_pred, scenario_name) {

  # Log-loss per species per plot
  logloss_matrix = calc_logloss(y_obs, y_pred)

  # Aggregate to species-level (mean across plots)
  logloss_species = colMeans(logloss_matrix, na.rm = TRUE)

  # Aggregate to plot-level (mean across species)
  logloss_plot = rowMeans(logloss_matrix, na.rm = TRUE)

  # Bray-Curtis per plot
  bc_plot = calc_bray_curtis(y_obs, y_pred)

  return(list(
    scenario = scenario_name,
    logloss_matrix = logloss_matrix,
    logloss_species = logloss_species,
    logloss_plot = logloss_plot,
    bc_plot = bc_plot
  ))
}

# Compute metrics for all scenarios
cat("Computing metrics for all scenarios...\n")

# Training data metrics (full dataset)
metrics$s0a_full = compute_metrics(Y[, model_sjsdm$species], predictions$s0a_train_marginal, "S0a: training marginal (full)")
metrics$s0b_full = compute_metrics(Y[, model_sjsdm$species], predictions$s0b_train_conditional, "S0b: training conditional (full)")

# Function to identify outliers using IQR method
identify_outliers = function(x, k = 1.5) {
  q1 = quantile(x, 0.25, na.rm = TRUE)
  q3 = quantile(x, 0.75, na.rm = TRUE)
  iqr = q3 - q1
  lower_bound = q1 - k * iqr
  upper_bound = q3 + k * iqr
  outliers = x < lower_bound | x > upper_bound
  return(outliers)
}

# Identify outliers in training log-loss
outliers_s0a = identify_outliers(metrics$s0a_full$logloss_plot)
outliers_s0b = identify_outliers(metrics$s0b_full$logloss_plot)

# Combine outliers from both scenarios
outlier_plots = outliers_s0a | outliers_s0b

cat(sprintf("Identified %d outlier plots out of %d (%.1f%%)\n",
            sum(outlier_plots),
            length(outlier_plots),
            100 * sum(outlier_plots) / length(outlier_plots)))

# Filter out outliers and recompute metrics
Y_train_clean = Y[!outlier_plots, model_sjsdm$species]
pred_s0a_clean = predictions$s0a_train_marginal[!outlier_plots, ]
pred_s0b_clean = predictions$s0b_train_conditional[!outlier_plots, ]

metrics$s0a = compute_metrics(Y_train_clean, pred_s0a_clean, "S0a: training marginal")
metrics$s0b = compute_metrics(Y_train_clean, pred_s0b_clean, "S0b: training conditional")

# Transplant experiment metrics
metrics$s1 = compute_metrics(Y_high, predictions$s1_high_elevation_marginal, "S1: high-elevation marginal")
metrics$s2 = compute_metrics(Y_high, predictions$s2_high_elevation_conditional, "S2: high-elevation conditional")
metrics$s3 = compute_metrics(Y_low, predictions$s3_low_elevation_marginal, "S3: low-elevation marginal")
metrics$s4 = compute_metrics(Y_low, predictions$s4_low_elevation_conditional, "S4: low-elevation conditional")
metrics$s5 = compute_metrics(Y_warm, predictions$s5_warmed_marginal, "S5: warmed marginal")
metrics$s6 = compute_metrics(Y_warm, predictions$s6_warmed_residents, "S6: warmed residents conditional")
metrics$s7 = compute_metrics(Y_warm, predictions$s7_warmed_colonizers, "S7: warmed colonizers conditional")
metrics$s8 = compute_metrics(Y_warm, predictions$s8_warmed_all, "S8: warmed conditional all")

cat("Metrics calculated for all scenarios\n")

# Function to summarize metrics with bootstrap CI
summarize_metric = function(metric_values, n_boot = 1000) {
  boot_fun = function(data, indices) {
    median(data[indices], na.rm = TRUE)
  }

  boot_results = boot(metric_values, boot_fun, R = n_boot)
  ci = boot.ci(boot_results, type = "perc")

  return(list(
    median = median(metric_values, na.rm = TRUE),
    mean = mean(metric_values, na.rm = TRUE),
    ci_lower = ci$percent[4],
    ci_upper = ci$percent[5]
  ))
}

# Summarize metrics across scenarios
summary_table = map_dfr(metrics, function(m) {

  logloss_plot_summary = summarize_metric(m$logloss_plot)
  bc_summary = summarize_metric(m$bc_plot)

  tibble(
    scenario = m$scenario,
    logloss_median = logloss_plot_summary$median,
    logloss_ci_lower = logloss_plot_summary$ci_lower,
    logloss_ci_upper = logloss_plot_summary$ci_upper,
    bc_median = bc_summary$median,
    bc_ci_lower = bc_summary$ci_lower,
    bc_ci_upper = bc_summary$ci_upper
  )
})

cat("\nMetric Summary:\n")
print(summary_table)

# Save summary table
write.csv(summary_table, "output/jsdm_prediction_metrics_summary.csv", row.names = FALSE)

# ==============================================================================
# 9. COMPUTE CONTRASTS
# ==============================================================================
cat("\n=== Computing contrasts ===\n")

contrasts = list()

# Contrast A: Conditional Gain (Δ)
cat("\nContrast A: Conditional Gain\n")

# For each plot set, calculate ΔM = M_conditional - M_marginal
contrasts$A_training = metrics$s0b$logloss_plot - metrics$s0a$logloss_plot
contrasts$A_high_elevation = metrics$s2$logloss_plot - metrics$s1$logloss_plot
contrasts$A_low_elevation = metrics$s4$logloss_plot - metrics$s3$logloss_plot
contrasts$A_warmed_residents = metrics$s6$logloss_plot - metrics$s5$logloss_plot
contrasts$A_warmed_colonizers = metrics$s7$logloss_plot - metrics$s5$logloss_plot
contrasts$A_warmed_all = metrics$s8$logloss_plot - metrics$s5$logloss_plot

# Summarize
contrast_A_summary = tibble(
  contrast = c("training", "high_elevation", "low_elevation", "warmed_residents", "warmed_colonizers", "warmed_all"),
  delta_median = c(
    median(contrasts$A_training),
    median(contrasts$A_high_elevation),
    median(contrasts$A_low_elevation),
    median(contrasts$A_warmed_residents),
    median(contrasts$A_warmed_colonizers),
    median(contrasts$A_warmed_all)
  )
)

cat("Conditional Gain (negative = improvement):\n")
print(contrast_A_summary)

# Create long-format dataframe for contrasts (similar to metrics_long)
cat("\n=== Creating contrasts long dataframe ===\n")

contrasts_long = bind_rows(
  # Contrast A: Conditional Gain (per-plot values)
  tibble(
    contrast_type = "A: Conditional Gain",
    scenario = "Training",
    plot_id = 1:length(contrasts$A_training),
    value = contrasts$A_training,
    treatment = "Training"
  ),
  tibble(
    contrast_type = "A: Conditional Gain",
    scenario = "High elevation",
    plot_id = 1:length(contrasts$A_high_elevation),
    value = contrasts$A_high_elevation,
    treatment = "High-elevation"
  ),
  tibble(
    contrast_type = "A: Conditional Gain",
    scenario = "Low elevation",
    plot_id = 1:length(contrasts$A_low_elevation),
    value = contrasts$A_low_elevation,
    treatment = "Low-elevation"
  ),
  tibble(
    contrast_type = "A: Conditional Gain",
    scenario = "Warmed - Residents",
    plot_id = 1:length(contrasts$A_warmed_residents),
    value = contrasts$A_warmed_residents,
    treatment = "Warmed"
  ),
  tibble(
    contrast_type = "A: Conditional Gain",
    scenario = "Warmed - Colonizers",
    plot_id = 1:length(contrasts$A_warmed_colonizers),
    value = contrasts$A_warmed_colonizers,
    treatment = "Warmed"
  ),
  tibble(
    contrast_type = "A: Conditional Gain",
    scenario = "Warmed - All",
    plot_id = 1:length(contrasts$A_warmed_all),
    value = contrasts$A_warmed_all,
    treatment = "Warmed"
  )
)

cat(sprintf("Created contrasts_long with %d rows\n", nrow(contrasts_long)))

# ==============================================================================
# 10. VISUALIZATIONS
# ==============================================================================
cat("\n=== Creating visualizations ===\n")

# Figure 1: Scenario Comparison - Boxplots
cat("Creating Figure 1: Scenario Comparison\n")

metrics_long = map_dfr(names(metrics), function(s) {
  tibble(
    scenario = metrics[[s]]$scenario,
    plot_id = 1:length(metrics[[s]]$logloss_plot),
    logloss = metrics[[s]]$logloss_plot,
    bray_curtis = metrics[[s]]$bc_plot,
    treatment = case_when(
      grepl("Training", scenario) ~ "Training",
      grepl("high-elevation", scenario) ~ "High-elevation",
      grepl("low-elevation", scenario) ~ "Low-elevation",
      TRUE ~ "Warmed"
    ),
    prediction_type = case_when(
      grepl("Marginal", scenario) ~ "Marginal",
      grepl("Conditional|Residents|Colonizers|All", scenario) ~ "Conditional"
    )
  )
})


# Create combined dataframe with marginal, conditional, and delta
combined_data = bind_rows(
  # Training
  tibble(
    treatment = "Training",
    scenario = "Training",
    metric_type = "Marginal",
    plot_id = 1:length(metrics$s0a$logloss_plot),
    value = metrics$s0a$logloss_plot
  ),
  tibble(
    treatment = "Training",
    scenario = "Training",
    metric_type = "Conditional",
    plot_id = 1:length(metrics$s0b$logloss_plot),
    value = metrics$s0b$logloss_plot
  ),
  tibble(
    treatment = "Training",
    scenario = "Training",
    metric_type = "ΔM (Gain)",
    plot_id = 1:length(contrasts$A_training),
    value = contrasts$A_training
  ),
  # High elevation
  tibble(
    treatment = "High-elevation",
    scenario = "High elevation",
    metric_type = "Marginal",
    plot_id = 1:length(metrics$s1$logloss_plot),
    value = metrics$s1$logloss_plot
  ),
  tibble(
    treatment = "High-elevation",
    scenario = "High elevation",
    metric_type = "Conditional",
    plot_id = 1:length(metrics$s2$logloss_plot),
    value = metrics$s2$logloss_plot
  ),
  tibble(
    treatment = "High-elevation",
    scenario = "High elevation",
    metric_type = "ΔM (Gain)",
    plot_id = 1:length(contrasts$A_high_elevation),
    value = contrasts$A_high_elevation
  ),
  # Low elevation
  tibble(
    treatment = "Low-elevation",
    scenario = "Low elevation",
    metric_type = "Marginal",
    plot_id = 1:length(metrics$s3$logloss_plot),
    value = metrics$s3$logloss_plot
  ),
  tibble(
    treatment = "Low-elevation",
    scenario = "Low elevation",
    metric_type = "Conditional",
    plot_id = 1:length(metrics$s4$logloss_plot),
    value = metrics$s4$logloss_plot
  ),
  tibble(
    treatment = "Low-elevation",
    scenario = "Low elevation",
    metric_type = "ΔM (Gain)",
    plot_id = 1:length(contrasts$A_low_elevation),
    value = contrasts$A_low_elevation
  ),
  # Warmed - Residents
  tibble(
    treatment = "Warmed",
    scenario = "Residents",
    metric_type = "Marginal",
    plot_id = 1:length(metrics$s5$logloss_plot),
    value = metrics$s5$logloss_plot
  ),
  tibble(
    treatment = "Warmed",
    scenario = "Residents",
    metric_type = "Conditional",
    plot_id = 1:length(metrics$s6$logloss_plot),
    value = metrics$s6$logloss_plot
  ),
  tibble(
    treatment = "Warmed",
    scenario = "Residents",
    metric_type = "ΔM (Gain)",
    plot_id = 1:length(contrasts$A_warmed_residents),
    value = contrasts$A_warmed_residents
  ),
  # Warmed - Colonizers
  tibble(
    treatment = "Warmed",
    scenario = "Colonizers",
    metric_type = "Marginal",
    plot_id = 1:length(metrics$s5$logloss_plot),
    value = metrics$s5$logloss_plot
  ),
  tibble(
    treatment = "Warmed",
    scenario = "Colonizers",
    metric_type = "Conditional",
    plot_id = 1:length(metrics$s7$logloss_plot),
    value = metrics$s7$logloss_plot
  ),
  tibble(
    treatment = "Warmed",
    scenario = "Colonizers",
    metric_type = "ΔM (Gain)",
    plot_id = 1:length(contrasts$A_warmed_colonizers),
    value = contrasts$A_warmed_colonizers
  ),
  # Warmed - All
  tibble(
    treatment = "Warmed",
    scenario = "All",
    metric_type = "Marginal",
    plot_id = 1:length(metrics$s5$logloss_plot),
    value = metrics$s5$logloss_plot
  ),
  tibble(
    treatment = "Warmed",
    scenario = "All",
    metric_type = "Conditional",
    plot_id = 1:length(metrics$s8$logloss_plot),
    value = metrics$s8$logloss_plot
  ),
  tibble(
    treatment = "Warmed",
    scenario = "All",
    metric_type = "ΔM (Gain)",
    plot_id = 1:length(contrasts$A_warmed_all),
    value = contrasts$A_warmed_all
  )
) %>%
  mutate(
    metric_type = factor(metric_type, levels = c("Marginal", "Conditional", "ΔM (Gain)")),
    treatment = factor(treatment, levels = c("Training", "High-elevation", "Low-elevation", "Warmed"))
  )

# Calculate medians for overlay
median_data = combined_data %>%
  group_by(treatment, scenario, metric_type) %>%
  summarize(
    median_value = median(value, na.rm = TRUE),
    .groups = "drop"
  )

p1 = ggplot(combined_data, aes(x = scenario, y = value, fill = metric_type)) +
  geom_violin(position = position_dodge(width = 0.8), alpha = 0.7) +
  # Add median points
  geom_point(
    data = median_data,
    aes(x = scenario, y = median_value, group = metric_type),
    position = position_dodge(width = 0.8),
    size = 3,
    shape = 18,
    color = "black"
  ) +
  # Add median horizontal lines
  geom_crossbar(
    data = median_data,
    aes(x = scenario, y = median_value, ymin = median_value, ymax = median_value, group = metric_type),
    position = position_dodge(width = 0.8),
    width = 0.15,
    color = "black",
    linewidth = 0.5,
    fatten = 0
  ) +
  # Add median text labels
  geom_text(
    data = median_data,
    aes(x = scenario, y = median_value, label = round(median_value, 2), group = metric_type),
    position = position_dodge(width = 0.8),
    vjust = -0.8,
    hjust = -0.3,
    size = 3,
    fontface = "bold"
  ) +
  facet_wrap(~treatment, scales = "free_x", ncol = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  scale_fill_manual(
    values = c(
      "Marginal" = "#E69F00",
      "Conditional" = "#56B4E9",
      "ΔM (Gain)" = "#009E73"
    ),
    name = NULL
  ) +
  labs(
    title = "Log-Loss and Conditional Gain by Treatment",
    subtitle = "ΔM = Conditional - Marginal (negative ΔM = improvement) | Numbers show medians",
    x = "Scenario",
    y = "Value"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

print(p1)

# Figure 1c: Combined Bray-Curtis and Differences
cat("Creating Figure 1c: Combined Bray-Curtis Dissimilarity\n")

# Create combined dataframe for Bray-Curtis with marginal, conditional, and delta
combined_data_bc = bind_rows(
  # Training
  tibble(
    treatment = "Training",
    scenario = "Training",
    metric_type = "Marginal",
    plot_id = 1:length(metrics$s0a$bc_plot),
    value = metrics$s0a$bc_plot
  ),
  tibble(
    treatment = "Training",
    scenario = "Training",
    metric_type = "Conditional",
    plot_id = 1:length(metrics$s0b$bc_plot),
    value = metrics$s0b$bc_plot
  ),
  tibble(
    treatment = "Training",
    scenario = "Training",
    metric_type = "ΔBC (Gain)",
    plot_id = 1:length(metrics$s0b$bc_plot),
    value = metrics$s0b$bc_plot - metrics$s0a$bc_plot
  ),
  # High elevation
  tibble(
    treatment = "High-elevation",
    scenario = "High elevation",
    metric_type = "Marginal",
    plot_id = 1:length(metrics$s1$bc_plot),
    value = metrics$s1$bc_plot
  ),
  tibble(
    treatment = "High-elevation",
    scenario = "High elevation",
    metric_type = "Conditional",
    plot_id = 1:length(metrics$s2$bc_plot),
    value = metrics$s2$bc_plot
  ),
  tibble(
    treatment = "High-elevation",
    scenario = "High elevation",
    metric_type = "ΔBC (Gain)",
    plot_id = 1:length(metrics$s2$bc_plot),
    value = metrics$s2$bc_plot - metrics$s1$bc_plot
  ),
  # Low elevation
  tibble(
    treatment = "Low-elevation",
    scenario = "Low elevation",
    metric_type = "Marginal",
    plot_id = 1:length(metrics$s3$bc_plot),
    value = metrics$s3$bc_plot
  ),
  tibble(
    treatment = "Low-elevation",
    scenario = "Low elevation",
    metric_type = "Conditional",
    plot_id = 1:length(metrics$s4$bc_plot),
    value = metrics$s4$bc_plot
  ),
  tibble(
    treatment = "Low-elevation",
    scenario = "Low elevation",
    metric_type = "ΔBC (Gain)",
    plot_id = 1:length(metrics$s4$bc_plot),
    value = metrics$s4$bc_plot - metrics$s3$bc_plot
  ),
  # Warmed - Residents
  tibble(
    treatment = "Warmed",
    scenario = "Residents",
    metric_type = "Marginal",
    plot_id = 1:length(metrics$s5$bc_plot),
    value = metrics$s5$bc_plot
  ),
  tibble(
    treatment = "Warmed",
    scenario = "Residents",
    metric_type = "Conditional",
    plot_id = 1:length(metrics$s6$bc_plot),
    value = metrics$s6$bc_plot
  ),
  tibble(
    treatment = "Warmed",
    scenario = "Residents",
    metric_type = "ΔBC (Gain)",
    plot_id = 1:length(metrics$s6$bc_plot),
    value = metrics$s6$bc_plot - metrics$s5$bc_plot
  ),
  # Warmed - Colonizers
  tibble(
    treatment = "Warmed",
    scenario = "Colonizers",
    metric_type = "Marginal",
    plot_id = 1:length(metrics$s5$bc_plot),
    value = metrics$s5$bc_plot
  ),
  tibble(
    treatment = "Warmed",
    scenario = "Colonizers",
    metric_type = "Conditional",
    plot_id = 1:length(metrics$s7$bc_plot),
    value = metrics$s7$bc_plot
  ),
  tibble(
    treatment = "Warmed",
    scenario = "Colonizers",
    metric_type = "ΔBC (Gain)",
    plot_id = 1:length(metrics$s7$bc_plot),
    value = metrics$s7$bc_plot - metrics$s5$bc_plot
  ),
  # Warmed - All
  tibble(
    treatment = "Warmed",
    scenario = "All",
    metric_type = "Marginal",
    plot_id = 1:length(metrics$s5$bc_plot),
    value = metrics$s5$bc_plot
  ),
  tibble(
    treatment = "Warmed",
    scenario = "All",
    metric_type = "Conditional",
    plot_id = 1:length(metrics$s8$bc_plot),
    value = metrics$s8$bc_plot
  ),
  tibble(
    treatment = "Warmed",
    scenario = "All",
    metric_type = "ΔBC (Gain)",
    plot_id = 1:length(metrics$s8$bc_plot),
    value = metrics$s8$bc_plot - metrics$s5$bc_plot
  )
) %>%
  mutate(
    metric_type = factor(metric_type, levels = c("Marginal", "Conditional", "ΔBC (Gain)")),
    treatment = factor(treatment, levels = c("Training", "High-elevation", "Low-elevation", "Warmed"))
  )

# Calculate medians for Bray-Curtis
median_data_bc = combined_data_bc %>%
  group_by(treatment, scenario, metric_type) %>%
  summarize(
    median_value = median(value, na.rm = TRUE),
    .groups = "drop"
  )

p1c = ggplot(combined_data_bc, aes(x = scenario, y = value, fill = metric_type)) +
  geom_violin(position = position_dodge(width = 0.8), alpha = 0.7) +
  # Add median points
  geom_point(
    data = median_data_bc,
    aes(x = scenario, y = median_value, group = metric_type),
    position = position_dodge(width = 0.8),
    size = 3,
    shape = 18,
    color = "black"
  ) +
  # Add median horizontal lines
  geom_crossbar(
    data = median_data_bc,
    aes(x = scenario, y = median_value, ymin = median_value, ymax = median_value, group = metric_type),
    position = position_dodge(width = 0.8),
    width = 0.15,
    color = "black",
    linewidth = 0.5,
    fatten = 0
  ) +
  # Add median text labels
  geom_text(
    data = median_data_bc,
    aes(x = scenario, y = median_value, label = round(median_value, 2), group = metric_type),
    position = position_dodge(width = 0.8),
    vjust = -0.8,
    hjust = -0.3,
    size = 3,
    fontface = "bold"
  ) +
  facet_wrap(~treatment, scales = "free_x", ncol = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  scale_fill_manual(
    values = c(
      "Marginal" = "#E69F00",
      "Conditional" = "#56B4E9",
      "ΔBC (Gain)" = "#009E73"
    ),
    name = NULL
  ) +
  labs(
    title = "Bray-Curtis Dissimilarity and Conditional Gain by Treatment",
    subtitle = "ΔBC = Conditional - Marginal (negative ΔBC = improvement) | Numbers show medians",
    x = "Scenario",
    y = "Bray-Curtis Dissimilarity"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

print(p1c)

# Get species-level log-loss for each scenario
species_logloss = map_dfr(names(metrics), function(s) {
  tibble(
    scenario = metrics[[s]]$scenario,
    species = names(metrics[[s]]$logloss_species),
    logloss = metrics[[s]]$logloss_species
  )
}) %>%
  left_join(species_pools, by = "species")

p6 = ggplot(species_logloss, aes(x = scenario, y = reorder(species, logloss), fill = logloss)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = median(species_logloss$logloss)) +
  facet_wrap(~pool, scales = "free_y", ncol = 1) +
  labs(
    title = "Species-Level Log-Loss Across Scenarios",
    x = "Scenario",
    y = "Species",
    fill = "Log-Loss"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# Save all figures to PDF
cat("Saving figures to PDF\n")

pdf("plot/predictions/jsdm_prediction_analysis.pdf", width = 14, height = 10)
print(p1)
print(p1b)
print(p1c)
print(p2)
print(p2b)
print(p3)
print(p4)
print(p5)
print(p6)
dev.off()

cat("\n=== Analysis complete! ===\n")
cat("\nOutputs saved:\n")
cat("  - output/jsdm_prediction_metrics_summary.csv\n")
cat("  - output/jsdm_prediction_contrasts.rds\n")
cat("  - plot/predictions/jsdm_prediction_analysis.pdf\n")

### Assess the colonizers vs. residents ----
cat("\n=== Extracting species associations ===\n")

# Extract correlation matrix from the model
# sjSDM stores species associations as correlations in the covariance matrix
association_matrix = getCor(model_sjsdm)
rownames(association_matrix) = colnames(association_matrix) = model_sjsdm$species

cat(sprintf("Association matrix: %d × %d species\n", nrow(association_matrix), ncol(association_matrix)))

# Subset to species present in warmed plots
warmed_species = union(warmed_colonizers, warmed_residents)
association_warmed = association_matrix[warmed_species, warmed_species]

# Further subset to show colonizers vs residents
# Rows = colonizers, Columns = residents (or vice versa)
association_colonizers_residents = association_matrix[warmed_colonizers, warmed_residents]

cat(sprintf("Colonizers vs Residents matrix: %d colonizers × %d residents\n",
            nrow(association_colonizers_residents),
            ncol(association_colonizers_residents)))

# Convert to long format for ggplot
association_long = association_colonizers_residents %>%
  as.data.frame() %>%
  rownames_to_column("colonizer") %>%
  pivot_longer(cols = -colonizer, names_to = "resident", values_to = "correlation")

# Create heatmap
cat("\n=== Creating association heatmap ===\n")

p_assoc = ggplot(association_long, aes(x = resident, y = colonizer, fill = correlation)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(
    low = "#2166AC",    # blue for negative
    mid = "white",      # white for zero
    high = "#B2182B",   # red for positive
    midpoint = 0,
    limits = c(-1, 1),
    name = "Association"
  ) +
  labs(
    title = "Species Associations: Warmed Colonizers vs. Warmed Residents",
    subtitle = "Biotic associations from JSDM residual covariance matrix",
    x = "Resident Species (from high elevation)",
    y = "Colonizer Species (from low elevation)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right"
  )

print(p_assoc)

# Save heatmap
ggsave("plot/predictions/association_heatmap_colonizers_residents.pdf",
       plot = p_assoc,
       width = 12,
       height = 10)

ggsave("plot/predictions/association_heatmap_colonizers_residents.png",
       plot = p_assoc,
       width = 12,
       height = 10,
       dpi = 300)

# Summary statistics
cat("\n=== Association Summary Statistics ===\n")
cat(sprintf("Mean association: %.3f\n", mean(association_colonizers_residents)))
cat(sprintf("Median association: %.3f\n", median(association_colonizers_residents)))
cat(sprintf("SD association: %.3f\n", sd(association_colonizers_residents)))
cat(sprintf("Range: [%.3f, %.3f]\n",
            min(association_colonizers_residents),
            max(association_colonizers_residents)))

# Count positive and negative associations
n_positive = sum(association_colonizers_residents > 0)
n_negative = sum(association_colonizers_residents < 0)
n_total = length(association_colonizers_residents)

cat(sprintf("\nPositive associations: %d (%.1f%%)\n",
            n_positive, 100 * n_positive / n_total))
cat(sprintf("Negative associations: %d (%.1f%%)\n",
            n_negative, 100 * n_negative / n_total))

# Find strongest associations
association_summary = association_long %>%
  mutate(abs_correlation = abs(correlation)) %>%
  arrange(desc(abs_correlation))

cat("\nTop 10 strongest associations:\n")
print(association_summary %>% head(10))

# Save association matrix
write.csv(association_colonizers_residents,
          "output/association_matrix_colonizers_residents.csv")

cat("\nAssociation analysis complete!\n")
cat("Outputs saved:\n")
cat("  - output/association_matrix_colonizers_residents.csv\n")
cat("  - plot/predictions/association_heatmap_colonizers_residents.pdf\n")
cat("  - plot/predictions/association_heatmap_colonizers_residents.png\n")

# ==============================================================================
# MISSING SPECIES ANALYSIS
# ==============================================================================
cat("\n=== Analyzing species absent from JSDM model ===\n")

# Get all species from transplant data (before filtering to common species)
transplant_species_all = colnames(Y_all)

# Identify species present in transplants but absent from model
missing_species = setdiff(transplant_species_all, model_sjsdm$species)

cat(sprintf("Total species in transplant data: %d\n", length(transplant_species_all)))
cat(sprintf("Species in JSDM model: %d\n", length(model_sjsdm$species)))
cat(sprintf("Missing species (in transplants, not in model): %d\n", length(missing_species)))

if(length(missing_species) > 0) {
  cat("\nMissing species:\n")
  print(missing_species)

  # Calculate cover contribution of missing species in each treatment group
  # Sum relative cover of missing species per plot (community)

  missing_relcover_per_plot = calanda_data %>%
    group_by(plot_id, plot_type) %>%
    summarize(
      total_relcover_all = sum(Rel_Cover, na.rm = TRUE),
      missing_relcover = sum(Rel_Cover[SpeciesName %in% missing_species], na.rm = TRUE),
      modeled_relcover = sum(Rel_Cover[!SpeciesName %in% missing_species], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      percent_missing = 100 * missing_relcover,
      percent_modeled = 100 * modeled_relcover
    )

  # Summary statistics by treatment
  cover_summary = missing_relcover_per_plot %>%
    group_by(plot_type) %>%
    summarize(
      n_plots = n(),
      mean_percent_missing = mean(percent_missing, na.rm = TRUE),
      median_percent_missing = median(percent_missing, na.rm = TRUE),
      sd_percent_missing = sd(percent_missing, na.rm = TRUE),
      min_percent_missing = min(percent_missing, na.rm = TRUE),
      max_percent_missing = max(percent_missing, na.rm = TRUE),
      .groups = "drop"
    )

  cat("\nRelative cover contribution by treatment:\n")
  print(cover_summary)

  # Count missing species occurrences by treatment
  missing_species_counts = calanda_data %>%
    filter(SpeciesName %in% missing_species) %>%
    group_by(plot_type, SpeciesName) %>%
    summarize(
      n_occurrences = n(),
      total_cover = sum(Cover, na.rm = TRUE),
      mean_cover = mean(Cover, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(plot_type, desc(total_cover))

  cat("\nTop missing species by cover in each treatment:\n")
  missing_species_counts %>%
    group_by(plot_type) %>%
    slice_head(n = 5) %>%
    print()

  # Create visualization 1: Violin plot of summed relative cover of missing species per plot
  cat("\n=== Creating missing species visualizations ===\n")

  # Prepare data for plotting
  plot_data = missing_relcover_per_plot %>%
    mutate(
      plot_type = factor(plot_type,
                        levels = c("high_elevation", "low_elevation", "warmed"),
                        labels = c("High elevation", "Low elevation", "Warmed"))
    )

  # Calculate summary statistics for annotation
  summary_stats = cover_summary %>%
    mutate(
      plot_type = factor(plot_type,
                        levels = c("high_elevation", "low_elevation", "warmed"),
                        labels = c("High elevation", "Low elevation", "Warmed"))
    )

  p_cover = ggplot(plot_data, aes(x = plot_type, y = percent_missing, fill = plot_type)) +
    geom_violin(alpha = 0.6, trim = FALSE) +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 1.5) +
    # Add median text labels
    geom_text(
      data = summary_stats,
      aes(x = plot_type, y = max_percent_missing + 2,
          label = sprintf("Median: %.1f%%\nMean: %.1f%%",
                         median_percent_missing, mean_percent_missing)),
      size = 3.5,
      fontface = "bold",
      inherit.aes = FALSE
    ) +
    scale_fill_manual(
      values = c("High elevation" = "#1B9E77",
                 "Low elevation" = "#D95F02",
                 "Warmed" = "#7570B3")
    ) +
    labs(
      title = "Relative Cover of Species Absent from JSDM Model",
      subtitle = sprintf("%d species present in transplants but absent from model (summed per plot)",
                        length(missing_species)),
      x = "Treatment Group",
      y = "Summed Relative Cover of Missing Species (%)"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 12)
    )

  print(p_cover)

  # Create visualization 2: Missing species richness and abundance
  missing_richness = calanda_data %>%
    filter(SpeciesName %in% missing_species) %>%
    group_by(plot_type, plot_id) %>%
    summarize(
      n_missing_species = n_distinct(SpeciesName),
      total_missing_cover = sum(Cover, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      plot_type = factor(plot_type, levels = c("high_elevation", "low_elevation", "warmed"),
                        labels = c("High elevation", "Low elevation", "Warmed"))
    )

  p_richness = ggplot(missing_richness, aes(x = plot_type, y = n_missing_species, fill = plot_type)) +
    geom_violin(alpha = 0.5, show.legend = FALSE) +
    geom_boxplot(width = 0.2, alpha = 0.8, show.legend = FALSE) +
    scale_fill_manual(values = c("High elevation" = "#1B9E77",
                                   "Low elevation" = "#D95F02",
                                   "Warmed" = "#7570B3")) +
    labs(
      title = "Missing Species Richness per Plot",
      subtitle = "Number of species absent from JSDM model",
      x = "Treatment Group",
      y = "Number of Missing Species per Plot"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 12)
    )

  print(p_richness)

  # Create visualization 3: Top missing species by treatment
  top_missing = missing_species_counts %>%
    group_by(plot_type) %>%
    slice_max(order_by = total_cover, n = 10) %>%
    ungroup() %>%
    mutate(
      plot_type = factor(plot_type, levels = c("high_elevation", "low_elevation", "warmed"),
                        labels = c("High elevation", "Low elevation", "Warmed"))
    )

  p_top_missing = ggplot(top_missing, aes(x = reorder(SpeciesName, total_cover),
                                           y = total_cover,
                                           fill = plot_type)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~plot_type, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c("High elevation" = "#1B9E77",
                                   "Low elevation" = "#D95F02",
                                   "Warmed" = "#7570B3")) +
    labs(
      title = "Top 10 Missing Species by Total Cover",
      subtitle = "Species present in transplants but absent from JSDM model",
      x = "Species",
      y = "Total Cover"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 12),
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 11)
    )

  print(p_top_missing)

  # Combined figure
  p_missing_combined = (p_cover / p_richness) +
    plot_annotation(
      title = "Species Absent from JSDM Model",
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )

  print(p_missing_combined)

  # Save figures
  ggsave("plot/predictions/missing_species_cover_contribution.pdf",
         plot = p_cover,
         width = 10,
         height = 7)

  ggsave("plot/predictions/missing_species_richness.pdf",
         plot = p_richness,
         width = 10,
         height = 7)

  ggsave("plot/predictions/missing_species_top10.pdf",
         plot = p_top_missing,
         width = 10,
         height = 12)

  ggsave("plot/predictions/missing_species_combined.pdf",
         plot = p_missing_combined,
         width = 12,
         height = 10)

  # Save summary tables
  write.csv(cover_summary, "output/missing_species_cover_summary.csv", row.names = FALSE)
  write.csv(missing_relcover_per_plot, "output/missing_species_relcover_per_plot.csv", row.names = FALSE)
  write.csv(missing_species_counts, "output/missing_species_by_treatment.csv", row.names = FALSE)
  write.csv(data.frame(missing_species = missing_species),
            "output/missing_species_list.csv", row.names = FALSE)

  cat("\nMissing species analysis complete!\n")
  cat("Outputs saved:\n")
  cat("  - output/missing_species_cover_summary.csv\n")
  cat("  - output/missing_species_relcover_per_plot.csv\n")
  cat("  - output/missing_species_by_treatment.csv\n")
  cat("  - output/missing_species_list.csv\n")
  cat("  - plot/predictions/missing_species_cover_contribution.pdf\n")
  cat("  - plot/predictions/missing_species_richness.pdf\n")
  cat("  - plot/predictions/missing_species_top10.pdf\n")
  cat("  - plot/predictions/missing_species_combined.pdf\n")

} else {
  cat("\nAll transplant species are included in the JSDM model.\n")
}

cat("\n=== SCRIPT COMPLETE ===\n")



