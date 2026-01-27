
# Load required libraries ----
setwd(here("Calanda_JSDM"))

# Load data ----
load("output/starter_data_25.04.25.RData")
an = readRDS("results_from_Max/an_sjsdm_calanda.rds")
model_jsdm = readRDS("results_from_Max/model_sjsdm_calanda.rds")
environment(model_jsdm$get_model)$device = "cpu"
res = readRDS("results_from_Max/res_sjsdm_calanda.rds")

# ==============================================================================
# VARIANCE PARTITIONING: VENN DIAGRAM + TERNARY PLOTS ----
# ==============================================================================
cat("\n=== Creating variance partitioning visualizations ===\n")

# Define custom colors
color_codist = "#00bd89"  # Green for species associations
color_spa = "#d00000"      # Red for space
color_env = "#81caf3"      # Blue for environment

# Get Venn
pdf(file = "plot/venn.pdf", height = 8, width = 8)
p_venn = plot.anova.custom(an)
dev.off()

# Calculate species altitude ranges
cat("Calculating species altitude ranges...\n")
species_altitude_ranges = veg.clim %>%
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

# Create combined plot function using custom ternary functions from functions_calanda.R
p1= plot_tern_sites(
    res = res,
    veg.clim = veg.clim,
    color_env = color_env,
    color_codist = color_codist,
    color_spa = color_spa
  )
  
# 2. Species ternary (colored by altitude range)
p2=plot_tern_species(
    res = res,
    veg = veg,
    veg.clim = veg.clim,
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

# ==============================================================================
# SPECIES-LEVEL ENVIRONMENTAL EFFECT SIZES  ----
# ==============================================================================
cat("\n=== Extracting species-level beta coefficients ===\n")

# Extract environmental coefficients and standard errors from sjSDM model
model_summary = summary(model_jsdm)
coef_env = model_summary$coefs[-1, ]  # Remove intercept
p_values = model_summary$P[-1, ]       # Standard errors

# Check structure
cat(sprintf("Extracted coefficients for %d species and %d environmental variables\n",
            ncol(coef_env), nrow(coef_env)))

rownames(p_values) = rownames(coef_env)
colnames(p_values) = colnames(coef_env)

# Create data frame with species as rows and environmental variables as columns
species_betas = as.data.frame(t(coef_env)) %>%
  rownames_to_column("species")

species_pvalues = as.data.frame(t(p_values)) %>%
  rownames_to_column("species")

# Get species variance components
species_variance = res$internals$Species %>%
  rownames_to_column("species")

# Prepare data for heatmap
heatmap_data = species_betas %>%
  pivot_longer(cols = -species, names_to = "variable", values_to = "beta") %>%
  left_join(
    species_pvalues %>%
      pivot_longer(cols = -species, names_to = "variable", values_to = "pvalue"),
    by = c("species", "variable")
  ) %>%
  left_join(species_variance %>% select(species, codist), by = "species") %>%
  mutate(
    significant = case_when(
      pvalue < 0.001 ~ "***",
      pvalue < 0.01 ~ "**",
      pvalue < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# Order species by total environmental variance explained
species_order = species_variance %>%
  arrange(desc(codist)) %>%
  pull(species)

heatmap_data$species = factor(heatmap_data$species, levels = species_order)

# Create double heatmap with beta coefficients and codist indicator
p_species_heatmap = ggplot(heatmap_data, aes(x = variable, y = species)) +
  # Background tiles colored by beta coefficient
  geom_tile(aes(fill = beta), color = "white", linewidth = 0.5) +
  scale_fill_gradient2(
    low = rgb(0, 0, 1),      # Blue for negative
    mid = "white",            # White for zero
    high = rgb(1, 0, 0),      # Red for positive
    midpoint = 0,
    name = "Beta\ncoefficient"
  ) +
  # Add points colored by codist R²
  geom_point(aes(color = sqrt(codist)), size = 2, alpha = 0.8) +
  scale_color_gradient(
    low = "white",
    high = rgb(0, 1, 0),     # Green for codist
    name = "Species assoc.\nR²"
  ) +
  # Add significance stars
  geom_text(aes(label = significant), size = 3, fontface = "bold",
            color = "black", vjust = 0.75) +
  # Formatting
  labs(
    x = "Environmental variable",
    y = "Species",
    title = "Species-specific environmental effect sizes",
    subtitle = "Tile color = beta coefficient | Point color = species associations R² | * p<0.05, ** p<0.01, *** p<0.001"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 6),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10),
    legend.position = "right",
    panel.grid = element_blank()
  ) +
  coord_flip()
print(p_species_heatmap)
cat("\nSpecies beta coefficient heatmap created\n")

# ==============================================================================
# REGRESSION: CODIST R² vs BETA COEFFICIENTS ----
# ==============================================================================
cat("\n=== Analyzing relationship between species associations and environmental effects ===\n")

# Prepare data: merge betas with codist variance
regression_data = species_betas %>%
  left_join(species_variance %>% select(species, codist, env, spa), by = "species") %>%
  drop_na()

cat(sprintf("Starting with %d species\n", nrow(regression_data)))

# Outlier detection on codist R² using IQR method
Q1 = quantile(regression_data$codist, 0.001)
Q3 = quantile(regression_data$codist, 0.999)
IQR = Q3 - Q1
lower_bound = Q1 - 1.5 * IQR
upper_bound = Q3 + 1.5 * IQR

# Identify outliers
outliers = regression_data %>%
  filter(codist < lower_bound | codist > upper_bound)

cat(sprintf("Detected %d outliers in codist R²\n", nrow(outliers)))
if(nrow(outliers) > 0) {
  cat("Outlier species:\n")
  print(outliers %>% select(species, codist))
}

# Remove outliers
regression_data_clean = regression_data %>%
  filter(codist >= lower_bound & codist <= upper_bound)

cat(sprintf("After outlier removal: %d species retained\n", nrow(regression_data_clean)))

# Fit multiple regression with second-degree polynomials
# Select only specific predictors
predictor_names = c("summer_temp", "fdd", "et.annual", "land_use")

# Check which predictors are available in the data
available_predictors = predictor_names[predictor_names %in% names(regression_data_clean)]

if(length(available_predictors) < length(predictor_names)) {
  missing_preds = setdiff(predictor_names, available_predictors)
  cat("\nWarning: The following predictors are not in the data:\n")
  print(missing_preds)
}

cat("\nPredictors included:\n")
print(available_predictors)

# Build formula with linear and quadratic terms
linear_terms = paste(available_predictors, collapse = " + ")
quadratic_terms = paste("I(", available_predictors, "^2)", sep = "", collapse = " + ")
regression_formula = as.formula(paste("codist ~", linear_terms, "+", quadratic_terms))

cat("\nFitting full model with second-degree polynomials:\n")
print(regression_formula)

# Fit the full model
multi_model_full = lm(regression_formula, data = regression_data_clean)

cat("\nFull model summary:\n")
cat(sprintf("Full model AIC: %.2f\n", AIC(multi_model_full)))
cat(sprintf("Full model R²: %.3f\n", summary(multi_model_full)$r.squared))
cat(sprintf("Full model Adjusted R²: %.3f\n", summary(multi_model_full)$adj.r.squared))

# ==============================================================================
# PERFORM BACKWARD SELECTION USING step() AND CAPTURE OUTPUT ----
# ==============================================================================
cat("\n=== Performing backward selection based on AIC ===\n")

# Capture step() output to extract the process
step_output = capture.output({
  multi_model = step(multi_model_full, direction = "backward", trace = 2, k = 2)
})

cat("\n" , rep("=", 80), "\n", sep = "")
cat("STEPWISE SELECTION COMPLETE\n")
cat(rep("=", 80), "\n", sep = "")
cat(sprintf("Final model AIC: %.3f\n", AIC(multi_model)))
cat(sprintf("Final number of predictors: %d\n", length(names(coef(multi_model))[-1])))
cat(sprintf("Final formula: %s\n", as.character(formula(multi_model))[3]))

# ==============================================================================
# CREATE COMPREHENSIVE STEP-BY-STEP SUMMARY TABLE BY REPLICATING step() LOGIC
# ==============================================================================
cat("\n=== Creating step-by-step summary table ===\n")

# Now manually replicate the stepwise process to create detailed tables
current_model_manual = multi_model_full
step_history = list()
iteration = 0

repeat {
  iteration = iteration + 1

  # Get current model AIC and formula
  current_aic = AIC(current_model_manual)
  current_vars = names(coef(current_model_manual))[-1]  # Exclude intercept
  current_formula = as.character(formula(current_model_manual))[3]

  # Try removing each variable
  removal_results = data.frame(
    variable_removed = character(),
    new_aic = numeric(),
    aic_change = numeric(),
    stringsAsFactors = FALSE
  )

  for(var in current_vars) {
    # Create formula without this variable
    remaining_vars = setdiff(current_vars, var)

    if(length(remaining_vars) == 0) {
      # Can't remove last variable
      next
    }

    formula_reduced = as.formula(paste("codist ~", paste(remaining_vars, collapse = " + ")))

    # Fit reduced model
    model_reduced = lm(formula_reduced, data = regression_data_clean)
    new_aic = AIC(model_reduced)
    aic_change = new_aic - current_aic

    removal_results = rbind(
      removal_results,
      data.frame(
        variable_removed = var,
        new_aic = new_aic,
        aic_change = aic_change,
        stringsAsFactors = FALSE
      )
    )
  }

  # Sort by AIC (best = lowest)
  removal_results = removal_results %>%
    arrange(new_aic) %>%
    mutate(
      var_clean = gsub("I\\(", "", variable_removed),
      var_clean = gsub("\\^2\\)", "²", var_clean),
      decision = ifelse(row_number() == 1 & aic_change < 0, "REMOVE", "Keep")
    )

  # Store iteration history
  step_history[[iteration]] = list(
    iteration = iteration,
    current_formula = paste("codist ~", current_formula),
    current_aic = current_aic,
    n_predictors = length(current_vars),
    removal_options = removal_results,
    action = NA,
    variable_removed = NA
  )

  # Decision: remove variable with lowest AIC if it improves the model
  best_removal = removal_results[1, ]

  if(best_removal$aic_change < 0) {
    # Removal improves model - do it
    var_to_remove = best_removal$variable_removed
    remaining_vars = setdiff(current_vars, var_to_remove)
    new_formula = as.formula(paste("codist ~", paste(remaining_vars, collapse = " + ")))
    current_model_manual = lm(new_formula, data = regression_data_clean)

    step_history[[iteration]]$action = "REMOVED"
    step_history[[iteration]]$variable_removed = best_removal$var_clean
  } else {
    # No improvement possible - stop
    step_history[[iteration]]$action = "STOP"
    break
  }
}

# Table 1: Iteration Summary
iteration_summary = data.frame(
  step = integer(),
  model_formula = character(),
  n_predictors = integer(),
  aic = numeric(),
  action = character(),
  variable_removed = character(),
  aic_change = numeric(),
  stringsAsFactors = FALSE
)

for(i in 1:length(step_history)) {
  h = step_history[[i]]

  if(h$action == "REMOVED") {
    # Get AIC change from removal options
    removed_var_info = h$removal_options %>%
      filter(decision == "REMOVE")
    aic_change = removed_var_info$aic_change[1]
  } else {
    aic_change = NA
  }

  iteration_summary = rbind(
    iteration_summary,
    data.frame(
      step = i,
      model_formula = h$current_formula,
      n_predictors = h$n_predictors,
      aic = h$current_aic,
      action = h$action,
      variable_removed = ifelse(is.na(h$variable_removed), "-", h$variable_removed),
      aic_change = aic_change,
      stringsAsFactors = FALSE
    )
  )
}

cat("\n=== STEPWISE SELECTION SUMMARY ===\n")
print(iteration_summary, row.names = FALSE)

# Table 2: All variables tested at each step
all_tests = data.frame(
  step = integer(),
  variable_tested = character(),
  aic_if_removed = numeric(),
  aic_change = numeric(),
  decision = character(),
  stringsAsFactors = FALSE
)

for(i in 1:length(step_history)) {
  h = step_history[[i]]

  if(nrow(h$removal_options) > 0) {
    step_tests = h$removal_options %>%
      select(var_clean, new_aic, aic_change, decision) %>%
      mutate(step = i) %>%
      select(step, var_clean, new_aic, aic_change, decision)

    colnames(step_tests) = c("step", "variable_tested", "aic_if_removed", "aic_change", "decision")

    all_tests = rbind(all_tests, step_tests)
  }
}

cat("\n=== ALL VARIABLES TESTED AT EACH STEP ===\n")
print(all_tests, row.names = FALSE)

# Verify that manual process matches step() result
cat("\n=== VERIFICATION ===\n")
cat("step() function result:\n")
cat(sprintf("  Formula: %s\n", as.character(formula(multi_model))[3]))
cat(sprintf("  AIC: %.3f\n", AIC(multi_model)))
cat(sprintf("  N predictors: %d\n", length(names(coef(multi_model))[-1])))

cat("\nManual replication result:\n")
cat(sprintf("  Formula: %s\n", as.character(formula(current_model_manual))[3]))
cat(sprintf("  AIC: %.3f\n", AIC(current_model_manual)))
cat(sprintf("  N predictors: %d\n", length(names(coef(current_model_manual))[-1])))

cat("\nModels match:", identical(formula(multi_model), formula(current_model_manual)), "\n")

# Save tables
write.csv(iteration_summary,
          "output/species_stepwise_iteration_summary.csv",
          row.names = FALSE)
write.csv(all_tests,
          "output/species_stepwise_all_tests.csv",
          row.names = FALSE)

cat("\nTables saved:\n")
cat("  - output/species_stepwise_iteration_summary.csv\n")
cat("  - output/species_stepwise_all_tests.csv\n")

# Get summary of selected model
multi_summary = summary(multi_model)

cat("\n=== Selected model summary ===\n")
print(multi_summary)

cat(sprintf("\nSelected model AIC: %.2f\n", AIC(multi_model)))
cat(sprintf("Selected model R²: %.3f\n", multi_summary$r.squared))
cat(sprintf("Adjusted R²: %.3f\n", multi_summary$adj.r.squared))
cat(sprintf("F-statistic: %.3f (p = %.4f)\n",
            multi_summary$fstatistic[1],
            pf(multi_summary$fstatistic[1],
               multi_summary$fstatistic[2],
               multi_summary$fstatistic[3],
               lower.tail = FALSE)))

# Extract coefficients from selected model
regression_coefs = as.data.frame(coef(multi_summary)) %>%
  rownames_to_column("env_variable") %>%
  filter(env_variable != "(Intercept)") %>%
  rename(
    estimate = Estimate,
    std_error = `Std. Error`,
    t_value = `t value`,
    pvalue = `Pr(>|t|)`
  ) %>%
  mutate(
    lower_ci = estimate - 1.96 * std_error,
    upper_ci = estimate + 1.96 * std_error,
    significant = pvalue < 0.05,
    # Clean variable names for plotting
    var_type = ifelse(grepl("\\^2\\)", env_variable), "Quadratic", "Linear"),
    var_clean = gsub("I\\(", "", env_variable),
    var_clean = gsub("\\^2\\)", "²", var_clean)
  ) %>%
  arrange(desc(abs(estimate)))

cat("\n=== Selected regression coefficients ===\n")
print(regression_coefs %>% select(env_variable, var_type, estimate, pvalue))

cat(sprintf("\nTotal predictors in selected model: %d\n", nrow(regression_coefs)))
cat(sprintf("Linear terms: %d\n", sum(regression_coefs$var_type == "Linear")))
cat(sprintf("Quadratic terms: %d\n", sum(regression_coefs$var_type == "Quadratic")))
cat(sprintf("Significant predictors (p < 0.05): %d\n", sum(regression_coefs$significant)))

# Create effect size plot with polynomial terms distinguished
p_codist_beta_effects = ggplot(regression_coefs,
                               aes(x = reorder(var_clean, abs(estimate)),
                                   y = estimate,
                                   color = significant,
                                   shape = var_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),
                width = 0.3, linewidth = 1) +
  coord_flip() +
  scale_color_manual(values = c("TRUE" = rgb(0, 1, 0), "FALSE" = "gray70"),
                     labels = c("TRUE" = "p < 0.05", "FALSE" = "n.s."),
                     name = "Significance") +
  scale_shape_manual(values = c("Linear" = 16, "Quadratic" = 17),
                     name = "Term type") +
  labs(
    x = "Environmental variable (beta coefficient)",
    y = "Effect on species associations R²",
    title = "Multiple regression: Environmental effects → Species associations",
    subtitle = paste0("AIC-selected model with polynomials | Multiple R² = ",
                     round(multi_summary$r.squared, 3),
                     " | Adj. R² = ", round(multi_summary$adj.r.squared, 3),
                     " | n = ", nrow(regression_data_clean), " species")
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )
print(p_codist_beta_effects)

cat("\nMultiple regression analysis completed\n")
cat(sprintf("Significant predictors (p < 0.05): %d\n",
            sum(regression_coefs$significant)))

# Create scatter plots for significant predictors with predicted values
# Group by base variable (combine linear and quadratic terms)
sig_base_vars = regression_coefs %>%
  filter(significant == TRUE) %>%
  mutate(base_var = gsub("²", "", var_clean)) %>%
  pull(base_var) %>%
  unique()

if(length(sig_base_vars) > 0) {
  cat(sprintf("\nCreating scatter plots for %d significant base variables\n", length(sig_base_vars)))

  # For each significant base variable, create a scatter plot showing the polynomial fit
  partial_data_list = list()

  for(base_var in sig_base_vars) {
    # Check which terms are in the model for this variable
    var_terms = regression_coefs %>%
      filter(grepl(base_var, var_clean)) %>%
      select(env_variable, estimate, var_type)

    if(nrow(var_terms) > 0) {
      # Get beta values
      beta_values = regression_data_clean[[base_var]]

      # Calculate partial prediction for this variable
      # Start with intercept
      partial_pred = coef(multi_model)[1]

      # Add contribution from each term (linear and/or quadratic)
      for(i in 1:nrow(var_terms)) {
        term_name = var_terms$env_variable[i]
        coef_val = var_terms$estimate[i]

        if(var_terms$var_type[i] == "Linear") {
          partial_pred = partial_pred + coef_val * beta_values
        } else if(var_terms$var_type[i] == "Quadratic") {
          partial_pred = partial_pred + coef_val * beta_values^2
        }
      }

      # Center predictions around mean for visualization
      partial_pred = partial_pred + mean(regression_data_clean$codist) - mean(partial_pred)

      # Create sequence for smooth curve
      beta_seq = seq(min(beta_values), max(beta_values), length.out = 100)
      pred_seq = coef(multi_model)[1]

      for(i in 1:nrow(var_terms)) {
        term_name = var_terms$env_variable[i]
        coef_val = var_terms$estimate[i]

        if(var_terms$var_type[i] == "Linear") {
          pred_seq = pred_seq + coef_val * beta_seq
        } else if(var_terms$var_type[i] == "Quadratic") {
          pred_seq = pred_seq + coef_val * beta_seq^2
        }
      }

      pred_seq = pred_seq + mean(regression_data_clean$codist) - mean(pred_seq)

      # Determine plot label
      if(nrow(var_terms) == 2) {
        plot_label = paste0(base_var, " (linear + quadratic)")
      } else if(var_terms$var_type[1] == "Quadratic") {
        plot_label = paste0(base_var, " (quadratic only)")
      } else {
        plot_label = paste0(base_var, " (linear only)")
      }

      partial_data_list[[base_var]] = list(
        points = data.frame(
          beta = beta_values,
          codist = regression_data_clean$codist,
          env = regression_data_clean$env,
          env_variable = plot_label
        ),
        curve = data.frame(
          beta = beta_seq,
          predicted = pred_seq,
          env_variable = plot_label
        )
      )
    }
  }

  # Combine data for plotting
  points_data = do.call(rbind, lapply(partial_data_list, function(x) x$points))
  curve_data = do.call(rbind, lapply(partial_data_list, function(x) x$curve))

  # Update facet labels with descriptive names
  points_data = points_data %>%
    mutate(env_variable = case_when(
      grepl("et.annual", env_variable) ~ gsub("et.annual", "Annual evapotranspiration", env_variable),
      grepl("land_use", env_variable) ~ gsub("land_use", "Land use intensity", env_variable),
      grepl("summer_temp", env_variable) ~ gsub("summer_temp", "Summer temperature", env_variable),
      TRUE ~ env_variable
    ))

  curve_data = curve_data %>%
    mutate(env_variable = case_when(
      grepl("et.annual", env_variable) ~ gsub("et.annual", "Annual evapotranspiration", env_variable),
      grepl("land_use", env_variable) ~ gsub("land_use", "Land use intensity", env_variable),
      grepl("summer_temp", env_variable) ~ gsub("summer_temp", "Summer temperature", env_variable),
      TRUE ~ env_variable
    ))

  # Set factor order for facets: Summer temperature, Annual evapotranspiration, Land use intensity
  facet_order = c("Summer temperature (linear + quadratic)",
                  "Summer temperature (linear only)",
                  "Summer temperature (quadratic only)",
                  "Annual evapotranspiration (linear + quadratic)",
                  "Annual evapotranspiration (linear only)",
                  "Annual evapotranspiration (quadratic only)",
                  "Land use intensity (linear + quadratic)",
                  "Land use intensity (linear only)",
                  "Land use intensity (quadratic only)")

  points_data = points_data %>%
    mutate(env_variable = factor(env_variable,
                                 levels = facet_order[facet_order %in% unique(env_variable)]))

  curve_data = curve_data %>%
    mutate(env_variable = factor(env_variable,
                                 levels = facet_order[facet_order %in% unique(env_variable)]))

  # Create effect size annotations for each facet
  # Extract coefficient estimates and confidence intervals for each variable
  effect_annotations = regression_coefs %>%
    mutate(base_var = gsub("²", "", var_clean)) %>%
    group_by(base_var) %>%
    summarize(
      effect_text = paste0(
        paste(
          # Create intuitive labels for terms
          ifelse(grepl("²", var_clean), "Quadratic", "Linear"),
          ": ",
              sprintf("%.3f", estimate),
              " [", sprintf("%.3f", lower_ci), ", ", sprintf("%.3f", upper_ci), "]",
              sep = "", collapse = "\n")
      ),
      .groups = "drop"
    ) %>%
    mutate(
      env_variable = case_when(
        grepl("et.annual", base_var) ~ {
          term_count = sum(grepl("et.annual", regression_coefs$var_clean))
          if(term_count == 2) "Annual evapotranspiration (linear + quadratic)"
          else if(grepl("²", regression_coefs$var_clean[grepl("et.annual", regression_coefs$var_clean)][1]))
            "Annual evapotranspiration (quadratic only)"
          else "Annual evapotranspiration (linear only)"
        },
        grepl("land_use", base_var) ~ {
          term_count = sum(grepl("land_use", regression_coefs$var_clean))
          if(term_count == 2) "Land use intensity (linear + quadratic)"
          else if(grepl("²", regression_coefs$var_clean[grepl("land_use", regression_coefs$var_clean)][1]))
            "Land use intensity (quadratic only)"
          else "Land use intensity (linear only)"
        },
        grepl("summer_temp", base_var) ~ {
          term_count = sum(grepl("summer_temp", regression_coefs$var_clean))
          if(term_count == 2) "Summer temperature (linear + quadratic)"
          else if(grepl("²", regression_coefs$var_clean[grepl("summer_temp", regression_coefs$var_clean)][1]))
            "Summer temperature (quadratic only)"
          else "Summer temperature (linear only)"
        },
        TRUE ~ base_var
      )
    ) %>%
    mutate(env_variable = factor(env_variable,
                                 levels = facet_order[facet_order %in% unique(env_variable)]))

  # Create scatter plots
  p_codist_beta_scatter = ggplot() +
    geom_point(data = points_data, aes(x = beta, y = codist),
               alpha = 0.5, size = 2.5, color = "gray40") +
    geom_line(data = curve_data, aes(x = beta, y = predicted),
              color = "#00bd89", linewidth = 1.2, alpha = 0.8) +
    geom_text(data = effect_annotations,
              aes(x = -Inf, y = Inf, label = effect_text),
              hjust = -0.1, vjust = 2, size = 3.5, fontface = "bold",
              color = "black") +
    facet_wrap(~env_variable, scales = "free_x") +
    labs(
      x = "Regression coefficient\n(Summer temperature)",
      y = "Species associations\n(Variance explained)"
    ) +
    theme_bw() +
    theme(
      text = element_text(size = 14),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 13),
      strip.text = element_text(size = 15, face = "bold"),
      strip.background = element_rect(fill = "white"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 17, face = "bold"),
      plot.subtitle = element_text(size = 14),
      legend.position = "none"
    )
  print(p_codist_beta_scatter)
  cat("Scatter plots created\n")
} else {
  cat("\nNo significant predictors found\n")
  p_codist_beta_scatter = NULL
}

# ==============================================================================
# PREPARE DATA FOR MULTIVARIATE ANALYSIS
# ==============================================================================
cat("\n=== Preparing data for multivariate regression ===\n")

# Reload sites with full environmental data for regression
sites =
  res$internals$Sites %>%
  rownames_to_column("plot_id_releve") %>%
  left_join(veg.env %>% as.data.frame() %>% rownames_to_column("plot_id_releve"))

# Calculate species environmental means 
veg.env.raw = 
  veg.clim %>% 
  group_by(plot_id_releve, Longitude, Latitude, altitude)%>%
  summarize(across(slope:land_use, ~mean(.)))%>%
  ungroup()%>%
  column_to_rownames(var = "plot_id_releve")%>%
  na.omit(.)%>%
  select(-Longitude, -Latitude)

sp_env = calculate_species_env_means(species_data = Y, env_data = veg.env.raw)

species =
  res$internals$Species %>%
  rownames_to_column("species")%>%
  left_join(sp_env)


# ==============================================================================
# PCA OF ENVIRONMENTAL VARIABLES
# ==============================================================================
cat("\n=== Principal Component Analysis of Environmental Variables ===\n")

# Prepare environmental data for PCA (select only numeric environmental variables)
env_for_pca = veg.env.raw %>%
  select(slope, summer_temp, fdd, et.annual, soil_depth_mean, soil_depth_var,
         trees_cover, shrubs_cover, rocks_cover, flowdir, tpi, roughness,
         snow_sum, land_use) %>%
  drop_na()

# Perform PCA with scaling
pca_temp = PCA(env_for_pca %>% select(summer_temp, fdd, et.annual))
pca_m.clim = PCA(env_for_pca %>% select(soil_depth_var, soil_depth_mean, trees_cover, shrubs_cover, rocks_cover))
pca_topo = PCA(env_for_pca %>% select(slope, roughness, tpi, flowdir))

p_pca_temp = fviz_pca_var(pca_temp, col.var = "contrib",
                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                  repel = TRUE # Avoid text overlapping
)+labs(color = "Contribution")+ggtitle("Macroclimate")
p_pca_mclim = fviz_pca_var(pca_m.clim, col.var = "contrib",
                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                  repel = TRUE # Avoid text overlapping
)+labs(color = "Contribution")+ggtitle("Microclimate")
p_pca_topo = fviz_pca_var(pca_topo, col.var = "contrib",
                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                  repel = TRUE # Avoid text overlapping
)+labs(color = "Contribution")+ggtitle("Topography")

# Combine all three PCAs
p_pca_combined = p_pca_temp + p_pca_mclim + p_pca_topo +
  plot_annotation(title = "PCA of Environmental Variables",
                  theme = theme(plot.title = element_text(size = 16, face = "bold")))

pc_temp = as.data.frame(pca_temp$ind$coord) %>%
  select(Dim.1)%>%
  rename(pc_temp = Dim.1)%>%
  rownames_to_column("plot_id_releve")

pc_mclim = as.data.frame(pca_m.clim$ind$coord) %>%
  select(Dim.1)%>%
  rename(pc_mclim = Dim.1)%>%
  rownames_to_column("plot_id_releve")

pc_topo = as.data.frame(pca_topo$ind$coord) %>%
  select(Dim.1)%>%
  rename(pc_topo = Dim.1)%>%
  rownames_to_column("plot_id_releve")
pcs = left_join(pc_temp, pc_mclim) %>% left_join(pc_topo)

# ==============================================================================
# MULTIVARIATE REGRESSION ANALYSIS: SITES + ENVIRONMENT
# ==============================================================================
cat("\n=== Model 1a: Sites ~ Environment ===\n")

# Prepare data for sites + environment
sites_env_data = sites %>%
  drop_na()%>%
  left_join(veg %>% select(plot_id_releve, open_close))%>%
  left_join(pcs)

# Fit qgam models for each response variable (env, spa, codist) and habitat type
cat("\n=== Fitting QGAM models (median regression) ===\n")

# List to store models
qgam_models = list()
response_vars = c("env", "spa", "codist")
habitat_types = c("open", "close")

# Fit model for each response and habitat type
for(habitat in habitat_types) {
  cat(sprintf("\n--- Fitting models for %s habitat ---\n", habitat))

  # Subset data by habitat
  data_subset = sites_env_data %>% filter(open_close == habitat)

  for(resp in response_vars) {
    cat(sprintf("Fitting model for %s in %s habitat...\n", resp, habitat))

    model_name = paste(resp, habitat, sep = "_")

    qgam_models[[model_name]] = qgam::qgam(
      as.formula(paste(resp, "~ pc_temp + pc_topo + pc_mclim + land_use")),
      data = data_subset,
      qu = 0.5,
      control = list(progress = "none")
    )

    cat(sprintf("Model for %s in %s habitat fitted successfully.\n", resp, habitat))
  }
}

# Extract and plot effect sizes
cat("\n--- Extracting coefficients and confidence intervals ---\n")

# Extract coefficients, SEs, and p-values from each qgam model
coef_list = list()

for(habitat in habitat_types) {
  for(resp in response_vars) {
    model_name = paste(resp, habitat, sep = "_")
    model = qgam_models[[model_name]]

    # Get summary table with coefficients, SEs, and p-values
    model_summary = summary(model)
    coef_table = model_summary$p.table

    # Extract values
    coef_df = data.frame(
      var = rownames(coef_table),
      response = resp,
      habitat = habitat,
      Estimate = coef_table[, "Estimate"],
      Std.Error = coef_table[, "Std. Error"],
      pvalue = coef_table[, "Pr(>|z|)"],
      row.names = NULL
    )

    # Calculate 95% confidence intervals manually
    coef_df$lower_ci = coef_df$Estimate - 1.96 * coef_df$Std.Error
    coef_df$upper_ci = coef_df$Estimate + 1.96 * coef_df$Std.Error

    coef_list[[model_name]] = coef_df
  }
}

# Combine all responses
coef_data = do.call(rbind, coef_list) %>%
  filter(var != "(Intercept)") %>%
  mutate(
    significant = pvalue < 0.05,
    color_group = ifelse(significant, response, "nonsig"),
    response_label = factor(response,
                           levels = c("env", "spa", "codist"),
                           labels = c("Environment", "Space", "Species associations")),
    habitat_label = factor(habitat,
                          levels = c("open", "close"),
                          labels = c("Open", "Close"))
  ) %>%
  group_by(var, habitat) %>%
  filter(any(significant == TRUE)) %>%
  ungroup()%>%
  mutate(var = case_when(var == "pc_temp"~ "Environmental gradient",
                         var == "pc_mclim" ~ "Microclimate",
                         var == "pc_topo" ~ "Topography",
                         var == "land_use"~ "Land use"
                         ))
  

cat("\n--- Creating effect size plot with confidence intervals ---\n")

p_effect_sizes = ggplot(coef_data %>% filter(response_label != "Space"), aes(x = reorder(var, abs(Estimate)), y = Estimate, color = color_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),
                width = 0.3, position = position_dodge(width = 0.6)) +
  facet_grid(habitat_label ~ response_label) +
  coord_flip() +
  scale_color_manual(values = c(
    "env" = rgb(red = 0,blue =1, green =0),
    "spa" = rgb(red = 1,blue =0, green =0),
    "codist" = rgb(red = 0,blue =0, green =1),
    "nonsig" = "grey70"
  )) +
  labs(
    y = "Coefficient Estimate (± 95% CI)",
    x = "",
    title = "Effect Sizes from QGAM Models by Habitat Type",
    subtitle = "Grey = non-significant (p > 0.05); Colored = significant"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank(),
    text = element_text(size = 15)
  )

print(p_effect_sizes)

# Extract fitted values vs predictors
cat("\n--- Creating fitted vs predictors plot ---\n")

# Get fitted values from each qgam model with unique names
fitted_vals_list = list()

for(habitat in habitat_types) {
  data_subset = sites_env_data %>% filter(open_close == habitat)

  fitted_vals_list[[habitat]] = data.frame(
    fitted_env = fitted(qgam_models[[paste("env", habitat, sep = "_")]]),
    fitted_spa = fitted(qgam_models[[paste("spa", habitat, sep = "_")]]),
    fitted_codist = fitted(qgam_models[[paste("codist", habitat, sep = "_")]]),
    habitat = habitat,
    data_subset %>% select(pc_temp, pc_mclim, pc_topo, land_use)
  )
}

fitted_vals = do.call(rbind, fitted_vals_list)

# Combine fitted values with environmental predictors
fitted_pred_data = fitted_vals %>%
  pivot_longer(cols = starts_with("fitted_"),
               names_to = "response",
               names_prefix = "fitted_",
               values_to = "fitted") %>%
  pivot_longer(cols = c(pc_temp, pc_mclim, pc_topo, land_use),
               names_to = "predictor",
               values_to = "predictor_value") %>%
  mutate(
    response = factor(response,
                     levels = c("env", "spa", "codist"),
                     labels = c("Environment", "Space", "Species associations")),
    habitat_label = factor(habitat,
                          levels = c("open", "close"),
                          labels = c("Open", "Close"))
  )%>%
  mutate(predictor = case_when(predictor == "pc_temp"~ "Environmental gradient",
                               predictor == "pc_mclim" ~ "Microclimate",
                               predictor== "pc_topo" ~ "Topography",
                               predictor == "land_use"~ "Land use"
  ))

# Filter to only show significant predictors for cleaner visualization
sig_predictors = coef_data %>%
  filter(significant == TRUE) %>%
  pull(var) %>%
  unique()

fitted_pred_data_sig = fitted_pred_data %>%
  filter(predictor %in% sig_predictors)

# Create fitted vs predictors plot
p_fitted_pred = ggplot(fitted_pred_data_sig, aes(x = predictor_value, y = fitted, color = response)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  facet_grid(habitat_label + response ~ predictor, scales = "free") +
  scale_color_manual(values = c(
    "Environment" = rgb(red = 0,blue =1, green =0),
    "Space" = rgb(red = 1,blue =0, green =0),
    "Species associations" = rgb(red = 0,blue =0, green =1)
  )) +
  labs(
    x = "Predictor value",
    y = "Fitted values",
    title = "Fitted Values vs Environmental Predictors by Habitat Type",
    subtitle = "Showing only significant predictors"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8)
  )

print(p_fitted_pred)

cat("QGAM analysis completed\n")
# ==============================================================================
# SAVE ALL PLOTS TO SINGLE PDF
# ==============================================================================
cat("\n=== Saving all plots to single PDF ===\n")

pdf("plot/calanda_jsdm_communities_results.pdf", height = 10, width = 14)

# 1. Combined variance partitioning (Venn + ternary plots)
create_variance_combined()

# 2. Species beta coefficients heatmap
print(p_species_heatmap)

# 2. Multiple regression: Codist R² vs beta effect sizes
print(p_codist_beta_effects)

# 3. Scatter plots for significant predictors (if any exist)
if(!is.null(p_codist_beta_scatter)) {
  print(p_codist_beta_scatter)
}

# 4. Ternary plot
print(p_ternary)

# 5. Contour plot with DEM
print(p_contour)

# 6. Variance partitioning
print(p_variance)

# 7. Histogram
print(p_histogram)

# 8. Combined PCA plots (all three together)
print(p_pca_combined)

# 9. QGAM effect sizes by habitat
print(p_effect_sizes)

# 10. Fitted vs predictors
print(p_fitted_pred)

# 11. Multivariate regression results (if traits were analyzed)
if(exists("p_sites_traits")) {
  print(p_sites_traits)
}
if(exists("p_species_env")) {
  print(p_species_env)
}
if(exists("p_species_traits")) {
  print(p_species_traits)
}

dev.off()

cat("\n=== Analysis complete! ===\n")
cat("All plots saved to: plot/calanda_jsdm_communities_results.pdf\n")
cat("\nPlots included:\n")
cat("  1. Combined variance partitioning (Venn diagram + species/community ternary plots)\n")
cat("  2. Species beta coefficients heatmap\n")
cat("  3. Multiple regression: Codist R² vs beta coefficients (effect sizes)\n")
if(!is.null(p_codist_beta_scatter)) {
  cat("  4. Scatter plots for significant predictors (with partial effects)\n")
  plot_num = 5
} else {
  plot_num = 4
}
cat(sprintf("  %d. Ternary plot of variance components\n", plot_num))
plot_num = plot_num + 1
cat(sprintf("  %d. DEM contour plot with RGB points\n", plot_num))
plot_num = plot_num + 1
cat(sprintf("  %d. Variance partitioning by habitat type\n", plot_num))
plot_num = plot_num + 1
cat(sprintf("  %d. Histogram of variance components\n", plot_num))
plot_num = plot_num + 1
cat(sprintf("  %d. Combined PCA of environmental variables\n", plot_num))
plot_num = plot_num + 1
cat(sprintf("  %d. QGAM effect sizes by habitat type\n", plot_num))
plot_num = plot_num + 1
cat(sprintf("  %d. Fitted values vs predictors\n", plot_num))
plot_num = plot_num + 1
if(exists("p_sites_traits")) {
  cat(sprintf("  %d. Sites ~ Traits model\n", plot_num))
  plot_num = plot_num + 1
}
if(exists("p_species_env")) {
  cat(sprintf("  %d. Species ~ Environment model\n", plot_num))
  plot_num = plot_num + 1
}
if(exists("p_species_traits")) {
  cat(sprintf("  %d. Species ~ Traits model\n", plot_num))
}

# # Traits ----
# if(file.exists("output/traits.csv")){
#   traits = read_csv("output/traits.csv")%>%
#     select(-`...1`)
# }else{
#   try =  read_csv("data/traits/try.quantitative_traits_2025-03-17.csv")
#   indicators = read_csv("data/traits/indicators_cleaned_calanda_2025-03-17.csv")
#   dispersal = read_csv("data/traits/dispersal_cleaned_calanda_2025-03-17.csv")
#   
#   traits =
#     try %>%
#     filter(Trait %in% c("N_percent", "LDMC", "vegetative_height", "SLA", "LA"))%>%
#     filter(species %in% colnames(Y))%>%
#     filter(!is.na(Value), !is.na(Trait))%>%
#     filter(!Climate_code %in% c("Af", "Am", "Aw", "BWh", "BWk", "BSh", "BSk"))%>%
#     group_by(species, species_TNRS, Trait)%>%
#     summarize(Value = mean(Value, na.rm = TRUE))%>%
#     pivot_wider(names_from = Trait, values_from = Value)%>%
#     ungroup()
#   
#   traits = 
#     left_join(dispersal %>% select(-`...1`), traits)%>%
#     left_join(indicators %>% select(-`...1`))%>%
#     select(species, seed_mass, dispersal_distance_class, LA, SLA, LDMC, vegetative_height, N_percent,
#            Light, Moisture, Nutrients, Disturbance.Severity
#     )%>%
#     rename(disturbance = Disturbance.Severity,
#            dispersal = dispersal_distance_class)
#   
#   imp.traits = impute_functional_traits(
#     traits,
#     variables_to_impute = c("seed_mass", "dispersal", "LA", "SLA", "LDMC", "vegetative_height", 
#                             "N_percent", "Light", "Moisture", "Nutrients", "disturbance"),
#     m = 30,
#     maxiter = 100,
#     num.trees = 500,
#     seed = 123,
#     validation_fraction = 0.3
#   )
#   
#   # Check the percentage of missing values for each trait and imputation performance
#   traits %>%
#     summarise(across(everything(), ~mean(is.na(.)) * 100)) %>%
#     pivot_longer(cols = everything(), 
#                  names_to = "variable", 
#                  values_to = "percent_missing") %>%
#     arrange(desc(percent_missing)) %>%
#     right_join(imp.traits$performance) %>%
#     mutate(percent_missing = round(percent_missing),
#            r_squared = round(r_squared, 2)) %>%
#     select(variable, percent_missing, r_squared) %>%
#     gt()
#   
#   traits = 
#     imp.traits$imputed_data %>%
#     select(species, contains("final"))%>%
#     rename_with(~str_remove(., "_final$"), ends_with("_final"))
#   write.csv(traits, file = "output/traits.csv")
# }
# 
# # Get CM traits
# traits_scaled = 
#   traits %>%
#   mutate(across(c(seed_mass, LA, SLA, LDMC, vegetative_height, N_percent), ~as.numeric(scale(log(.))))) %>%
#   mutate(across(c(dispersal, Light, Nutrients, Moisture, disturbance), ~as.numeric(scale(.))))
# 
# community_traits = calculate_community_traits_extended(
#   community_data = Y,
#   traits_data = traits_scaled,
#   veg_data = veg.clim,
#   species_col = "species",
#   abundance_col = "rel_cover"
# )
# 
# sp_trait_diff = community_traits$species_abs_trait_differences %>%
#   group_by(species, trait)%>%
#   summarize(trait_diff = mean(abs_trait_difference))%>%
#   ungroup()%>%
#   pivot_wider(names_from = trait, values_from = trait_diff, names_glue = "trait_diff_{trait}")
# 
# # ==============================================================================
# # MULTIVARIATE REGRESSION ANALYSIS: SITES + TRAITS
# # ==============================================================================
# cat("\n=== Model 1b: Sites + Traits (CWM + Variance) ===\n")
# 
# # Prepare sites data with community traits
# sites_traits =
#   res$internals$Sites %>%
#   rownames_to_column("plot_id_releve") %>%
#   left_join(community_traits$community_means)
# 
# # Prepare data for sites + traits
# sites_traits_data = sites_traits %>%
#   select(env, spa, codist,
#          seed_mass, vegetative_height, dispersal, SLA, LDMC, N_percent,
#          Moisture, Nutrients, disturbance, Light,
#          seed_mass_var, vegetative_height_var, dispersal_var, SLA_var,
#          LDMC_var, N_percent_var, Moisture_var, Nutrients_var,
#          disturbance_var, Light_var) %>%
#   drop_na()
# 
# pca_trait = PCA(sites_traits_data %>% select(seed_mass, vegetative_height, dispersal, SLA, LDMC, N_percent, disturbance, Light))
# 
# pca_trait_var = PCA(sites_traits_data %>% select(seed_mass_var, vegetative_height_var,
#                                                  LDMC_var, N_percent_var))
# 
# # Create mvabund object with raw variance components
# Y_sites_traits = mvabund(sites_traits_data[, c("env", "spa", "codist")])
# 
# # Fit full model with CWM and variance
# fit_sites_traits_full = manylm(Y_sites_traits ~ seed_mass + vegetative_height + dispersal +
#                                                  SLA + LDMC + N_percent + Moisture + Nutrients +
#                                                  disturbance + Light + seed_mass_var +
#                                                  vegetative_height_var + dispersal_var + SLA_var +
#                                                  LDMC_var + N_percent_var + Moisture_var +
#                                                  Nutrients_var + disturbance_var + Light_var,
#                                data = sites_traits_data)
# 
# cat("\n--- ANOVA with adjusted p-values ---\n")
# anova_sites_traits = anova.manylm(fit_sites_traits_full, p.uni = "adjusted")
# print(anova_sites_traits)
# 
# # ==============================================================================
# # MULTIVARIATE REGRESSION ANALYSIS: SPECIES + ENVIRONMENT
# # ==============================================================================
# cat("\n=== Model 2a: Species + Environment ===\n")
# 
# # Prepare species data with environmental means
# species_env =
#   res$internals$Species %>%
#   rownames_to_column("species") %>%
#   left_join(sp_env)
# 
# # Prepare data for species + environment
# species_env_data = species_env %>%
#   select(env, spa, codist, slope:land_use) %>%
#   drop_na()
# 
# # Create mvabund object with raw variance components
# Y_species_env = mvabund(species_env_data[, c("env", "spa", "codist")])
# 
# # Fit full model
# fit_species_env_full = manylm(Y_species_env ~ slope + summer_temp + fdd + et.annual +
#                                               soil_depth_mean + soil_depth_var + trees_cover +
#                                               shrubs_cover + rocks_cover + flowdir + tpi +
#                                               roughness + snow_sum + land_use,
#                               data = species_env_data)
# 
# cat("\n--- ANOVA with adjusted p-values ---\n")
# anova_species_env = anova.manylm(fit_species_env_full, p.uni = "adjusted")
# print(anova_species_env)
# 
# 
# # ==============================================================================
# # MULTIVARIATE REGRESSION ANALYSIS: SPECIES + TRAITS
# # ==============================================================================
# cat("\n=== Model 2b: Species + Traits ===\n")
# 
# # Prepare species data with traits
# species_traits =
#   res$internals$Species %>%
#   rownames_to_column("species") %>%
#   left_join(traits_scaled)
# 
# # Prepare data for species + traits
# species_traits_data = species_traits %>%
#   select(env, spa, codist,
#          seed_mass, vegetative_height, dispersal, SLA, LDMC, N_percent,
#          Moisture, Nutrients, disturbance, Light) %>%
#   drop_na()
# 
# pca_trait_sp = PCA(species_traits_data %>% select(seed_mass, vegetative_height, SLA, LDMC, N_percent, disturbance, dispersal))
# 
# # Create mvabund object with raw variance components
# Y_species_traits = mvabund(species_traits_data[, c("env", "spa", "codist")])
# 
# # Fit full model with species traits only
# fit_species_traits_full = manylm(Y_species_traits ~ seed_mass + vegetative_height + dispersal +
#                                                      SLA + LDMC + N_percent + Moisture + Nutrients +
#                                                      disturbance + Light,
#                                  data = species_traits_data)
# 
# # Backward selection
# cat("\nRunning backward selection...\n")
# fit_species_traits_step = stepMvabund(fit_species_traits_full, direction = "backward", trace = 2)
# 
# # Final model summary
# cat("\n--- Final Model Summary ---\n")
# print(summary(fit_species_traits_step))
# 
# cat("\n--- ANOVA with adjusted p-values ---\n")
# anova_species_traits = anova.manylm(fit_species_traits_step, p.uni = "adjusted")
# print(anova_species_traits)
# 
# # ==============================================================================
# # VISUALIZATION OF MULTIVARIATE RESULTS
# # ==============================================================================
# 
# # Function to extract and plot coefficients from manylm models
# extract_coef_plot = function(model, title_text) {
#   # Get coefficients
#   coef_summary = summary(model)
# 
#   # Extract coefficients for each response variable
#   coef_list = lapply(names(coef_summary$coefficients), function(response) {
#     coef_df = as.data.frame(coef_summary$coefficients[[response]])
#     coef_df$var = rownames(coef_df)
#     coef_df$response = response
#     return(coef_df)
#   })
# 
#   # Combine into one data frame
#   coef_data = do.call(rbind, coef_list)
#   colnames(coef_data)[1:4] = c("Estimate", "Std.Error", "t.value", "Pr")
# 
#   # Remove intercept
#   coef_data = coef_data %>%
#     filter(var != "(Intercept)") %>%
#     mutate(
#       lower_ci = Estimate - 1.96 * Std.Error,
#       upper_ci = Estimate + 1.96 * Std.Error,
#       significant = Pr < 0.05
#     )
# 
#   # Create plot
#   p = ggplot(coef_data, aes(x = var, y = Estimate, color = response)) +
#     geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
#     geom_point(size = 3, aes(alpha = significant), position = position_dodge(width = 0.6)) +
#     geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci, alpha = significant),
#                   width = 0.3, position = position_dodge(width = 0.6)) +
#     facet_wrap(~response, scales = "free_x") +
#     coord_flip() +
#     scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3)) +
#     scale_color_manual(values = c(
#       "env" = rgb(0.3, 0.3, 1),
#       "spa" = rgb(1, 0.3, 0.3),
#       "codist" = rgb(0.3, 1, 0.3)
#     )) +
#     labs(
#       y = "Coefficient Estimate",
#       x = "Predictor Variable",
#       title = title_text,
#       subtitle = "Transparent points/bars indicate non-significant effects (p > 0.05)"
#     ) +
#     theme_bw() +
#     theme(
#       legend.position = "bottom",
#       strip.background = element_rect(fill = "white"),
#       panel.grid.minor = element_blank()
#     ) +
#     guides(alpha = "none")
# 
#   return(p)
# }
# 
# # Create plots for each model
# cat("\nCreating visualization plots...\n")
# 
# # Note: Using full models since stepwise selection was only done for species_traits
# if(exists("fit_sites_traits_full")) {
#   p_sites_traits = extract_coef_plot(fit_sites_traits_full,
#                                       "Model 1b: Sites ~ Traits (CWM + Variance)")
# }
# 
# if(exists("fit_species_env_full")) {
#   p_species_env = extract_coef_plot(fit_species_env_full,
#                                      "Model 2a: Species ~ Environment")
# }
# 
# if(exists("fit_species_traits_step")) {
#   p_species_traits = extract_coef_plot(fit_species_traits_step,
#                                         "Model 2b: Species ~ Traits")
# }
# 
# # ==============================================================================
# # SAVE ALL PLOTS TO SINGLE PDF
# # ==============================================================================
# cat("\n=== Saving all plots to single PDF ===\n")
# 
# pdf("plot/calanda_jsdm_mvabund_results.pdf", height = 10, width = 14)
# 
# # 1. Ternary plot
# print(p_ternary)
# 
# # 2. Contour plot with DEM
# print(p_contour)
# 
# # 3. Variance partitioning
# print(p_variance)
# 
# # 4. Histogram
# print(p_histogram)
# 
# # 5. Combined PCA plots (all three together)
# print(p_pca_combined)
# 
# # 6. QGAM effect sizes by habitat
# print(p_effect_sizes)
# 
# # 7. Fitted vs predictors
# print(p_fitted_pred)
# 
# # 8. Multivariate regression results (if traits were analyzed)
# if(exists("p_sites_traits")) {
#   print(p_sites_traits)
# }
# if(exists("p_species_env")) {
#   print(p_species_env)
# }
# if(exists("p_species_traits")) {
#   print(p_species_traits)
# }
# 
# dev.off()
# 
# cat("\n=== Analysis complete! ===\n")
# cat("All plots saved to: plot/calanda_jsdm_mvabund_results.pdf\n")
# cat("\nPlots included:\n")
# cat("  1. Ternary plot of variance components\n")
# cat("  2. DEM contour plot with RGB points\n")
# cat("  3. Variance partitioning by habitat type\n")
# cat("  4. Histogram of variance components\n")
# cat("  5. Combined PCA of environmental variables\n")
# cat("  6. QGAM effect sizes by habitat type\n")
# cat("  7. Fitted values vs predictors\n")
# if(exists("p_sites_traits")) cat("  8. Sites ~ Traits model\n")
# if(exists("p_species_env")) cat("  9. Species ~ Environment model\n")
# if(exists("p_species_traits")) cat(" 10. Species ~ Traits model\n")
