# ==============================================================================
# Script: 05_community_postJSDM.R
# Purpose: Community-level environmental regressions on variance components
#
# Inputs:
#   - output/starter_data_25.04.25.RData (res, veg.env)
#   - R/00_setup/functions_calanda.R (label_env_var)
#
# Outputs:
#   - plot/community_regressions.pdf
#   - plot/spatial_topography_regressions.pdf
# ==============================================================================

library(tidyverse)
library(ggrepel)
library(here)

source(here("Calanda_JSDM", "R", "00_setup", "functions_calanda.R"))

# Load data ----
load(here("Calanda_JSDM", "output", "starter_data_25.04.25.RData"))
res = readRDS(here("Calanda_JSDM", "results_from_Max", "res_sjsdm_calanda.rds"))

# Rename dotted objects from RData to snake_case
veg_env = veg.env
rm(veg.env)

# ==============================================================================
# COMMUNITY-LEVEL ENVIRONMENTAL EFFECTS POST-JSDM
# ==============================================================================
cat("\n=== Community-level environmental effects analysis ===\n")

# Prepare regression data with proportions
regression_data = res$internals$Sites %>%
  rownames_to_column("plot_id_releve") %>%
  left_join(veg_env %>% as.data.frame() %>% rownames_to_column("plot_id_releve")) %>%
  select(plot_id_releve, env, codist, spa, summer_temp, fdd, et_annual, land_use) %>%
  drop_na() %>%
  mutate(
    total = env + codist + spa,
    env_prop = env / total,
    codist_prop = codist / total,
    spa_prop = spa / total
  )

cat(sprintf("\nSites with complete data: %d\n", nrow(regression_data)))

# Fit all models with stepwise AIC selection
models = list(
  env_prop = step(lm(env_prop ~ summer_temp + I(summer_temp^2) + fdd + I(fdd^2) +
                              et_annual + I(et_annual^2) + land_use + I(land_use^2),
                    data = regression_data), direction = "backward", trace = 0),
  codist_prop = step(lm(codist_prop ~ summer_temp + I(summer_temp^2) + fdd + I(fdd^2) +
                                   et_annual + I(et_annual^2) + land_use + I(land_use^2),
                       data = regression_data), direction = "backward", trace = 0),
  env_raw = step(lm(env ~ summer_temp + I(summer_temp^2) + fdd + I(fdd^2) +
                       et_annual + I(et_annual^2) + land_use + I(land_use^2),
                   data = regression_data), direction = "backward", trace = 0),
  codist_raw = step(lm(codist ~ summer_temp + I(summer_temp^2) + fdd + I(fdd^2) +
                             et_annual + I(et_annual^2) + land_use + I(land_use^2),
                      data = regression_data), direction = "backward", trace = 0)
)

# Helper: generate predictions for a model
predict_model = function(model, env_var, data) {
  coef_summary = summary(model)$coefficients
  linear = env_var
  quad = paste0("I(", env_var, "^2)")

  if (!linear %in% rownames(coef_summary) && !quad %in% rownames(coef_summary)) return(NULL)

  pred_df = data.frame(
    summer_temp = mean(data$summer_temp),
    fdd = mean(data$fdd),
    et_annual = mean(data$et_annual),
    land_use = mean(data$land_use)
  )[rep(1, 100), ]

  pred_df[[env_var]] = seq(min(data[[env_var]]), max(data[[env_var]]), length.out = 100)
  pred_df$y = predict(model, pred_df)
  pred_df$x = pred_df[[env_var]]
  pred_df$env_var = env_var
  pred_df$significant = (linear %in% rownames(coef_summary) && coef_summary[linear, 4] < 0.05) ||
                        (quad %in% rownames(coef_summary) && coef_summary[quad, 4] < 0.05)
  pred_df
}

# Generate all predictions
env_vars = c("summer_temp", "fdd", "et_annual", "land_use")
predictions = bind_rows(lapply(names(models), function(m) {
  bind_rows(lapply(env_vars, function(v) {
    pred = predict_model(models[[m]], v, regression_data)
    if (is.null(pred)) return(NULL)
    pred$model = m
    pred
  }))
})) %>%
  separate(model, c("variance", "type"), sep = "_") %>%
  mutate(
    variance_type = ifelse(variance == "env", "Environment", "Species associations"),
    model_type = ifelse(type == "prop", "Dominance", "Raw"),
    env_variable_label = label_env_var(env_var)
  )

# Filter: keep only env_var + variance_type with at least one significant model
sig_combos = predictions %>%
  group_by(env_variable_label, variance_type, model_type) %>%
  filter(any(significant)) %>%
  ungroup()

# Prepare plot data
plot_data = regression_data %>%
  pivot_longer(c(env_prop, codist_prop, env, codist), names_to = "measure", values_to = "value") %>%
  pivot_longer(all_of(env_vars), names_to = "env_var", values_to = "env_value") %>%
  mutate(
    variance_type = ifelse(grepl("env", measure), "Environment", "Species associations"),
    model_type = ifelse(grepl("prop", measure), "Dominance", "Raw"),
    env_variable_label = label_env_var(env_var)
  ) %>%
  filter(!is.nan(spa_prop)) %>%
  semi_join(sig_combos %>% select(env_variable_label, variance_type, model_type) %>% distinct(),
            by = c("env_variable_label", "variance_type", "model_type"))

# Extract coefficients and create effect size labels
extract_coefs = function(model_name) {
  coefs = summary(models[[model_name]])$coefficients %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    filter(term != "(Intercept)")

  if (nrow(coefs) == 0) {
    return(data.frame(base_var = character(), var_type = character(),
                     estimate = numeric(), lower_ci = numeric(),
                     upper_ci = numeric(), significant = logical(),
                     model = character()))
  }

  coefs %>%
    mutate(
      base_var = case_when(
        grepl("summer_temp", term) ~ "summer_temp",
        grepl("fdd", term) ~ "fdd",
        grepl("et_annual", term) ~ "et_annual",
        grepl("land_use", term) ~ "land_use"
      ),
      var_type = as.character(ifelse(grepl("\\^2|I\\(", term), "Q", "L")),
      estimate = Estimate,
      lower_ci = Estimate - 1.96 * `Std. Error`,
      upper_ci = Estimate + 1.96 * `Std. Error`,
      significant = `Pr(>|t|)` < 0.05,
      model = model_name
    ) %>%
    select(base_var, var_type, estimate, lower_ci, upper_ci, significant, model)
}

effect_labels = bind_rows(lapply(names(models), extract_coefs)) %>%
  separate(model, c("variance", "type"), sep = "_") %>%
  mutate(
    variance_type = ifelse(variance == "env", "Environment", "Species associations"),
    model_type = ifelse(type == "prop", "Dominance", "Raw"),
    env_variable_label = label_env_var(base_var)
  ) %>%
  semi_join(sig_combos %>% select(env_variable_label, variance_type, model_type) %>% distinct(),
            by = c("env_variable_label", "variance_type", "model_type")) %>%
  group_by(env_variable_label, variance_type, model_type) %>%
  summarize(
    effect_text = paste(paste0(var_type, ": ", sprintf("%.3f", estimate),
                               ifelse(significant, "*", "")),
                        collapse = ", "),
    .groups = "drop"
  )

# Position labels on lines
effect_pos = sig_combos %>%
  group_by(env_variable_label, variance_type, model_type) %>%
  summarize(x_pos = min(x), y_pos = y[which.min(x)], .groups = "drop") %>%
  left_join(effect_labels, by = c("env_variable_label", "variance_type", "model_type"))

model_pos = sig_combos %>%
  group_by(env_variable_label, variance_type, model_type) %>%
  summarize(x_pos = max(x), y_pos = y[which.max(x)], .groups = "drop")

# Plot
p_community = ggplot() +
  geom_point(data = plot_data, aes(x = env_value, y = value, color = variance_type),
             alpha = 0.5, size = 0.3) +
  geom_line(data = sig_combos, aes(x = x, y = y, color = variance_type,
                                    linetype = model_type, alpha = significant),
            linewidth = 1) +
  geom_text_repel(data = effect_pos, aes(x = x_pos, y = y_pos, label = effect_text,
                                          color = variance_type),
                  size = 4, fontface = "bold", show.legend = FALSE,
                  box.padding = 0.5, max.overlaps = 20) +
  geom_text_repel(data = model_pos, aes(x = x_pos, y = y_pos, label = model_type,
                                         color = variance_type),
                  size = 4, fontface = "italic", show.legend = FALSE,
                  box.padding = 0.5, max.overlaps = 20) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5), guide = "none") +
  scale_linetype_manual(values = c("Dominance" = "solid", "Raw" = "dashed")) +
  scale_color_manual(values = c("Environment" = "#81caf3", "Species associations" = "#00bd89")) +
  facet_wrap(~ env_variable_label, scales = "free") +
  labs(x = "Environmental value (scaled)",
       y = "Variance explained\n(raw or dominance)",
       color = "Component", linetype = "Variance type") +
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 11, face = "bold"),
        legend.position = "bottom")

pdf(here("Calanda_JSDM", "plot", "community_regressions.pdf"), width = 10, height = 6)
print(p_community)
dev.off()
cat("Community-level scatter plots created\n")

# ==============================================================================
# SPATIAL VARIANCE WITH TOPOGRAPHY
# ==============================================================================
cat("\n=== Spatial variance vs topography analysis ===\n")

# Prepare spatial regression data with topography
regression_data_spa = res$internals$Sites %>%
  rownames_to_column("plot_id_releve") %>%
  left_join(veg_env %>% as.data.frame() %>% rownames_to_column("plot_id_releve")) %>%
  select(plot_id_releve, env, codist, spa, slope, tpi, roughness) %>%
  drop_na() %>%
  mutate(
    total = env + codist + spa,
    spa_prop = spa / total
  )

cat(sprintf("\nSites with complete topography data: %d\n", nrow(regression_data_spa)))

# Fit spatial models with stepwise AIC selection
models_spa = list(
  spa_prop = step(lm(spa_prop ~ slope + I(slope^2) + tpi + I(tpi^2) +
                              roughness + I(roughness^2),
                    data = regression_data_spa), direction = "backward", trace = 0),
  spa_raw = step(lm(spa ~ slope + I(slope^2) + tpi + I(tpi^2) +
                       roughness + I(roughness^2),
                   data = regression_data_spa), direction = "backward", trace = 0)
)

# Helper: label topography variables
label_topo_var = function(var) {
  factor(case_when(
    var == "slope" ~ "Slope",
    var == "tpi" ~ "TPI",
    var == "roughness" ~ "Roughness"
  ), levels = c("Slope", "TPI", "Roughness"))
}

# Helper: generate predictions for topography models
predict_topo_model = function(model, topo_var, data) {
  coef_summary = summary(model)$coefficients
  linear = topo_var
  quad = paste0("I(", topo_var, "^2)")

  if (!linear %in% rownames(coef_summary) && !quad %in% rownames(coef_summary)) return(NULL)

  pred_df = data.frame(
    slope = mean(data$slope),
    tpi = mean(data$tpi),
    roughness = mean(data$roughness)
  )[rep(1, 100), ]

  pred_df[[topo_var]] = seq(min(data[[topo_var]]), max(data[[topo_var]]), length.out = 100)
  pred_df$y = predict(model, pred_df)
  pred_df$x = pred_df[[topo_var]]
  pred_df$topo_var = topo_var
  pred_df$significant = (linear %in% rownames(coef_summary) && coef_summary[linear, 4] < 0.05) ||
                        (quad %in% rownames(coef_summary) && coef_summary[quad, 4] < 0.05)
  pred_df
}

# Generate all predictions
topo_vars = c("slope", "tpi", "roughness")
predictions_spa = bind_rows(lapply(names(models_spa), function(m) {
  bind_rows(lapply(topo_vars, function(v) {
    pred = predict_topo_model(models_spa[[m]], v, regression_data_spa)
    if (is.null(pred)) return(NULL)
    pred$model = m
    pred
  }))
})) %>%
  separate(model, c("variance", "type"), sep = "_") %>%
  mutate(
    model_type = ifelse(type == "prop", "Dominance", "Raw"),
    topo_variable_label = label_topo_var(topo_var)
  )

# Keep all predictions (both significant and non-significant)
plot_combos_spa = predictions_spa

# Prepare plot data
plot_data_spa = regression_data_spa %>%
  pivot_longer(c(spa_prop, spa), names_to = "measure", values_to = "value") %>%
  pivot_longer(all_of(topo_vars), names_to = "topo_var", values_to = "topo_value") %>%
  mutate(
    model_type = ifelse(grepl("prop", measure), "Dominance", "Raw"),
    topo_variable_label = label_topo_var(topo_var)
  )

# Extract coefficients for topography models
extract_topo_coefs = function(model_name) {
  coefs = summary(models_spa[[model_name]])$coefficients %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    filter(term != "(Intercept)")

  if (nrow(coefs) == 0) {
    return(data.frame(base_var = character(), var_type = character(),
                     estimate = numeric(), lower_ci = numeric(),
                     upper_ci = numeric(), significant = logical(),
                     model = character()))
  }

  coefs %>%
    mutate(
      base_var = case_when(
        grepl("slope", term) ~ "slope",
        grepl("tpi", term) ~ "tpi",
        grepl("roughness", term) ~ "roughness"
      ),
      var_type = as.character(ifelse(grepl("\\^2|I\\(", term), "Q", "L")),
      estimate = Estimate,
      lower_ci = Estimate - 1.96 * `Std. Error`,
      upper_ci = Estimate + 1.96 * `Std. Error`,
      significant = `Pr(>|t|)` < 0.05,
      model = model_name
    ) %>%
    select(base_var, var_type, estimate, lower_ci, upper_ci, significant, model)
}

effect_labels_spa = bind_rows(lapply(names(models_spa), extract_topo_coefs)) %>%
  separate(model, c("variance", "type"), sep = "_") %>%
  mutate(
    model_type = ifelse(type == "prop", "Dominance", "Raw"),
    topo_variable_label = label_topo_var(base_var)
  ) %>%
  group_by(topo_variable_label, model_type) %>%
  summarize(
    effect_text = paste(paste0(var_type, ": ", sprintf("%.3f", estimate),
                               ifelse(significant, "*", "")),
                        collapse = ", "),
    has_sig = any(significant),
    .groups = "drop"
  )

# Position labels
effect_pos_spa = plot_combos_spa %>%
  group_by(topo_variable_label, model_type) %>%
  summarize(x_pos = min(x), y_pos = y[which.min(x)], .groups = "drop") %>%
  left_join(effect_labels_spa, by = c("topo_variable_label", "model_type"))

model_pos_spa = plot_combos_spa %>%
  group_by(topo_variable_label, model_type) %>%
  summarize(x_pos = max(x), y_pos = y[which.max(x)], .groups = "drop")

# Plot
p_spatial = ggplot() +
  geom_point(data = plot_data_spa, aes(x = topo_value, y = value),
             alpha = 0.3, size = 0.3, color = "grey60") +
  geom_line(data = plot_combos_spa, aes(x = x, y = y, linetype = model_type,
                                         color = significant, alpha = significant),
            linewidth = 1) +
  geom_text_repel(data = effect_pos_spa %>% filter(!is.na(effect_text)),
                  aes(x = x_pos, y = y_pos, label = effect_text,
                      color = has_sig),
                  size = 3.5, fontface = "bold", show.legend = FALSE,
                  box.padding = 0.5, max.overlaps = 20) +
  geom_text_repel(data = model_pos_spa, aes(x = x_pos, y = y_pos, label = model_type),
                  size = 3.5, fontface = "italic", color = "grey30",
                  show.legend = FALSE, box.padding = 0.5, max.overlaps = 20) +
  scale_color_manual(values = c("TRUE" = "#ff6b6b", "FALSE" = "grey60"),
                     labels = c("TRUE" = "Significant", "FALSE" = "Non-significant"),
                     name = "Relationship") +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4), guide = "none") +
  scale_linetype_manual(values = c("Dominance" = "solid", "Raw" = "dashed")) +
  facet_wrap(~ topo_variable_label, scales = "free") +
  labs(x = "Topographic value (scaled)",
       y = "Spatial variance\n(raw or dominance)",
       linetype = "Variance type") +
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 11, face = "bold"),
        legend.position = "bottom")

pdf(here("Calanda_JSDM", "plot", "spatial_topography_regressions.pdf"), width = 10, height = 6)
print(p_spatial)
dev.off()
cat("Spatial-topography scatter plots created\n")

cat("\n=== Community-level analysis complete ===\n")
