# ==============================================================================
# Script: 06_species_postJSDM.R
# Purpose: Species-level environmental effect sizes on co-distribution variance
#
# Inputs:
#   - results_from_Max/model_sjsdm_calanda.rds
#   - results_from_Max/res_sjsdm_calanda.rds
#   - R/00_setup/functions_calanda.R (label_env_var)
#
# Outputs:
#   - plot/species_regressions.pdf
# ==============================================================================

library(tidyverse)
library(sjSDM)
library(ggrepel)
library(here)

source(here("Calanda_JSDM", "R", "00_setup", "functions_calanda.R"))

# Load data ----
model_jsdm = readRDS(here("Calanda_JSDM", "results_from_Max", "model_sjsdm_calanda.rds"))
environment(model_jsdm$get_model)$device = "cpu"
res = readRDS(here("Calanda_JSDM", "results_from_Max", "res_sjsdm_calanda.rds"))

# ==============================================================================
# SPECIES-LEVEL ENVIRONMENTAL EFFECT SIZES POST-JSDM
# ==============================================================================
cat("\n=== Species-level environmental effects analysis ===\n")

# Extract coefficients and variance components
model_summary = summary(model_jsdm)
species_betas = as.data.frame(t(model_summary$coefs[-1, ])) %>%
  rownames_to_column("species")

species_variance = res$internals$Species %>%
  rownames_to_column("species") %>%
  mutate(
    total = env + codist + spa,
    codist_prop = codist / total
  )

# Prepare regression data
regression_data = species_betas %>%
  left_join(species_variance %>% select(species, codist, codist_prop, env, spa), by = "species") %>%
  drop_na()

cat(sprintf("Species with complete data: %d\n", nrow(regression_data)))

# Remove outliers (0.1% and 99.9% quantiles)
remove_outliers = function(data, var) {
  Q1 = quantile(data[[var]], 0.001)
  Q3 = quantile(data[[var]], 0.999)
  IQR_val = Q3 - Q1
  data %>% filter(!!sym(var) >= Q1 - 1.5 * IQR_val & !!sym(var) <= Q3 + 1.5 * IQR_val)
}

regression_data_clean = regression_data %>%
  remove_outliers("codist") %>%
  remove_outliers("codist_prop")

cat(sprintf("After outlier removal: %d species\n", nrow(regression_data_clean)))

# Fit models with stepwise selection
env_vars = c("summer_temp", "fdd", "et_annual", "land_use")
available_vars = env_vars[env_vars %in% names(regression_data_clean)]

cat(sprintf("\nPredictors: %s\n", paste(available_vars, collapse = ", ")))

models = list(
  codist_raw = step(lm(as.formula(paste("codist ~",
                                         paste(c(available_vars,
                                                paste0("I(", available_vars, "^2)")),
                                              collapse = " + "))),
                      data = regression_data_clean),
                   direction = "backward", trace = 0),
  codist_prop = step(lm(as.formula(paste("codist_prop ~",
                                          paste(c(available_vars,
                                                 paste0("I(", available_vars, "^2)")),
                                               collapse = " + "))),
                        data = regression_data_clean),
                    direction = "backward", trace = 0)
)

cat("\nFinal models fitted\n")

# Helper: generate predictions
predict_model = function(model, env_var, data) {
  coef_summary = summary(model)$coefficients
  linear = env_var
  quad = paste0("I(", env_var, "^2)")

  if (!linear %in% rownames(coef_summary) && !quad %in% rownames(coef_summary)) return(NULL)

  beta_seq = seq(min(data[[env_var]]), max(data[[env_var]]), length.out = 100)
  pred = coef(model)[1]

  if (linear %in% rownames(coef_summary)) {
    pred = pred + coef(model)[linear] * beta_seq
  }
  if (quad %in% rownames(coef_summary)) {
    pred = pred + coef(model)[quad] * beta_seq^2
  }

  data.frame(
    x = beta_seq,
    y = pred,
    env_var = env_var,
    significant = (linear %in% rownames(coef_summary) && coef_summary[linear, 4] < 0.05) ||
                  (quad %in% rownames(coef_summary) && coef_summary[quad, 4] < 0.05)
  )
}

# Generate predictions
predictions = bind_rows(lapply(names(models), function(m) {
  bind_rows(lapply(available_vars, function(v) {
    pred = predict_model(models[[m]], v, regression_data_clean)
    if (is.null(pred)) return(NULL)
    pred$model = m
    pred
  }))
})) %>%
  separate(model, c("variance", "type"), sep = "_") %>%
  mutate(
    model_type = ifelse(type == "prop", "Dominance", "Raw"),
    env_variable_label = label_env_var(env_var)
  )

# Filter significant combinations
sig_combos = predictions %>%
  group_by(env_variable_label, model_type) %>%
  filter(any(significant)) %>%
  ungroup()

# Prepare plot data
plot_data = regression_data_clean %>%
  pivot_longer(c(codist, codist_prop), names_to = "measure", values_to = "value") %>%
  pivot_longer(all_of(available_vars), names_to = "env_var", values_to = "beta") %>%
  mutate(
    model_type = ifelse(grepl("prop", measure), "Dominance", "Raw"),
    env_variable_label = label_env_var(env_var)
  ) %>%
  semi_join(sig_combos %>% select(env_variable_label, model_type) %>% distinct(),
            by = c("env_variable_label", "model_type"))

# Extract coefficients for labels
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
    model_type = ifelse(type == "prop", "Dominance", "Raw"),
    env_variable_label = label_env_var(base_var)
  ) %>%
  semi_join(sig_combos %>% select(env_variable_label, model_type) %>% distinct(),
            by = c("env_variable_label", "model_type")) %>%
  group_by(env_variable_label, model_type) %>%
  summarize(
    effect_text = paste(paste0(var_type, ": ", sprintf("%.3f", estimate),
                               ifelse(significant, "*", "")),
                        collapse = ", "),
    .groups = "drop"
  )

# Position labels
effect_pos = sig_combos %>%
  group_by(env_variable_label, model_type) %>%
  summarize(x_pos = min(x), y_pos = y[which.min(x)], .groups = "drop") %>%
  left_join(effect_labels, by = c("env_variable_label", "model_type"))

model_pos = sig_combos %>%
  group_by(env_variable_label, model_type) %>%
  summarize(x_pos = max(x), y_pos = y[which.max(x)], .groups = "drop")

# Plot
p_species = ggplot() +
  geom_point(data = plot_data, aes(x = beta, y = value),
             alpha = 0.4, size = 1.5, color = "gray40") +
  geom_line(data = sig_combos, aes(x = x, y = y, linetype = model_type, alpha = significant),
            color = "#00bd89", linewidth = 1.2) +
  geom_text_repel(data = effect_pos, aes(x = x_pos, y = y_pos, label = effect_text),
                  size = 3.5, fontface = "bold", color = "#00bd89",
                  show.legend = FALSE, box.padding = 0.5, max.overlaps = 20) +
  geom_text_repel(data = model_pos, aes(x = x_pos, y = y_pos, label = model_type),
                  size = 3.5, fontface = "italic", color = "#00bd89",
                  show.legend = FALSE, box.padding = 0.5, max.overlaps = 20) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5), guide = "none") +
  scale_linetype_manual(values = c("Dominance" = "solid", "Raw" = "dashed")) +
  facet_wrap(~ env_variable_label, scales = "free") +
  labs(x = "Species regression coefficient\n(environmental effect size)",
       y = "Species associations\n(raw or dominance)",
       linetype = "Variance type") +
  theme_bw() +
  theme(text = element_text(size = 13),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 13, face = "bold"),
        legend.position = "bottom")

pdf(here("Calanda_JSDM", "plot", "species_regressions.pdf"), width = 10, height = 6)
print(p_species)
dev.off()
cat("Species-level scatter plots created\n")

cat("\n=== Species-level analysis complete ===\n")
