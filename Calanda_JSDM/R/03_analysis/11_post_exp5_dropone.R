# ==============================================================================
# Script: 11_post_exp5_dropone.R
# Purpose: Post-analysis of sjSDM Experiment 5 drop-one environmental variable
#          Compare each drop-one model to the BASE (full) model:
#            - Model-level anova R2 differences (all 8 Venn fractions)
#            - Species-level: mean of per-species differences + correlations
#            - Site-level: mean of per-site differences + correlations
#
# Inputs:
#   - output/results/dropone/dropone_BASE_<param_tag>.rds
#   - output/results/dropone/dropone_*_<param_tag>.rds (12 drop-one runs)
#   - param_tag format: a<alpha>_l<lambda>_le<lambda_env>
#
# Outputs:
#   - plot/exp5_anova_diff_<param_tag>.pdf
#   - plot/exp5_species_diff_<param_tag>.pdf
#   - plot/exp5_sites_diff_<param_tag>.pdf
#   - plot/exp5_correlations_<param_tag>.pdf
#   - output/results/exp5_dropone_summary_<param_tag>.csv
#
# Requires: R/00_setup/functions_calanda.R
# ==============================================================================

library(tidyverse)
library(conflicted)
library(here)
library(patchwork)

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::select)

source(here("Calanda_JSDM", "R", "00_setup", "functions_calanda.R"))

# Colors
color_env     = "#81caf3"
color_spa     = "#d00000"
color_codist  = "#00bd89"
color_overall = "black"

color_env_bio = "#41a4ae"
color_env_spa = "#b0656a"
color_bio_spa = "#6b5f45"
color_shared  = "#9b59b6"

anova_colors = c(
  "Overall"     = color_overall,
  "Environment" = color_env,
  "Spatial"     = color_spa,
  "Biotic"      = color_codist,
  "Env x Bio"   = color_env_bio,
  "Env x Spa"   = color_env_spa,
  "Bio x Spa"   = color_bio_spa,
  "Shared"      = color_shared
)

unit_colors = c(
  "Overall"     = color_overall,
  "Environment" = color_env,
  "Spatial"     = color_spa,
  "Biotic"      = color_codist
)

dir.create(here("Calanda_JSDM", "plot"), showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# STEP 1: LOAD DATA
# ==============================================================================
cat("\n=== Step 1: Loading drop-one results ===\n")

# Detect param_tag from available files
all_dropone_files = list.files(
  here("Calanda_JSDM", "output", "results", "dropone"),
  pattern = "^dropone_.*\\.rds$", full.names = TRUE
)

tag_match = regmatches(basename(all_dropone_files[1]),
                       regexpr("a[0-9.]+_l[0-9.]+_le[0-9.]+(?=\\.)", basename(all_dropone_files[1]), perl = TRUE))
if (length(tag_match) == 0 || tag_match == "") {
  stop("Could not detect param_tag from filenames. Expected format: dropone_*_a<alpha>_l<lambda>_le<lambda_env>.rds")
}
param_tag = tag_match
cat(sprintf("  Detected param_tag: %s\n", param_tag))

# Filter to matching param_tag
tagged_files = all_dropone_files[grepl(param_tag, basename(all_dropone_files), fixed = TRUE)]

base_file = tagged_files[grepl("dropone_BASE", basename(tagged_files))]
base = readRDS(base_file)
cat(sprintf("  BASE R2: %.6f\n", base$partition$R2))

# Extract hyperparameters
hp_alpha      = base$alpha
hp_lambda     = base$lambda
hp_lambda_env = base$lambda_env
hp_subtitle   = sprintf("alpha = %s, lambda = %s, lambda_env = %s | fit sampling = %s, anova samples = %s",
                         hp_alpha, hp_lambda, hp_lambda_env,
                         base$fit_sampling, base$anova_samples)
cat(sprintf("  Hyperparameters: %s\n", hp_subtitle))

drop_files = tagged_files[!grepl("dropone_BASE", basename(tagged_files))]
cat(sprintf("  Found %d drop-one files\n", length(drop_files)))

# Extract anova R2 from a partition object
get_anova_r2 = function(partition) {
  anova_res = partition$anova$results
  r2 = setNames(anova_res$`R2 McFadden`, anova_res$models)
  tibble(
    metric = c("Overall", "Environment", "Spatial", "Biotic",
               "Env x Bio", "Env x Spa", "Bio x Spa", "Shared"),
    value  = c(r2["Full"], r2["F_A"], r2["F_S"], r2["F_B"],
               r2["F_AB"], r2["F_AS"], r2["F_BS"], r2["F_ABS"])
  )
}

# Base anova
base_anova = get_anova_r2(base$partition)

# Process each drop-one run
anova_diff_list   = list()
species_diff_list = list()
sites_diff_list   = list()
corr_list         = list()

for (f in drop_files) {
  run = readRDS(f)
  dropped = run$dropped

  # --- Anova differences ---
  drop_anova = get_anova_r2(run$partition)
  anova_diff_list[[dropped]] = tibble(
    dropped  = dropped,
    metric   = base_anova$metric,
    base_val = base_anova$value,
    drop_val = drop_anova$value,
    diff     = drop_anova$value - base_anova$value
  )

  # --- Species-level: per-species differences then mean ---
  sp_base = base$partition$species  # 352 x 4
  sp_drop = run$partition$species
  sp_diff = sp_drop - sp_base  # per-species difference

  species_diff_list[[dropped]] = tibble(
    dropped   = dropped,
    component = c("Overall", "Environment", "Spatial", "Biotic"),
    col       = c("r2", "env", "spa", "codist"),
    mean_diff = c(mean(sp_diff[, "r2"], na.rm = TRUE),
                  mean(sp_diff[, "env"], na.rm = TRUE),
                  mean(sp_diff[, "spa"], na.rm = TRUE),
                  mean(sp_diff[, "codist"], na.rm = TRUE))
  )

  # --- Site-level: per-site differences then mean ---
  si_base = base$partition$sites  # 576 x 4
  si_drop = run$partition$sites
  si_diff = si_drop - si_base

  sites_diff_list[[dropped]] = tibble(
    dropped   = dropped,
    component = c("Overall", "Environment", "Spatial", "Biotic"),
    col       = c("r2", "env", "spa", "codist"),
    mean_diff = c(mean(si_diff[, "r2"], na.rm = TRUE),
                  mean(si_diff[, "env"], na.rm = TRUE),
                  mean(si_diff[, "spa"], na.rm = TRUE),
                  mean(si_diff[, "codist"], na.rm = TRUE))
  )

  # --- Correlations: base vs drop-one per species/site ---
  corr_rows = list()
  for (col in c("r2", "env", "spa", "codist")) {
    comp_label = recode(col, r2 = "Overall", env = "Environment",
                        spa = "Spatial", codist = "Biotic")
    corr_rows = c(corr_rows, list(tibble(
      dropped     = dropped,
      component   = comp_label,
      species_r   = cor(sp_base[, col], sp_drop[, col], use = "complete.obs"),
      sites_r     = cor(si_base[, col], si_drop[, col], use = "complete.obs")
    )))
  }
  corr_list[[dropped]] = bind_rows(corr_rows)

  cat(sprintf("  %s: R2 diff = %+.6f\n", dropped,
              drop_anova$value[drop_anova$metric == "Overall"] - base_anova$value[base_anova$metric == "Overall"]))
}

df_anova_diff   = bind_rows(anova_diff_list)
df_species_diff = bind_rows(species_diff_list)
df_sites_diff   = bind_rows(sites_diff_list)
df_corr         = bind_rows(corr_list)

# ==============================================================================
# STEP 2: SUMMARY TABLE
# ==============================================================================
cat("\n=== Step 2: Summary table ===\n")

# Combine into one summary
df_summary = df_anova_diff %>%
  select(dropped, metric, diff) %>%
  pivot_wider(names_from = metric, values_from = diff, names_prefix = "anova_")

df_sp_wide = df_species_diff %>%
  select(dropped, component, mean_diff) %>%
  pivot_wider(names_from = component, values_from = mean_diff, names_prefix = "species_mean_diff_")

df_si_wide = df_sites_diff %>%
  select(dropped, component, mean_diff) %>%
  pivot_wider(names_from = component, values_from = mean_diff, names_prefix = "sites_mean_diff_")

df_corr_wide = df_corr %>%
  pivot_wider(names_from = component,
              values_from = c(species_r, sites_r),
              names_glue = "{.value}_{component}")

df_full_summary = df_summary %>%
  left_join(df_sp_wide, by = "dropped") %>%
  left_join(df_si_wide, by = "dropped") %>%
  left_join(df_corr_wide, by = "dropped") %>%
  mutate(alpha = hp_alpha, lambda = hp_lambda, lambda_env = hp_lambda_env, .before = 1)

summary_file = paste0("exp5_dropone_summary_", param_tag, ".csv")
write_csv(df_full_summary,
          here("Calanda_JSDM", "output", "results", summary_file))
cat(sprintf("  Saved %s\n", summary_file))

# ==============================================================================
# STEP 3: ANOVA R2 DIFFERENCES PLOT
# ==============================================================================
cat("\n=== Step 3: Anova R2 differences plot ===\n")

df_anova_diff = df_anova_diff %>%
  mutate(metric = factor(metric, levels = names(anova_colors)))

# Order variables by overall R2 difference (most negative first = most important)
var_order = df_anova_diff %>%
  filter(metric == "Overall") %>%
  arrange(diff) %>%
  pull(dropped)

df_anova_diff = df_anova_diff %>%
  mutate(dropped = factor(dropped, levels = var_order))

p_anova = ggplot(df_anova_diff, aes(x = dropped, y = diff, fill = metric)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = anova_colors) +
  labs(
    title = "Experiment 5: Model-level anova R\u00B2 change when dropping each variable",
    subtitle = paste0("Difference = drop-one minus BASE (negative = variable contributes positively)\n", hp_subtitle),
    x = "Dropped variable", y = expression(Delta * " McFadden R\u00B2"),
    fill = "Component"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )

pdf(here("Calanda_JSDM", "plot", paste0("exp5_anova_diff_", param_tag, ".pdf")),
    width = 14, height = 7)
print(p_anova)
dev.off()
cat(sprintf("  Saved exp5_anova_diff_%s.pdf\n", param_tag))

# ==============================================================================
# STEP 4: SPECIES-LEVEL MEAN DIFFERENCES
# ==============================================================================
cat("\n=== Step 4: Species mean differences plot ===\n")

df_species_diff = df_species_diff %>%
  mutate(
    component = factor(component, levels = names(unit_colors)),
    dropped   = factor(dropped, levels = var_order)
  )

n_species = nrow(base$partition$species)

p_species = ggplot(df_species_diff, aes(x = dropped, y = mean_diff, fill = component)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = unit_colors) +
  labs(
    title = "Experiment 5: Mean per-species difference (drop-one minus BASE)",
    subtitle = paste0("Average across ", n_species, " species of (drop_i - base_i) per component\n", hp_subtitle),
    x = "Dropped variable", y = "Mean difference",
    fill = "Component"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

pdf(here("Calanda_JSDM", "plot", paste0("exp5_species_diff_", param_tag, ".pdf")),
    width = 14, height = 7)
print(p_species)
dev.off()
cat(sprintf("  Saved exp5_species_diff_%s.pdf\n", param_tag))

# ==============================================================================
# STEP 5: SITE-LEVEL MEAN DIFFERENCES
# ==============================================================================
cat("\n=== Step 5: Site mean differences plot ===\n")

df_sites_diff = df_sites_diff %>%
  mutate(
    component = factor(component, levels = names(unit_colors)),
    dropped   = factor(dropped, levels = var_order)
  )

n_sites = nrow(base$partition$sites)

p_sites = ggplot(df_sites_diff, aes(x = dropped, y = mean_diff, fill = component)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = unit_colors) +
  labs(
    title = "Experiment 5: Mean per-site difference (drop-one minus BASE)",
    subtitle = paste0("Average across ", n_sites, " sites of (drop_i - base_i) per component\n", hp_subtitle),
    x = "Dropped variable", y = "Mean difference",
    fill = "Component"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

pdf(here("Calanda_JSDM", "plot", paste0("exp5_sites_diff_", param_tag, ".pdf")),
    width = 14, height = 7)
print(p_sites)
dev.off()
cat(sprintf("  Saved exp5_sites_diff_%s.pdf\n", param_tag))

# ==============================================================================
# STEP 6: CORRELATIONS (BASE vs DROP-ONE)
# ==============================================================================
cat("\n=== Step 6: Correlations plot ===\n")

df_corr_long = df_corr %>%
  mutate(
    component = factor(component, levels = names(unit_colors)),
    dropped   = factor(dropped, levels = var_order)
  ) %>%
  pivot_longer(cols = c(species_r, sites_r),
               names_to = "level", values_to = "correlation") %>%
  mutate(level = recode(level, species_r = "Species", sites_r = "Sites"))

p_corr = ggplot(df_corr_long,
                aes(x = dropped, y = correlation, fill = component)) +
  geom_hline(yintercept = 1, linewidth = 0.3, linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = unit_colors) +
  facet_wrap(~level, ncol = 1) +
  labs(
    title = "Experiment 5: Correlation between BASE and drop-one (per species/site)",
    subtitle = paste0("Pearson r across ", n_species, " species or ", n_sites, " sites per component\n", hp_subtitle),
    x = "Dropped variable", y = "Pearson r",
    fill = "Component"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )

pdf(here("Calanda_JSDM", "plot", paste0("exp5_correlations_", param_tag, ".pdf")),
    width = 14, height = 10)
print(p_corr)
dev.off()
cat(sprintf("  Saved exp5_correlations_%s.pdf\n", param_tag))

cat("\n=== Post Experiment 5 analysis complete ===\n")
