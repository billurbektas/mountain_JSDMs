# ==============================================================================
# Script: vp_diagnostic_figures.R
# Author: Billur Bektas
# Claude (Anthropic) was used to assist with code refactoring, validation, and documentation.
#
# Purpose: Diagnostic figures comparing proportional vs discard variance
#          allocation at the per-unit level, and visualizing the distribution
#          of raw shared fractions. Uses one fold for illustration.
#
# Outputs:
#   - plot/vp_diagnostic_discard_vs_proportional.pdf
#   - plot/vp_diagnostic_raw_shared_fractions.pdf
#
# Requires: R/00_setup/functions_calanda.R
# ==============================================================================

library(tidyverse)
library(conflicted)
library(here)
library(patchwork)

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)

source(here("Calanda_JSDM", "R", "00_setup", "functions_calanda.R"))

# --- Colors ---
color_env    = "#81caf3"
color_spa    = "#d00000"
color_codist = "#00bd89"

# ==============================================================================
# LOAD DATA — use fold 1 as representative
# ==============================================================================
cat("\n=== VP diagnostic figures (fold 1) ===\n")

fold_file = sort(list.files(
  here("Calanda_JSDM", "output", "results", "cv"),
  pattern = "^fold_.*\\.rds$", full.names = TRUE
))[1]

fd = readRDS(fold_file)
cat(sprintf("  Loaded: %s\n", basename(fold_file)))

# ==============================================================================
# EXTRACT PER-UNIT DATA
# ==============================================================================

# --- Sites ---
raw_si = fd$partition$anova$sites$R2_McFadden
prop_si = fd$partition$anova$sites$R2_McFadden_shared$proportional
site_names = names(raw_si$F_A)
n_sites = length(site_names)

df_sites = tibble(
  unit = site_names,
  level = "Communities",
  # Raw 7 fractions
  F_A = raw_si$F_A, F_B = raw_si$F_B, F_S = raw_si$F_S,
  F_AB = raw_si$F_AB, F_AS = raw_si$F_AS, F_BS = raw_si$F_BS,
  F_ABS = raw_si$F_ABS, Full = raw_si$Full,
  # Proportional allocation
  prop_env = prop_si$F_A, prop_assoc = prop_si$F_B, prop_space = prop_si$F_S,
  # Discard = unique only, floored
  disc_env = pmax(raw_si$F_A, 0),
  disc_assoc = pmax(raw_si$F_B, 0),
  disc_space = pmax(raw_si$F_S, 0)
)

# --- Species ---
raw_sp = fd$partition$anova$species$R2_McFadden
prop_sp = fd$partition$anova$species$R2_McFadden_shared$proportional
sp_names = names(raw_sp$F_A)
n_species = length(sp_names)

df_species = tibble(
  unit = sp_names,
  level = "Species",
  F_A = raw_sp$F_A, F_B = raw_sp$F_B, F_S = raw_sp$F_S,
  F_AB = raw_sp$F_AB, F_AS = raw_sp$F_AS, F_BS = raw_sp$F_BS,
  F_ABS = raw_sp$F_ABS, Full = raw_sp$Full,
  prop_env = prop_sp$F_A, prop_assoc = prop_sp$F_B, prop_space = prop_sp$F_S,
  disc_env = pmax(raw_sp$F_A, 0),
  disc_assoc = pmax(raw_sp$F_B, 0),
  disc_space = pmax(raw_sp$F_S, 0)
)

df_all = bind_rows(df_sites, df_species)

cat(sprintf("  Sites: %d, Species: %d\n", n_sites, n_species))

# ==============================================================================
# FIGURE 1: Discard vs Proportional — scatter plots
# ==============================================================================
cat("  Creating discard vs proportional scatter plots...\n")

build_scatter = function(df, comp_disc, comp_prop, comp_label, comp_color) {
  r = cor(df[[comp_disc]], df[[comp_prop]], use = "complete.obs")

  ggplot(df, aes(x = .data[[comp_disc]], y = .data[[comp_prop]])) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_hline(yintercept = 0, linewidth = 0.2, color = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.2, color = "grey70") +
    geom_point(alpha = 0.5, size = 1.2, color = comp_color) +
    annotate("text", x = min(df[[comp_disc]], na.rm = TRUE),
             y = max(df[[comp_prop]], na.rm = TRUE),
             label = sprintf("r = %.2f", r),
             hjust = 0, vjust = 1, size = 4, fontface = "italic") +
    labs(x = paste0(comp_label, " (discard)"),
         y = paste0(comp_label, " (proportional)")) +
    theme_bw(base_size = 11) +
    theme(axis.title = element_text(size = 12))
}

# Species
p_sp_env = build_scatter(df_species, "disc_env", "prop_env", "Environment", color_env)
p_sp_assoc = build_scatter(df_species, "disc_assoc", "prop_assoc", "Species Associations", color_codist)
p_sp_space = build_scatter(df_species, "disc_space", "prop_space", "Space", color_spa)

# Sites
p_si_env = build_scatter(df_sites, "disc_env", "prop_env", "Environment", color_env)
p_si_assoc = build_scatter(df_sites, "disc_assoc", "prop_assoc", "Species Associations", color_codist)
p_si_space = build_scatter(df_sites, "disc_space", "prop_space", "Space", color_spa)

# Add titles to first row
p_sp_env = p_sp_env + ggtitle("A) Species") +
  theme(plot.title = element_text(size = 14, face = "bold"))
p_si_env = p_si_env + ggtitle("B) Communities") +
  theme(plot.title = element_text(size = 14, face = "bold"))

p_scatter = (p_sp_env | p_sp_assoc | p_sp_space) /
            (p_si_env | p_si_assoc | p_si_space)

pdf(here("Calanda_JSDM", "plot", "vp_diagnostic_discard_vs_proportional.pdf"),
    width = 14, height = 9)
print(p_scatter)
dev.off()
cat("  Saved vp_diagnostic_discard_vs_proportional.pdf\n")

# ==============================================================================
# FIGURE 2: Raw shared fractions — violins
# ==============================================================================
cat("  Creating raw shared fraction violins...\n")

# Build long format of raw 7 fractions
df_raw_long = df_all %>%
  select(unit, level, F_A, F_B, F_S, F_AB, F_AS, F_BS, F_ABS) %>%
  pivot_longer(cols = F_A:F_ABS, names_to = "fraction", values_to = "value") %>%
  mutate(
    fraction_type = case_when(
      fraction %in% c("F_A", "F_B", "F_S") ~ "Unique",
      fraction %in% c("F_AB", "F_AS", "F_BS") ~ "Pairwise shared",
      fraction == "F_ABS" ~ "Triple shared"
    ),
    fraction_label = case_when(
      fraction == "F_A" ~ "Env",
      fraction == "F_B" ~ "Assoc",
      fraction == "F_S" ~ "Space",
      fraction == "F_AB" ~ "Env x Assoc",
      fraction == "F_AS" ~ "Env x Space",
      fraction == "F_BS" ~ "Assoc x Space",
      fraction == "F_ABS" ~ "All three"
    ),
    fraction_label = factor(fraction_label,
      levels = c("Env", "Assoc", "Space", "Env x Assoc", "Env x Space", "Assoc x Space", "All three")),
    fraction_color = case_when(
      fraction == "F_A" ~ color_env,
      fraction == "F_B" ~ color_codist,
      fraction == "F_S" ~ color_spa,
      fraction == "F_AB" ~ "#5090a0",
      fraction == "F_AS" ~ "#b06060",
      fraction == "F_BS" ~ "#608060",
      fraction == "F_ABS" ~ "grey50"
    )
  )

# Compute % negative per fraction per level
neg_stats = df_raw_long %>%
  group_by(level, fraction_label) %>%
  summarise(
    n_neg = sum(value < 0),
    n_total = n(),
    pct_neg = 100 * n_neg / n_total,
    median_val = median(value),
    .groups = "drop"
  )

# Violin + jitter
frac_colors = c("Env" = color_env, "Assoc" = color_codist, "Space" = color_spa,
                "Env x Assoc" = "#5090a0", "Env x Space" = "#b06060",
                "Assoc x Space" = "#608060", "All three" = "grey50")

p_violins = ggplot(df_raw_long, aes(x = fraction_label, y = value, fill = fraction_label)) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "black", linetype = "dashed") +
  geom_violin(alpha = 0.3, color = NA) +
  geom_jitter(aes(color = fraction_label), width = 0.15, size = 0.3, alpha = 0.3) +
  # Add % negative annotation
  geom_text(data = neg_stats %>% filter(pct_neg > 0),
            aes(x = fraction_label, y = min(df_raw_long$value) - 0.05,
                label = sprintf("%.0f%% neg", pct_neg)),
            size = 3, color = "red3", inherit.aes = FALSE) +
  scale_fill_manual(values = frac_colors, guide = "none") +
  scale_color_manual(values = frac_colors, guide = "none") +
  facet_wrap(~ level, ncol = 1, scales = "free_y") +
  labs(x = NULL, y = "McFadden R\u00B2 fraction",
       title = "Raw variance fractions before proportional allocation",
       subtitle = "Dashed line = 0. Red annotations = % of units with negative values.") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 10),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10)
  )

pdf(here("Calanda_JSDM", "plot", "vp_diagnostic_raw_shared_fractions.pdf"),
    width = 10, height = 10)
print(p_violins)
dev.off()
cat("  Saved vp_diagnostic_raw_shared_fractions.pdf\n")

# ==============================================================================
# PRINT SUMMARY TABLE
# ==============================================================================
cat("\n=== Summary: Global anova vs per-unit means ===\n")

an = fd$partition$anova$results
r2_global = setNames(an[["R2 McFadden"]], an$models)

cat(sprintf("\n  %-20s %10s %10s %10s %10s\n", "", "Global", "Per-site", "Per-species", ""))
cat(sprintf("  %-20s %10s %10s %10s %10s\n", "", "anova", "prop mean", "prop mean", ""))

# Global proportional allocation
FA = r2_global["F_A"]; FB = r2_global["F_B"]; FS = r2_global["F_S"]
FAB = r2_global["F_AB"]; FAS = r2_global["F_AS"]; FBS = r2_global["F_BS"]; FABS = r2_global["F_ABS"]
st = abs(FAB) + abs(FBS) + abs(FAS)
FAB2 = FAB + FABS * abs(FAB)/st; FBS2 = FBS + FABS * abs(FBS)/st; FAS2 = FAS + FABS * abs(FAS)/st
g_env = FA + FAB2 * abs(FA)/(abs(FA)+abs(FB)) + FAS2 * abs(FA)/(abs(FS)+abs(FA))
g_assoc = FB + FAB2 * abs(FB)/(abs(FA)+abs(FB)) + FBS2 * abs(FB)/(abs(FS)+abs(FB))
g_space = FS + FAS2 * abs(FS)/(abs(FS)+abs(FA)) + FBS2 * abs(FS)/(abs(FS)+abs(FB))

cat(sprintf("  %-20s %10.4f %10.4f %10.4f\n", "Env (proportional):", g_env, mean(prop_si$F_A), mean(prop_sp$F_A)))
cat(sprintf("  %-20s %10.4f %10.4f %10.4f\n", "Assoc (proportional):", g_assoc, mean(prop_si$F_B), mean(prop_sp$F_B)))
cat(sprintf("  %-20s %10.4f %10.4f %10.4f\n", "Space (proportional):", g_space, mean(prop_si$F_S), mean(prop_sp$F_S)))
cat(sprintf("  %-20s %10.4f %10.4f %10.4f\n", "Env (discard):", r2_global["F_A"], mean(df_sites$disc_env), mean(df_species$disc_env)))
cat(sprintf("  %-20s %10.4f %10.4f %10.4f\n", "Assoc (discard):", r2_global["F_B"], mean(df_sites$disc_assoc), mean(df_species$disc_assoc)))
cat(sprintf("  %-20s %10.4f %10.4f %10.4f\n", "Space (discard):", r2_global["F_S"], mean(df_sites$disc_space), mean(df_species$disc_space)))

# ==============================================================================
# FIGURE 3 & 4: ALL 10 FOLDS AGGREGATED
# ==============================================================================
cat("\n=== Aggregating across all 10 folds ===\n")

fold_files = sort(list.files(
  here("Calanda_JSDM", "output", "results", "cv"),
  pattern = "^fold_.*\\.rds$", full.names = TRUE
))

all_sites_list = list()
all_species_list = list()

for (fi in seq_along(fold_files)) {
  fd_i = readRDS(fold_files[fi])

  raw_si_i = fd_i$partition$anova$sites$R2_McFadden
  prop_si_i = fd_i$partition$anova$sites$R2_McFadden_shared$proportional
  si_names_i = names(raw_si_i$F_A)

  all_sites_list[[fi]] = tibble(
    unit = si_names_i, fold = fi, level = "Communities",
    F_A = raw_si_i$F_A, F_B = raw_si_i$F_B, F_S = raw_si_i$F_S,
    F_AB = raw_si_i$F_AB, F_AS = raw_si_i$F_AS, F_BS = raw_si_i$F_BS,
    F_ABS = raw_si_i$F_ABS, Full = raw_si_i$Full,
    prop_env = prop_si_i$F_A, prop_assoc = prop_si_i$F_B, prop_space = prop_si_i$F_S,
    disc_env = pmax(raw_si_i$F_A, 0),
    disc_assoc = pmax(raw_si_i$F_B, 0),
    disc_space = pmax(raw_si_i$F_S, 0)
  )

  raw_sp_i = fd_i$partition$anova$species$R2_McFadden
  prop_sp_i = fd_i$partition$anova$species$R2_McFadden_shared$proportional
  sp_names_i = names(raw_sp_i$F_A)

  all_species_list[[fi]] = tibble(
    unit = sp_names_i, fold = fi, level = "Species",
    F_A = raw_sp_i$F_A, F_B = raw_sp_i$F_B, F_S = raw_sp_i$F_S,
    F_AB = raw_sp_i$F_AB, F_AS = raw_sp_i$F_AS, F_BS = raw_sp_i$F_BS,
    F_ABS = raw_sp_i$F_ABS, Full = raw_sp_i$Full,
    prop_env = prop_sp_i$F_A, prop_assoc = prop_sp_i$F_B, prop_space = prop_sp_i$F_S,
    disc_env = pmax(raw_sp_i$F_A, 0),
    disc_assoc = pmax(raw_sp_i$F_B, 0),
    disc_space = pmax(raw_sp_i$F_S, 0)
  )
}

df_all_sites = bind_rows(all_sites_list)
df_all_species = bind_rows(all_species_list)

cat(sprintf("  Aggregated: %d site-fold observations, %d species-fold observations\n",
            nrow(df_all_sites), nrow(df_all_species)))

# Average across folds per unit
df_sites_avg = df_all_sites %>%
  group_by(unit, level) %>%
  summarise(across(F_A:disc_space, ~ mean(.x, na.rm = TRUE)), .groups = "drop")

df_species_avg = df_all_species %>%
  group_by(unit, level) %>%
  summarise(across(F_A:disc_space, ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# ==============================================================================
# FIGURE 3: Discard vs Proportional — all folds averaged
# ==============================================================================
cat("  Creating all-folds scatter plots...\n")

p_sp_env_all = build_scatter(df_species_avg, "disc_env", "prop_env", "Environment", color_env)
p_sp_assoc_all = build_scatter(df_species_avg, "disc_assoc", "prop_assoc", "Species Associations", color_codist)
p_sp_space_all = build_scatter(df_species_avg, "disc_space", "prop_space", "Space", color_spa)

p_si_env_all = build_scatter(df_sites_avg, "disc_env", "prop_env", "Environment", color_env)
p_si_assoc_all = build_scatter(df_sites_avg, "disc_assoc", "prop_assoc", "Species Associations", color_codist)
p_si_space_all = build_scatter(df_sites_avg, "disc_space", "prop_space", "Space", color_spa)

p_sp_env_all = p_sp_env_all + ggtitle("A) Species (10-fold mean)") +
  theme(plot.title = element_text(size = 14, face = "bold"))
p_si_env_all = p_si_env_all + ggtitle("B) Communities (10-fold mean)") +
  theme(plot.title = element_text(size = 14, face = "bold"))

p_scatter_all = (p_sp_env_all | p_sp_assoc_all | p_sp_space_all) /
                (p_si_env_all | p_si_assoc_all | p_si_space_all)

pdf(here("Calanda_JSDM", "plot", "vp_diagnostic_discard_vs_proportional_allfolds.pdf"),
    width = 14, height = 9)
print(p_scatter_all)
dev.off()
cat("  Saved vp_diagnostic_discard_vs_proportional_allfolds.pdf\n")

# ==============================================================================
# FIGURE 4: Raw shared fractions — all folds
# ==============================================================================
cat("  Creating all-folds raw shared fraction violins...\n")

df_all_combined = bind_rows(df_all_sites, df_all_species)

df_raw_long_all = df_all_combined %>%
  select(unit, fold, level, F_A, F_B, F_S, F_AB, F_AS, F_BS, F_ABS) %>%
  pivot_longer(cols = F_A:F_ABS, names_to = "fraction", values_to = "value") %>%
  mutate(
    fraction_label = case_when(
      fraction == "F_A" ~ "Env",
      fraction == "F_B" ~ "Assoc",
      fraction == "F_S" ~ "Space",
      fraction == "F_AB" ~ "Env x Assoc",
      fraction == "F_AS" ~ "Env x Space",
      fraction == "F_BS" ~ "Assoc x Space",
      fraction == "F_ABS" ~ "All three"
    ),
    fraction_label = factor(fraction_label,
      levels = c("Env", "Assoc", "Space", "Env x Assoc", "Env x Space", "Assoc x Space", "All three"))
  )

neg_stats_all = df_raw_long_all %>%
  group_by(level, fraction_label) %>%
  summarise(
    n_neg = sum(value < 0),
    n_total = n(),
    pct_neg = 100 * n_neg / n_total,
    .groups = "drop"
  )

p_violins_all = ggplot(df_raw_long_all, aes(x = fraction_label, y = value, fill = fraction_label)) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "black", linetype = "dashed") +
  geom_violin(alpha = 0.3, color = NA) +
  geom_jitter(aes(color = fraction_label), width = 0.15, size = 0.15, alpha = 0.15) +
  geom_text(data = neg_stats_all %>% filter(pct_neg > 0),
            aes(x = fraction_label, y = min(df_raw_long_all$value) - 0.05,
                label = sprintf("%.0f%% neg", pct_neg)),
            size = 3, color = "red3", inherit.aes = FALSE) +
  scale_fill_manual(values = frac_colors, guide = "none") +
  scale_color_manual(values = frac_colors, guide = "none") +
  facet_wrap(~ level, ncol = 1, scales = "free_y") +
  labs(x = NULL, y = "McFadden R\u00B2 fraction",
       title = "Raw variance fractions before proportional allocation (all 10 folds)",
       subtitle = "Each unit appears 9-10 times (once per fold it was in the training set).") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 10),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10)
  )

pdf(here("Calanda_JSDM", "plot", "vp_diagnostic_raw_shared_fractions_allfolds.pdf"),
    width = 10, height = 10)
print(p_violins_all)
dev.off()
cat("  Saved vp_diagnostic_raw_shared_fractions_allfolds.pdf\n")

# ==============================================================================
# PRINT ALL-FOLDS SUMMARY
# ==============================================================================
cat("\n=== All-folds summary: per-unit means (averaged across folds) ===\n")
cat(sprintf("  %-25s %10s %10s\n", "", "Per-site", "Per-species"))
cat(sprintf("  %-25s %10.4f %10.4f\n", "Env (proportional):", mean(df_sites_avg$prop_env), mean(df_species_avg$prop_env)))
cat(sprintf("  %-25s %10.4f %10.4f\n", "Assoc (proportional):", mean(df_sites_avg$prop_assoc), mean(df_species_avg$prop_assoc)))
cat(sprintf("  %-25s %10.4f %10.4f\n", "Space (proportional):", mean(df_sites_avg$prop_space), mean(df_species_avg$prop_space)))
cat(sprintf("  %-25s %10.4f %10.4f\n", "Env (discard):", mean(df_sites_avg$disc_env), mean(df_species_avg$disc_env)))
cat(sprintf("  %-25s %10.4f %10.4f\n", "Assoc (discard):", mean(df_sites_avg$disc_assoc), mean(df_species_avg$disc_assoc)))
cat(sprintf("  %-25s %10.4f %10.4f\n", "Space (discard):", mean(df_sites_avg$disc_space), mean(df_species_avg$disc_space)))

cat(sprintf("\n  Raw fractions %% negative (all folds):\n"))
neg_stats_all %>%
  filter(pct_neg > 0) %>%
  arrange(level, fraction_label) %>%
  mutate(label = sprintf("    %s - %-15s %5.1f%% negative (%d/%d)", level, fraction_label, pct_neg, n_neg, n_total)) %>%
  pull(label) %>%
  cat(sep = "\n")

cat("\n\n=== VP diagnostic figures complete ===\n")
