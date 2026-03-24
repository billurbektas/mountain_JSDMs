# ==============================================================================
# Script: 03_merge_and_assess_traits.R
# Author: Billur Bektas
# Claude (Anthropic) was used to assist with code refactoring, validation, and documentation.
#
# Purpose: Merge TRY and field trait data, calculate species-level medians and
#          kCV, and assess trait coverage per community and bias in variance
#          components.
#
# Description:
#   Combines individual-level trait observations from:
#     1. TRY database (output/try_traits_individual.csv from 01)
#     2. Field measurements (output/field_traits_clean.csv from 02)
#   Harmonizes units (cm, g/m2, mg/g, ug), converts TRY SLA to LMA (1/SLA),
#   field heights (mm→cm) and N (% → mg/g).
#   Retains traits: vegetative_height, LNC, LCC, LDMC, LMA, seed_mass, LA.
#   Note: LCC is included in species medians and median imputation but excluded
#   from kCV calculation and kCV imputation.
#   Calculates species-level medians (all observations) and kCV (= CV/(1+CV)
#   on log-transformed data, for species with >= 5 observations per trait).
#   Computes species-level trait correlations, merges with Nutrients
#   indicator and dispersal traits.
#
# Input files:
#   - output/try_traits_individual.csv (from 01_prepare_TRY_traits.R)
#   - output/field_traits_clean.csv (from 02_prepare_field_trait_data.R)
#   - output/indicators.csv (from 01_prepare_TRY_traits.R)
#   - output/dispersal.csv (from 01_prepare_TRY_traits.R)
#
# Output files:
#   - output/all_traits_individual.csv (combined individual-level data)
#   - output/traits_raw.csv (species summaries with medians, kCV, indicators)
#   - output/traits_medians_imputed.csv (imputed species medians + indicators)
#   - output/traits_medians_imputation_flags.csv (Original/Imputed per species-trait)
#   - output/traits_kcv_imputed.csv (imputed species kCV values)
#   - output/traits_kcv_imputation_flags.csv (Original/Imputed per species-trait)
#   - output/species_trait_correlations.csv
#   - plot/unit_comparison_field_vs_try.pdf
#   - plot/trait_distributions_raw_vs_log.pdf
#   - plot/species_trait_correlations.pdf
#   - plot/imputation_distributions.pdf
# ==============================================================================

library(tidyverse)
library(here)
library(corrplot)
library(FactoMineR)
library(factoextra)
library(ggrepel)
library(gt)

source(here("Calanda_JSDM", "R", "00_setup", "functions_calanda.R"))

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================
cat("\n=== Loading trait data ===\n")

# Load TRY individual-level data
try_individual = read_csv(
  here("Calanda_JSDM", "output", "try_traits_individual.csv"),
  show_col_types = FALSE
)
cat("TRY data:", nrow(try_individual), "observations\n")

# Load field individual-level data
field_individual = read_csv(
  here("Calanda_JSDM", "output", "field_traits_clean.csv"),
  show_col_types = FALSE
)
cat("Field data:", nrow(field_individual), "observations\n")

# Load indicators and dispersal
indicators = read_csv(
  here("Calanda_JSDM", "output", "indicators.csv"),
  show_col_types = FALSE
)

dispersal = read_csv(
  here("Calanda_JSDM", "output", "dispersal.csv"),
  show_col_types = FALSE
)

# ==============================================================================
# 2. DEFINE TARGET UNITS AND RETAINED TRAITS
# ==============================================================================

# Continuous traits to retain (7 traits)
continuous_traits = c("vegetative_height", "LNC", "LCC", "LDMC", "LMA", "seed_mass", "LA")

# Target units for merged traits
target_units = tibble(
  trait = c("vegetative_height", "LMA", "LA", "LDMC", "LNC", "LCC", "seed_mass"),
  unit = c("cm", "g/m2", "cm2", "mg/g", "mg/g", "mg/g", "ug")
)

cat("\n=== Target units for merged traits ===\n")
print(target_units, n = Inf)

# ==============================================================================
# 3. HARMONIZE TRY TRAIT NAMES AND CONVERT SLA → LMA
# ==============================================================================
cat("\n=== Harmonizing TRY trait names ===\n")

# TRY data: standardize names and keep only target traits
try_long = try_individual %>%
  mutate(
    trait = case_when(
      TraitName == "N_percent" ~ "LNC",
      TraitName == "C_percent" ~ "LCC",
      TraitName == "SLA" ~ "SLA",
      TraitName == "LA" ~ "LA",
      TraitName == "LDMC" ~ "LDMC",
      TraitName == "vegetative_height" ~ "vegetative_height",
      TraitName == "seed_mass" ~ "seed_mass",
      TRUE ~ NA_character_
    ),
    source = "TRY"
  ) %>%
  filter(!is.na(trait)) %>%
  select(species_TNRS, trait, value = Value, source)

# Convert TRY units to target units:
#   SLA  (mm2/mg = m2/kg) -> LMA (g/m2): LMA = 1/SLA * 1000
#   LDMC (g/g)            -> mg/g:        * 1000
#   vegetative_height (m) -> cm:          * 100
#   seed_mass (mg)        -> ug:          * 1000
#   LA   (mm2)            -> mm2:         no conversion needed
#   LNC  (mg/g)           -> mg/g:        no conversion needed
#   LCC  (mg/g)           -> mg/g:        no conversion needed
n_sla = sum(try_long$trait == "SLA", na.rm = TRUE)
n_ldmc = sum(try_long$trait == "LDMC", na.rm = TRUE)
try_long = try_long %>%
  mutate(
    # SLA (mm2/mg = m2/kg) -> LMA (g/m2): 1/SLA gives kg/m2, * 1000 gives g/m2
    value = ifelse(trait == "SLA" & value > 0, (1 / value) * 1000, value),
    trait = ifelse(trait == "SLA", "LMA", trait),
    # LDMC: g/g -> mg/g
    value = ifelse(trait == "LDMC", value * 1000, value),
    # vegetative_height: m -> cm
    value = ifelse(trait == "vegetative_height", value * 100, value),
    # seed_mass: mg -> ug
    value = ifelse(trait == "seed_mass", value * 1000, value),
    #LA: mm2 -> cm2
    value = ifelse(trait == "LA", value/100, value)
  )
cat(sprintf("Converted %d TRY SLA observations to LMA (1/SLA * 1000, mm2/mg -> g/m2)\n", n_sla))
cat(sprintf("Converted %d TRY LDMC observations from g/g to mg/g (* 1000)\n", n_ldmc))
cat("Converted TRY vegetative_height from m to cm (* 100)\n")
cat("Converted TRY seed_mass from mg to ug (* 1000)\n")
cat("Converted TRY LA from mm2 to cm2 (/100)\n")

cat("TRY traits standardized:", n_distinct(try_long$trait), "traits\n")

# ==============================================================================
# 4. CONVERT FIELD UNITS AND RESHAPE TO LONG FORMAT
# ==============================================================================
cat("\n=== Converting field units to TRY standard ===\n")

# Field unit conversions:
#   vegetative_height:   mm -> cm   (/ 10)
#   N_content_corr:      % -> mg/g  (* 10)
#   C_content_corr:      % -> mg/g  (* 10)
#   LMA:                 kg/m2 -> g/m2 (* 1000)
#   LDMC:                mg/g       (no conversion)

field_long = field_individual %>%
  select(
    plant_species,
    vegetative_height,
    LMA,
    LDMC,
    N_content_corr,
    C_content_corr
  ) %>%
  mutate(
    vegetative_height = vegetative_height / 10,          # mm -> cm
    N_content_corr = N_content_corr * 10,                # % -> mg/g
    C_content_corr = C_content_corr * 10,                # % -> mg/g
    LMA = LMA * 1000                                     # kg/m2 -> g/m2
  ) %>%
  rename(LNC = N_content_corr, LCC = C_content_corr) %>%
  pivot_longer(
    cols = -plant_species,
    names_to = "trait",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%
  rename(species_TNRS = plant_species) %>%
  mutate(source = "Field")

cat("Field traits standardized:", n_distinct(field_long$trait), "traits\n")

# ==============================================================================
# 5. COMBINE TRY AND FIELD DATA + UNIT CHECK
# ==============================================================================
cat("\n=== Combining TRY and field data ===\n")

# Combine both sources
all_individual = bind_rows(try_long, field_long)

# --- Unit harmonization check ---
cat("\n=== Unit harmonization check ===\n")

# Verify all traits are in the target set
merged_traits = sort(unique(all_individual$trait))
expected_traits = sort(target_units$trait)
unexpected = setdiff(merged_traits, expected_traits)
if (length(unexpected) > 0) {
  cat("WARNING: unexpected traits found:", paste(unexpected, collapse = ", "), "\n")
} else {
  cat("All traits match target unit reference.\n")
}

# Show per-trait ranges by source to verify unit consistency
range_check = all_individual %>%
  group_by(trait, source) %>%
  summarize(
    n = n(),
    min = round(min(value, na.rm = TRUE), 4),
    median = round(median(value, na.rm = TRUE), 4),
    max = round(max(value, na.rm = TRUE), 4),
    .groups = "drop"
  ) %>%
  left_join(target_units, by = "trait") %>%
  arrange(trait, source)

cat("\nTrait value ranges by source (verify units are compatible):\n")
print(range_check, n = 30)

# --- Visual unit comparison: Field vs TRY per species and trait ---
cat("\n=== Generating unit comparison plot ===\n")

# Keep only traits present in both sources
shared_traits = all_individual %>%
  group_by(trait) %>%
  filter(n_distinct(source) == 2) %>%
  ungroup() %>%
  distinct(trait) %>%
  pull(trait)

if (length(shared_traits) > 0) {

  # Compute species-level summaries per source for shared traits
  species_source_summary = all_individual %>%
    filter(trait %in% shared_traits) %>%
    group_by(species_TNRS, trait, source) %>%
    summarize(
      median = median(value, na.rm = TRUE),
      q25 = quantile(value, 0.25, na.rm = TRUE),
      q75 = quantile(value, 0.75, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )

  # Add unit labels for facets
  species_source_summary = species_source_summary %>%
    left_join(target_units, by = "trait") %>%
    mutate(facet_label = paste0(trait, " (", unit, ")"))

  # Abbreviate species names: "Genus species" -> "G. species"
  species_source_summary = species_source_summary %>%
    mutate(species_short = str_replace(species_TNRS, "^(\\w)\\w+", "\\1."))

  # Keep only species with data from both sources (for direct comparison)
  both_sources = species_source_summary %>%
    group_by(species_TNRS, trait) %>%
    filter(n_distinct(source) == 2) %>%
    ungroup()

  if (nrow(both_sources) > 0) {

    # Order species within each trait by overall median (manual reorder per facet)
    both_sources = both_sources %>%
      group_by(species_short, facet_label) %>%
      mutate(overall_median = median(median)) %>%
      group_by(facet_label) %>%
      mutate(species_rank = dense_rank(overall_median)) %>%
      ungroup() %>%
      mutate(species_facet = paste0(species_short, "__", facet_label),
             species_facet = reorder(species_facet, species_rank))

    # Compute mean of species medians per trait and source (for vertical reference lines)
    trait_means = both_sources %>%
      group_by(trait, source, facet_label) %>%
      summarize(mean_value = mean(median, na.rm = TRUE), .groups = "drop")

    p_compare = ggplot(both_sources,
                       aes(x = median, y = species_facet)) +
      geom_vline(data = trait_means, aes(xintercept = mean_value, color = source),
                 linetype = "dashed", linewidth = 0.5, alpha = 0.7) +
      geom_pointrange(aes(xmin = q25, xmax = q75, color = source, fill = n),
                      position = position_dodge(width = 0.6),
                      size = 0.3, fatten = 3, shape = 21, stroke = 0.5) +
      facet_wrap(~ facet_label, scales = "free", ncol = 2) +
      scale_y_discrete(labels = function(x) str_remove(x, "__.*$")) +
      scale_color_manual(values = c("Field" = "#E76F51", "TRY" = "#2A9D8F")) +
      scale_fill_viridis_c(trans = "log10", name = "N obs",
                           option = "viridis", direction = 1) +
      guides(color = guide_legend(order = 1),
             fill = guide_colorbar(order = 2)) +
      labs(x = "Value (median with IQR, dashed lines = dataset means)",
           y = NULL,
           color = "Source",
           title = "Unit harmonization check: Field vs TRY trait distributions",
           subtitle = "Species with data from both sources; point fill = sample size") +
      theme_bw() +
      theme(
        text = element_text(size = 9),
        axis.text.y = element_text(size = 6),
        strip.text = element_text(face = "bold", size = 9),
        legend.position = "top",
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3)
      )

    pdf(here("Calanda_JSDM", "plot", "unit_comparison_field_vs_try.pdf"),
        width = 12, height = max(6, nrow(distinct(both_sources, species_short)) * 0.2))
    print(p_compare)
    dev.off()

    cat("Saved plot/unit_comparison_field_vs_try.pdf\n")

  } else {
    cat("No species with data from both sources — skipping comparison plot.\n")
  }
} else {
  cat("No shared traits between Field and TRY — skipping comparison plot.\n")
}

cat("Combined dataset:", nrow(all_individual), "observations\n")
cat("Species:", n_distinct(all_individual$species_TNRS), "\n")
cat("Traits:", n_distinct(all_individual$trait), "\n")

# Summary by source
cat("\nObservations by source:\n")
all_individual %>%
  group_by(source) %>%
  summarize(
    n_obs = n(),
    n_species = n_distinct(species_TNRS),
    .groups = "drop"
  ) %>%
  print()

# Summary by trait and source
cat("\nObservations by trait and source:\n")
all_individual %>%
  group_by(trait, source) %>%
  summarize(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = source, values_from = n, values_fill = 0) %>%
  mutate(total = TRY + Field) %>%
  arrange(desc(total)) %>%
  print()

# Save combined individual-level data
write_csv(all_individual, here("Calanda_JSDM", "output", "all_traits_individual.csv"))
cat("\nSaved combined individual data to output/all_traits_individual.csv\n")

# ==============================================================================
# 6. DISTRIBUTION PLOTS: RAW VS LOG-TRANSFORMED TRAITS
# ==============================================================================
cat("\n=== Plotting trait distributions (raw vs log-transformed) ===\n")

# Prepare raw and log-transformed data for plotting
dist_data = all_individual %>%
  mutate(log_value = log(value)) %>%
  filter(is.finite(log_value))  # Remove log(0) or log(negative)

# Raw distributions
p_raw = ggplot(dist_data, aes(x = value, fill = source)) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
  facet_wrap(~ trait, scales = "free", ncol = 3) +
  scale_fill_manual(values = c("Field" = "#E76F51", "TRY" = "#2A9D8F")) +
  labs(x = "Value (original scale)", y = "Count",
       title = "Trait distributions (raw values)",
       fill = "Source") +
  theme_bw() +
  theme(legend.position = "top", strip.text = element_text(face = "bold"))

# Log-transformed distributions
p_log = ggplot(dist_data, aes(x = log_value, fill = source)) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
  facet_wrap(~ trait, scales = "free", ncol = 3) +
  scale_fill_manual(values = c("Field" = "#E76F51", "TRY" = "#2A9D8F")) +
  labs(x = "Value (log-transformed)", y = "Count",
       title = "Trait distributions (log-transformed)",
       fill = "Source") +
  theme_bw() +
  theme(legend.position = "top", strip.text = element_text(face = "bold"))

pdf(here("Calanda_JSDM", "plot", "trait_distributions_raw_vs_log.pdf"),
    width = 12, height = 14)
gridExtra::grid.arrange(p_raw, p_log, ncol = 1)
dev.off()

cat("Saved plot/trait_distributions_raw_vs_log.pdf\n")

# ==============================================================================
# 7. CALCULATE SPECIES-LEVEL MEDIANS AND kCV
# ==============================================================================
cat("\n=== Calculating species-level medians and kCV ===\n")

# Species medians: computed from ALL available observations
species_medians = all_individual %>%
  group_by(species_TNRS, trait) %>%
  summarize(
    Median = median(value, na.rm = TRUE),
    n_obs = n(),
    n_sources = n_distinct(source),
    .groups = "drop"
  )

# kCV: computed only for species-trait combos with >= 5 observations
# kCV = CV / (1 + CV), where CV = sd(log(x)) / mean(log(x))
species_kcv = all_individual %>%
  group_by(species_TNRS, trait) %>%
  filter(n() >= 5) %>%
  summarize(
    log_mean = mean(log(value), na.rm = TRUE),
    log_sd = sd(log(value), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    CV = ifelse(log_mean != 0, log_sd / log_mean, NA_real_),
    kCV = ifelse(!is.na(CV), CV / (1 + CV), NA_real_)
  ) %>%
  select(species_TNRS, trait, kCV)

cat(sprintf("Species-trait combinations with >= 5 obs for kCV: %d / %d\n",
            nrow(species_kcv), nrow(species_medians)))

# Combine medians and kCV
species_traits = species_medians %>%
  left_join(species_kcv, by = c("species_TNRS", "trait"))

cat("Species-trait summaries calculated\n")

# Pivot to wide format
species_traits_wide = species_traits %>%
  pivot_wider(
    names_from = trait,
    values_from = c(Median, kCV, n_obs, n_sources),
    names_glue = "{.value}_{trait}"
  )

cat("Species summaries calculated for", nrow(species_traits_wide), "species\n")

# Show summary of kCV values
cat("\nkCV summary per trait:\n")
species_traits %>%
  group_by(trait) %>%
  summarize(
    n_species = n(),
    kCV_mean = round(mean(kCV, na.rm = TRUE), 3),
    kCV_median = round(median(kCV, na.rm = TRUE), 3),
    kCV_min = round(min(kCV, na.rm = TRUE), 3),
    kCV_max = round(max(kCV, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  print()

# ==============================================================================
# 8. SPECIES-LEVEL TRAIT CORRELATIONS
# ==============================================================================
cat("\n=== Species-level trait correlations ===\n")

# Build correlation matrix from species means
cor_data = species_traits_wide %>%
  select(starts_with("Median_")) %>%
  rename_with(~ str_remove(., "^Median_"))

# Need at least some complete pairwise observations
cor_matrix = cor(cor_data, use = "pairwise.complete.obs")

cat("Correlation matrix (pairwise complete observations):\n")
print(round(cor_matrix, 2))

# Save correlation matrix
write_csv(
  as.data.frame(cor_matrix) %>% rownames_to_column("trait"),
  here("Calanda_JSDM", "output", "species_trait_correlations.csv")
)

# Plot correlation matrix
pdf(here("Calanda_JSDM", "plot", "species_trait_correlations.pdf"), width = 8, height = 8)
corrplot(cor_matrix, method = "ellipse", type = "lower",
         tl.col = "black", tl.srt = 45, tl.cex = 0.9,
         addCoef.col = "black", number.cex = 0.7,
         col = colorRampPalette(c("#E76F51", "white", "#2A9D8F"))(200),
         title = "Species-level trait correlations (medians)",
         mar = c(0, 0, 2, 0))
dev.off()

cat("Saved plot/species_trait_correlations.pdf\n")
cat("Saved output/species_trait_correlations.csv\n")

# Build correlation matrix from species kCV values
cat("\n=== Species-level kCV correlations ===\n")

cor_kcv_data = species_traits_wide %>%
  select(starts_with("kCV_")) %>%
  rename_with(~ str_remove(., "^kCV_"))

cor_kcv_matrix = cor(cor_kcv_data, use = "pairwise.complete.obs")

cat("kCV correlation matrix (pairwise complete observations):\n")
print(round(cor_kcv_matrix, 2))

write_csv(
  as.data.frame(cor_kcv_matrix) %>% rownames_to_column("trait"),
  here("Calanda_JSDM", "output", "species_kcv_correlations.csv")
)

pdf(here("Calanda_JSDM", "plot", "species_kcv_correlations.pdf"), width = 8, height = 8)
corrplot(cor_kcv_matrix, method = "ellipse", type = "lower",
         tl.col = "black", tl.srt = 45, tl.cex = 0.9,
         addCoef.col = "black", number.cex = 0.7,
         col = colorRampPalette(c("#E76F51", "white", "#2A9D8F"))(200),
         title = "Species-level trait correlations (kCV)",
         mar = c(0, 0, 2, 0))
dev.off()

cat("Saved plot/species_kcv_correlations.pdf\n")
cat("Saved output/species_kcv_correlations.csv\n")

# ==============================================================================
# 9. MERGE WITH INDICATORS AND DISPERSAL
# ==============================================================================
cat("\n=== Merging with indicators (Nutrients, disturbance) and dispersal ===\n")

# Merge indicator traits (Nutrients + disturbance + dispersal)
traits = species_traits_wide %>%
  left_join(dispersal %>% select(species_TNRS, dispersal), by = "species_TNRS") %>%
  left_join(indicators %>% select(species_TNRS, Nutrients, disturbance),
            by = "species_TNRS")

cat("Merged dataset:", nrow(traits), "species\n")

# Select and order columns for output
traits_ordered = traits %>%
  select(
    species_TNRS,
    # Species means
    Median_vegetative_height,
    Median_LNC,
    Median_LCC,
    Median_LDMC,
    Median_LMA,
    Median_seed_mass,
    Median_LA,
    # kCV values (LCC excluded — too much difference between field and TRY observations)
    kCV_vegetative_height,
    kCV_LNC,
    kCV_LDMC,
    kCV_LMA,
    kCV_seed_mass,
    kCV_LA,
    # Indicator traits
    Nutrients,
    disturbance,
    dispersal,
    # Sample sizes (for QC)
    starts_with("n_obs_"),
    starts_with("n_sources_")
  )

# Save raw (un-imputed) traits
write_csv(traits_ordered, here("Calanda_JSDM", "output", "traits_raw.csv"))
cat("Saved raw traits to output/traits_raw.csv\n")

# ==============================================================================
# 10. IMPUTE SPECIES MEANS
# ==============================================================================
cat("\n=== Imputing missing species-level trait medians ===\n")

# Prepare means data for imputation (continuous means + indicators)
means_for_imputation = traits_ordered %>%
  select(
    species_TNRS,
    Median_vegetative_height,
    Median_LNC,
    Median_LCC,
    Median_LDMC,
    Median_LMA,
    Median_seed_mass,
    Median_LA,
    Nutrients,
    disturbance,
    dispersal
  )

means_vars_to_impute = c(
  "Median_vegetative_height", "Median_LNC", "Median_LCC", "Median_LDMC",
  "Median_LMA", "Median_seed_mass", "Median_LA",
  "Nutrients", "disturbance", "dispersal"
)

# Filter out species with fewer than 3 observed traits (too little information
# for reliable imputation — these would contribute more noise than signal)
n_observed_means = means_for_imputation %>%
  mutate(n_obs_traits = rowSums(!is.na(pick(all_of(means_vars_to_impute))))) %>%
  pull(n_obs_traits)

min_observed = 3
n_excluded = sum(n_observed_means < min_observed)
cat(sprintf("Excluding %d species with fewer than %d observed mean traits\n",
            n_excluded, min_observed))

means_for_imputation = means_for_imputation %>%
  filter(n_observed_means >= min_observed)

cat("\nMissing data in means before imputation:\n")
means_for_imputation %>%
  summarize(across(-species_TNRS, ~ round(100 * mean(is.na(.)), 1))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "percent_missing") %>%
  filter(percent_missing > 0) %>%
  arrange(desc(percent_missing)) %>%
  print(n = 20)

# Run imputation for means (m=10, maxiter=20, num_trees=200 for stability)
means_imputation = impute_functional_traits(
  data = means_for_imputation,
  variables_to_impute = means_vars_to_impute,
  m = 10, maxiter = 20, num_trees = 200, seed = 123
)

# Display imputation performance
cat("\nMeans imputation performance:\n")
print(means_imputation$performance)

if (nrow(means_imputation$performance) > 0) {
  means_imputation$performance %>%
    gt() %>%
    fmt_number(columns = c(r_squared, rmse), decimals = 3) %>%
    tab_header(title = "Means imputation validation (20% holdout)") %>%
    print()
}

# Extract final imputed means: use _final columns, keep species_TNRS
means_imputed = means_imputation$imputed_data %>%
  select(species_TNRS, ends_with("_final")) %>%
  rename_with(~ str_remove(., "_final$"), ends_with("_final"))

cat("\nMissing data in means after imputation:\n")
means_imputed %>%
  summarize(across(-species_TNRS, ~ round(100 * mean(is.na(.)), 1))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "percent_missing") %>%
  filter(percent_missing > 0) %>%
  arrange(desc(percent_missing)) %>%
  print(n = 20)

# Extract imputation flags for means (Original vs Imputed per species per trait)
means_flags = means_imputation$imputed_data %>%
  select(species_TNRS, ends_with("_flag")) %>%
  rename_with(~ str_remove(., "_flag$"), ends_with("_flag"))

write_csv(means_imputed, here("Calanda_JSDM", "output", "traits_medians_imputed.csv"))
write_csv(means_flags, here("Calanda_JSDM", "output", "traits_medians_imputation_flags.csv"))
cat("Saved imputed means to output/traits_medians_imputed.csv\n")
cat("Saved imputation flags to output/traits_medians_imputation_flags.csv\n")

# ==============================================================================
# 11. IMPUTE SPECIES kCV
# ==============================================================================
cat("\n=== Imputing missing species-level kCV values ===\n")

# Prepare kCV data for imputation
# NOTE: LCC (leaf carbon content) is excluded from kCV — it is only used for
# species means. LCC is not included in intraspecific variability analysis.
kcv_for_imputation = traits_ordered %>%
  select(
    species_TNRS,
    kCV_vegetative_height,
    kCV_LNC,
    kCV_LDMC,
    kCV_LMA,
    kCV_seed_mass,
    kCV_LA
  )

kcv_vars_to_impute = c(
  "kCV_vegetative_height", "kCV_LNC", "kCV_LDMC",
  "kCV_LMA", "kCV_seed_mass", "kCV_LA"
)

# Filter out species with fewer than 3 observed kCV traits
n_observed_kcv = kcv_for_imputation %>%
  mutate(n_obs_traits = rowSums(!is.na(pick(all_of(kcv_vars_to_impute))))) %>%
  pull(n_obs_traits)

n_excluded_kcv = sum(n_observed_kcv < min_observed)
cat(sprintf("Excluding %d species with fewer than %d observed kCV traits\n",
            n_excluded_kcv, min_observed))

kcv_for_imputation = kcv_for_imputation %>%
  filter(n_observed_kcv >= min_observed)

cat("\nMissing data in kCV before imputation:\n")
kcv_for_imputation %>%
  summarize(across(-species_TNRS, ~ round(100 * mean(is.na(.)), 1))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "percent_missing") %>%
  filter(percent_missing > 0) %>%
  arrange(desc(percent_missing)) %>%
  print(n = 20)

# Run imputation for kCV (m=10, maxiter=20, num_trees=200 for stability)
kcv_imputation = impute_functional_traits(
  data = kcv_for_imputation,
  variables_to_impute = kcv_vars_to_impute,
  m = 10, maxiter = 20, num_trees = 200, seed = 123
)

# Display imputation performance
cat("\nkCV imputation performance:\n")
print(kcv_imputation$performance)

if (nrow(kcv_imputation$performance) > 0) {
  kcv_imputation$performance %>%
    gt() %>%
    fmt_number(columns = c(r_squared, rmse), decimals = 3) %>%
    tab_header(title = "kCV imputation validation (20% holdout)") %>%
    print()
}

# Extract final imputed kCV: use _final columns, keep species_TNRS
kcv_imputed = kcv_imputation$imputed_data %>%
  select(species_TNRS, ends_with("_final")) %>%
  rename_with(~ str_remove(., "_final$"), ends_with("_final"))

cat("\nMissing data in kCV after imputation:\n")
kcv_imputed %>%
  summarize(across(-species_TNRS, ~ round(100 * mean(is.na(.)), 1))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "percent_missing") %>%
  filter(percent_missing > 0) %>%
  arrange(desc(percent_missing)) %>%
  print(n = 20)

# Extract imputation flags for kCV
kcv_flags = kcv_imputation$imputed_data %>%
  select(species_TNRS, ends_with("_flag")) %>%
  rename_with(~ str_remove(., "_flag$"), ends_with("_flag"))

write_csv(kcv_imputed, here("Calanda_JSDM", "output", "traits_kcv_imputed.csv"))
write_csv(kcv_flags, here("Calanda_JSDM", "output", "traits_kcv_imputation_flags.csv"))
cat("Saved imputed kCV to output/traits_kcv_imputed.csv\n")
cat("Saved kCV imputation flags to output/traits_kcv_imputation_flags.csv\n")

# ==============================================================================
# 12. PCA ON SPECIES MEDIANS (using imputed data)
# ==============================================================================
cat("\n=== PCA on species-level trait medians ===\n")

# abbrev_species is defined in functions_calanda.R

# Build wide matrix of trait means from imputed data (exclude indicators)
pca_means_data = means_imputed %>%
  select(species_TNRS, starts_with("Median_")) %>%
  column_to_rownames("species_TNRS")

# Drop species with any NA (PCA requires complete cases)
n_before = nrow(pca_means_data)
pca_means_data = pca_means_data[complete.cases(pca_means_data), , drop = FALSE]
cat(sprintf("Species with complete mean data: %d / %d\n", nrow(pca_means_data), n_before))

if (nrow(pca_means_data) >= 5 && ncol(pca_means_data) >= 2) {

  pca_means_res = PCA(pca_means_data, scale.unit = TRUE, graph = FALSE)

  # Prepare coordinates
  pca_means_ind = get_pca_ind(pca_means_res)$coord %>%
    as.data.frame() %>%
    rownames_to_column("species_TNRS") %>%
    mutate(species_label = abbrev_species(species_TNRS))

  pca_means_var = get_pca_var(pca_means_res)$coord %>%
    as.data.frame() %>%
    rownames_to_column("trait") %>%
    mutate(trait_label = str_remove(trait, "^Median_"))

  eig_means = get_eigenvalue(pca_means_res)
  pc1_var = round(eig_means[1, "variance.percent"], 1)
  pc2_var = round(eig_means[2, "variance.percent"], 1)
  cat(sprintf("Medians PCA — PC1: %.1f%%, PC2: %.1f%%\n", pc1_var, pc2_var))

  # Arrow scale factor
  arrow_scale = max(abs(pca_means_ind$Dim.1), abs(pca_means_ind$Dim.2)) /
    max(abs(pca_means_var$Dim.1), abs(pca_means_var$Dim.2)) * 0.8

  p_pca_means = ggplot() +
    geom_segment(data = pca_means_var,
                 aes(x = 0, y = 0, xend = Dim.1 * arrow_scale, yend = Dim.2 * arrow_scale),
                 arrow = arrow(length = unit(0.2, "cm")),
                 color = "grey40", linewidth = 0.5) +
    geom_text(data = pca_means_var,
              aes(x = Dim.1 * arrow_scale * 1.1, y = Dim.2 * arrow_scale * 1.1,
                  label = trait_label),
              color = "grey30", size = 3, fontface = "bold") +
    geom_point(data = pca_means_ind,
               aes(x = Dim.1, y = Dim.2),
               size = 2, alpha = 0.7, color = "#2A9D8F") +
    geom_text_repel(data = pca_means_ind,
                    aes(x = Dim.1, y = Dim.2, label = species_label),
                    size = 2, max.overlaps = 20, segment.size = 0.2) +
    labs(x = sprintf("PC1 (%.1f%%)", pc1_var),
         y = sprintf("PC2 (%.1f%%)", pc2_var),
         title = "PCA of species-level trait medians",
         subtitle = sprintf("%d species, %d traits", nrow(pca_means_data), ncol(pca_means_data))) +
    theme_bw() +
    theme(legend.position = "none")

  pdf(here("Calanda_JSDM", "plot", "pca_species_means.pdf"), width = 10, height = 8)
  print(p_pca_means)
  dev.off()

  write_csv(as.data.frame(eig_means) %>% rownames_to_column("PC"),
            here("Calanda_JSDM", "output", "pca_means_eigenvalues.csv"))
  cat("Saved plot/pca_species_means.pdf\n")

} else {
  cat("Not enough complete species for medians PCA — skipping.\n")
}

# ==============================================================================
# 13. PCA ON SPECIES kCV (using imputed data)
# ==============================================================================
cat("\n=== PCA on species-level kCV ===\n")

# Build wide matrix of kCV values from imputed data
pca_kcv_data = kcv_imputed %>%
  column_to_rownames("species_TNRS")

# Drop species with any NA
n_before = nrow(pca_kcv_data)
pca_kcv_data = pca_kcv_data[complete.cases(pca_kcv_data), , drop = FALSE]
cat(sprintf("Species with complete kCV data: %d / %d\n", nrow(pca_kcv_data), n_before))

if (nrow(pca_kcv_data) >= 5 && ncol(pca_kcv_data) >= 2) {

  pca_kcv_res = PCA(pca_kcv_data, scale.unit = TRUE, graph = FALSE)

  # Prepare coordinates
  pca_kcv_ind = get_pca_ind(pca_kcv_res)$coord %>%
    as.data.frame() %>%
    rownames_to_column("species_TNRS") %>%
    mutate(species_label = abbrev_species(species_TNRS))

  pca_kcv_var = get_pca_var(pca_kcv_res)$coord %>%
    as.data.frame() %>%
    rownames_to_column("trait") %>%
    mutate(trait_label = str_remove(trait, "^kCV_"))

  eig_kcv = get_eigenvalue(pca_kcv_res)
  pc1_var = round(eig_kcv[1, "variance.percent"], 1)
  pc2_var = round(eig_kcv[2, "variance.percent"], 1)
  cat(sprintf("kCV PCA — PC1: %.1f%%, PC2: %.1f%%\n", pc1_var, pc2_var))

  # Arrow scale factor
  arrow_scale = max(abs(pca_kcv_ind$Dim.1), abs(pca_kcv_ind$Dim.2)) /
    max(abs(pca_kcv_var$Dim.1), abs(pca_kcv_var$Dim.2)) * 0.8

  p_pca_kcv = ggplot() +
    geom_segment(data = pca_kcv_var,
                 aes(x = 0, y = 0, xend = Dim.1 * arrow_scale, yend = Dim.2 * arrow_scale),
                 arrow = arrow(length = unit(0.2, "cm")),
                 color = "grey40", linewidth = 0.5) +
    geom_text(data = pca_kcv_var,
              aes(x = Dim.1 * arrow_scale * 1.1, y = Dim.2 * arrow_scale * 1.1,
                  label = trait_label),
              color = "grey30", size = 3, fontface = "bold") +
    geom_point(data = pca_kcv_ind,
               aes(x = Dim.1, y = Dim.2),
               size = 2, alpha = 0.7, color = "#E76F51") +
    geom_text_repel(data = pca_kcv_ind,
                    aes(x = Dim.1, y = Dim.2, label = species_label),
                    size = 2, max.overlaps = 20, segment.size = 0.2) +
    labs(x = sprintf("PC1 (%.1f%%)", pc1_var),
         y = sprintf("PC2 (%.1f%%)", pc2_var),
         title = "PCA of species-level kCV (intraspecific variability)",
         subtitle = sprintf("%d species, %d traits", nrow(pca_kcv_data), ncol(pca_kcv_data))) +
    theme_bw() +
    theme(legend.position = "none")

  pdf(here("Calanda_JSDM", "plot", "pca_species_kcv.pdf"), width = 10, height = 8)
  print(p_pca_kcv)
  dev.off()

  write_csv(as.data.frame(eig_kcv) %>% rownames_to_column("PC"),
            here("Calanda_JSDM", "output", "pca_kcv_eigenvalues.csv"))
  cat("Saved plot/pca_species_kcv.pdf\n")

} else {
  cat("Not enough complete species for kCV PCA — skipping.\n")
}

# ==============================================================================
# 14. IMPUTATION DIAGNOSTICS: DISTRIBUTION COMPARISON
# ==============================================================================
cat("\n=== Plotting imputation distribution diagnostics ===\n")

# --- Build long-format data combining raw and imputed values with flags ---

# Means: raw vs imputed
means_raw_long = means_for_imputation %>%
  select(species_TNRS, all_of(means_vars_to_impute)) %>%
  pivot_longer(-species_TNRS, names_to = "variable", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(status = "Original")

means_imputed_only = means_imputation$imputed_data %>%
  select(species_TNRS, all_of(paste0(means_vars_to_impute, "_flag")),
         all_of(paste0(means_vars_to_impute, "_final")))

# Pivot flag and final value separately, then combine
means_imp_long = means_imputed_only %>%
  pivot_longer(
    cols = ends_with("_final"),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(variable = str_remove(variable, "_final$")) %>%
  select(species_TNRS, variable, value)

means_flag_long = means_imputed_only %>%
  pivot_longer(
    cols = ends_with("_flag"),
    names_to = "variable",
    values_to = "status"
  ) %>%
  mutate(variable = str_remove(variable, "_flag$")) %>%
  select(species_TNRS, variable, status)

means_combined = means_imp_long %>%
  left_join(means_flag_long, by = c("species_TNRS", "variable")) %>%
  filter(!is.na(value)) %>%
  mutate(
    trait = str_remove(variable, "^Median_"),
    measure = "Median"
  )

# kCV: raw vs imputed
kcv_raw_long = kcv_for_imputation %>%
  select(species_TNRS, all_of(kcv_vars_to_impute)) %>%
  pivot_longer(-species_TNRS, names_to = "variable", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(status = "Original")

kcv_imputed_only = kcv_imputation$imputed_data %>%
  select(species_TNRS, all_of(paste0(kcv_vars_to_impute, "_flag")),
         all_of(paste0(kcv_vars_to_impute, "_final")))

kcv_imp_long = kcv_imputed_only %>%
  pivot_longer(
    cols = ends_with("_final"),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(variable = str_remove(variable, "_final$")) %>%
  select(species_TNRS, variable, value)

kcv_flag_long = kcv_imputed_only %>%
  pivot_longer(
    cols = ends_with("_flag"),
    names_to = "variable",
    values_to = "status"
  ) %>%
  mutate(variable = str_remove(variable, "_flag$")) %>%
  select(species_TNRS, variable, status)

kcv_combined = kcv_imp_long %>%
  left_join(kcv_flag_long, by = c("species_TNRS", "variable")) %>%
  filter(!is.na(value)) %>%
  mutate(
    trait = str_remove(variable, "^kCV_"),
    measure = "kCV"
  )

# Combine means and kCV for plotting
all_imp_data = bind_rows(means_combined, kcv_combined) %>%
  mutate(
    # Order: means first, then kCV; same trait side by side
    facet_label = paste0(measure, " — ", trait),
    measure = factor(measure, levels = c("Median", "kCV")),
    trait = factor(trait, levels = c("vegetative_height", "LNC", "LCC", "LDMC",
                                     "LMA", "seed_mass", "LA",
                                     "Nutrients", "disturbance", "dispersal"))
  )

# Calculate summary statistics per facet and status
imp_stats = all_imp_data %>%
  group_by(trait, measure, status) %>%
  summarize(
    n = n(),
    mean_val = mean(value, na.rm = TRUE),
    median_val = median(value, na.rm = TRUE),
    min_val = min(value, na.rm = TRUE),
    max_val = max(value, na.rm = TRUE),
    .groups = "drop"
  )

# Build annotation labels
imp_labels = imp_stats %>%
  mutate(
    label = sprintf("%s (n=%d)\nmean=%.2f  med=%.2f\nmin=%.2f  max=%.2f",
                    status, n, mean_val, median_val, min_val, max_val)
  )

# Separate plots for means and kCV so each has truly free x scales

# Shared theme
imp_theme = theme_bw() +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold", size = 8),
    axis.text = element_text(size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# --- Means plot ---
means_plot_data = all_imp_data %>% filter(measure == "Median")
means_plot_stats = imp_stats %>% filter(measure == "Median")
means_plot_labels = imp_labels %>% filter(measure == "Median")

p_imp_means = ggplot(means_plot_data, aes(x = value, fill = status)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  geom_vline(data = means_plot_stats, aes(xintercept = mean_val, color = status),
             linetype = "dashed", linewidth = 0.6) +
  geom_vline(data = means_plot_stats, aes(xintercept = median_val, color = status),
             linetype = "solid", linewidth = 0.6) +
  facet_wrap(~ trait, scales = "free", nrow = 1) +
  scale_fill_manual(values = c("Original" = "#2A9D8F", "Imputed" = "#E76F51")) +
  scale_color_manual(values = c("Original" = "#2A9D8F", "Imputed" = "#E76F51")) +
  geom_text(data = means_plot_labels %>% filter(status == "Original"),
            aes(x = Inf, y = Inf, label = label),
            hjust = 1.05, vjust = 1.1, size = 1.8, color = "#2A9D8F",
            inherit.aes = FALSE) +
  geom_text(data = means_plot_labels %>% filter(status == "Imputed"),
            aes(x = Inf, y = Inf, label = label),
            hjust = 1.05, vjust = 4.0, size = 1.8, color = "#E76F51",
            inherit.aes = FALSE) +
  labs(x = "Value", y = "Count", fill = "Status", color = "Status",
       title = "Species medians: Original vs Imputed distributions",
       subtitle = "Dashed line = mean, solid line = median") +
  imp_theme

# --- kCV plot ---
kcv_plot_data = all_imp_data %>% filter(measure == "kCV")
kcv_plot_stats = imp_stats %>% filter(measure == "kCV")
kcv_plot_labels = imp_labels %>% filter(measure == "kCV")

p_imp_kcv = ggplot(kcv_plot_data, aes(x = value, fill = status)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  geom_vline(data = kcv_plot_stats, aes(xintercept = mean_val, color = status),
             linetype = "dashed", linewidth = 0.6) +
  geom_vline(data = kcv_plot_stats, aes(xintercept = median_val, color = status),
             linetype = "solid", linewidth = 0.6) +
  facet_wrap(~ trait, scales = "free", nrow = 1) +
  scale_fill_manual(values = c("Original" = "#2A9D8F", "Imputed" = "#E76F51")) +
  scale_color_manual(values = c("Original" = "#2A9D8F", "Imputed" = "#E76F51")) +
  geom_text(data = kcv_plot_labels %>% filter(status == "Original"),
            aes(x = Inf, y = Inf, label = label),
            hjust = 1.05, vjust = 1.1, size = 1.8, color = "#2A9D8F",
            inherit.aes = FALSE) +
  geom_text(data = kcv_plot_labels %>% filter(status == "Imputed"),
            aes(x = Inf, y = Inf, label = label),
            hjust = 1.05, vjust = 4.0, size = 1.8, color = "#E76F51",
            inherit.aes = FALSE) +
  labs(x = "Value", y = "Count", fill = "Status", color = "Status",
       title = "Species kCV: Original vs Imputed distributions",
       subtitle = "Dashed line = mean, solid line = median") +
  imp_theme

# Combine and save
pdf(here("Calanda_JSDM", "plot", "imputation_distributions.pdf"),
    width = 24, height = 10)
gridExtra::grid.arrange(p_imp_means, p_imp_kcv, ncol = 1)
dev.off()

cat("Saved plot/imputation_distributions.pdf\n")

# Print summary table to console
cat("\nImputation summary statistics:\n")
imp_stats %>%
  mutate(across(c(mean_val, median_val, min_val, max_val), ~ round(., 3))) %>%
  arrange(measure, trait, status) %>%
  print(n = 40)

# ==============================================================================
# 15. MERGE SUMMARY
# ==============================================================================
cat("\n=== Merge Summary ===\n")
cat("TRY individual observations:", nrow(try_individual), "\n")
cat("Field individual observations:", nrow(field_individual), "\n")
cat("Combined individual observations:", nrow(all_individual), "\n")
cat("Species-trait combos with kCV (>= 5 obs):", nrow(species_kcv), "\n")
cat("Species in final traits:", nrow(traits_ordered), "\n")

cat("\nMissing data in traits_raw.csv:\n")
traits_ordered %>%
  summarize(across(everything(), ~ round(100 * mean(is.na(.)), 1))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "percent_missing") %>%
  filter(percent_missing > 0) %>%
  arrange(desc(percent_missing)) %>%
  print(n = 30)

cat("\n=== Merge completed successfully ===\n")
