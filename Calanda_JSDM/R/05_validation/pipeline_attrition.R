# ==============================================================================
# Script: pipeline_attrition.R
# Purpose: Track species and site attrition from raw data through all analysis
#          stages. Reconstructs filtering logic READ-ONLY from raw data and
#          saved outputs — NEVER runs data prep scripts or overwrites any files.
#
# Outputs:
#   - output/species_attrition.csv (one row per species, columns for each filter)
#   - output/site_attrition.csv (one row per site, columns for each filter)
#   - Console: detailed attrition tables
#
# Usage: source("R/05_validation/pipeline_attrition.R")
# ==============================================================================

library(tidyverse)
library(here)
library(stringi)

# ==============================================================================
# HELPERS
# ==============================================================================

print_row = function(stage, n, n_removed, reason, n_start) {
  pct = sprintf("%.1f%%", 100 * n / n_start)
  cat(sprintf("  %-55s  %5d  %5d removed  %-35s  %6s remaining\n",
              stage, n, n_removed, paste0("(", reason, ")"), pct))
}

divider = function(title) {
  cat(sprintf("\n%s\n  %s\n%s\n",
              strrep("=", 80), title, strrep("=", 80)))
}

# ==============================================================================
# AUTO-DETECT param_tag
# ==============================================================================

fold_files = list.files(
  here("Calanda_JSDM", "output", "results", "cv"),
  pattern = "^fold_.*\\.rds$"
)
stopifnot(length(fold_files) > 0)
param_tag = regmatches(fold_files[1],
                       regexpr("a[0-9.]+_l[0-9.]+_le[0-9.]+(?=\\.)",
                               fold_files[1], perl = TRUE))
cat(sprintf("Detected param_tag: %s\n", param_tag))

# ==============================================================================
# LOAD DATA (read-only — no files are modified)
# ==============================================================================
cat("\n=== Loading data (read-only) ===\n")

# Raw vegetation
raw_veg = read_csv(
  here("Calanda_JSDM", "data", "vegetation",
       "2024_CAPHE_SpeDis_CleanData_20240214.csv"),
  show_col_types = FALSE
)

# Environmental data (already saved by script 04)
veg_clim = read_csv(here("Calanda_JSDM", "output", "veg_clim.csv"),
                     show_col_types = FALSE) %>%
  select(-any_of("...1"))

# Final JSDM data
data_calanda = readRDS(
  here("Calanda_JSDM", "output", "data_calanda_jsdm_2026-03-06.rds")
)
X = data_calanda$X
Y = data_calanda$Y

# Species and site CV
species_cv = read_csv(
  here("Calanda_JSDM", "output", "results",
       paste0("exp6_species_cv_", param_tag, ".csv")),
  show_col_types = FALSE
)
site_cv = read_csv(
  here("Calanda_JSDM", "output", "results",
       paste0("exp6_sites_cv_", param_tag, ".csv")),
  show_col_types = FALSE
)

# VP summaries
vp_si = read_csv(
  here("Calanda_JSDM", "output", "results",
       paste0("vp_sites_summary_", param_tag, ".csv")),
  show_col_types = FALSE
)

# Trait data
traits_medians = read_csv(
  here("Calanda_JSDM", "output", "traits_medians_imputed.csv"),
  show_col_types = FALSE
)
traits_kcv = read_csv(
  here("Calanda_JSDM", "output", "traits_kcv_imputed.csv"),
  show_col_types = FALSE
)
sp_distinct = read_csv(
  here("Calanda_JSDM", "output", "species_functional_distinctiveness.csv"),
  show_col_types = FALSE
)
comm_traits = read_csv(
  here("Calanda_JSDM", "output", "community_traits_unweighted_imputed.csv"),
  show_col_types = FALSE
)

cat("  All data loaded.\n")

# ==============================================================================
# CONSTANTS
# ==============================================================================

woody_species = c(
  "Fagus sylvatica", "Picea abies", "Sorbus aucuparia", "Acer pseudoplatanus",
  "Larix decidua", "Fraxinus excelsior", "Pinus sylvestris", "Corylus avellana",
  "Cornus sanguinea", "Ligustrum vulgare", "Lonicera xylosteum",
  "Rubus idaeus", "Rubus fruticosus", "Rubus caesius"
)

core_traits = c("vegetative_height", "LNC", "LA", "seed_mass", "LDMC")
sp_median_cols  = paste0("Median_", core_traits)
sp_kcv_cols     = paste0("kCV_", core_traits)
sp_distinct_cols = paste0("distinct_", core_traits)
sp_trait_cols   = c(sp_median_cols, sp_kcv_cols, sp_distinct_cols)


# ##############################################################################
#
#                  SPECIES ATTRITION — GRANULAR DATA PREP
#
# ##############################################################################
divider("SPECIES ATTRITION")
cat(sprintf("  %-55s  %5s  %13s  %-35s  %6s\n",
            "Stage", "N", "", "Reason", "Cumul."))
cat(sprintf("  %s\n", strrep("-", 130)))

# ---- Step 0: All raw species (including genus-only) ----
all_raw_species = raw_veg$taxon_global %>%
  stri_trans_general("Latin-ASCII") %>%
  unique() %>%
  sort()

# ---- Step 1a: Remove genus-only names ----
raw_species_binomial = all_raw_species[!is.na(word(all_raw_species, 2))]
removed_genus_only = setdiff(all_raw_species, raw_species_binomial)
n_raw_sp = length(raw_species_binomial)

print_row("1a. Raw veg: binomial species", n_raw_sp,
          length(removed_genus_only), "genus-only names removed", n_raw_sp)

# ---- Step 1b: Reconstruct 5% prevalence filter ----
# Replay script 04 logic in-memory (read-only, no files written)
veg_clean = raw_veg %>%
  mutate(plot_id_releve = paste0(plot_id, releve_id)) %>%
  mutate(species = stri_trans_general(taxon_global, "Latin-ASCII")) %>%
  filter(plot_id_releve != "dupl_108.OID2979333") %>%
  filter(!is.na(word(species, 2)))

veg_comm = veg_clean %>%
  select(plot_id_releve, species_cover, species) %>%
  rename(cover = species_cover) %>%
  distinct() %>%
  filter(!is.na(species)) %>%
  mutate(cover = ifelse(is.na(cover), 0.001, cover)) %>%
  group_by(plot_id_releve, species) %>%
  summarize(cover = sum(cover), .groups = "drop") %>%
  group_by(plot_id_releve) %>%
  mutate(total_cover = sum(cover), rel_cover = cover / total_cover) %>%
  select(-total_cover, -cover) %>%
  ungroup()

# PA matrix for prevalence calculation
veg_pa = veg_comm %>%
  mutate(pa = ifelse(rel_cover > 0, 1, 0)) %>%
  select(plot_id_releve, species, pa) %>%
  pivot_wider(names_from = species, values_from = pa, values_fill = 0) %>%
  column_to_rownames(var = "plot_id_releve")

# Species failing 5% prevalence
n_plots_for_prev = nrow(veg_pa)
sp_occurrences = colSums(veg_pa)
rare_species = names(sp_occurrences[sp_occurrences < n_plots_for_prev * 0.05])
sp_after_rare = setdiff(colnames(veg_pa), rare_species)

print_row("1b. After 5% prevalence filter", length(sp_after_rare),
          length(rare_species),
          sprintf("< 5%% of %d plots", n_plots_for_prev), n_raw_sp)

# ---- Step 1c: 70% coverage filter on plots ----
# This removes plots, not species directly. But some species may lose all plots.
veg_abund = veg_comm %>%
  filter(!species %in% rare_species) %>%
  group_by(plot_id_releve) %>%
  mutate(tot_rel_cover = sum(rel_cover)) %>%
  filter(tot_rel_cover >= 0.70) %>%
  select(-tot_rel_cover) %>%
  pivot_wider(names_from = species, values_from = rel_cover, values_fill = 0)

plots_after_coverage = unique(veg_abund$plot_id_releve)
sp_after_coverage = setdiff(colnames(veg_abund), "plot_id_releve")
sp_lost_coverage = setdiff(sp_after_rare, sp_after_coverage)

print_row("1c. After 70% coverage filter on plots", length(sp_after_coverage),
          length(sp_lost_coverage),
          sprintf("%d plots removed", n_plots_for_prev - length(plots_after_coverage)),
          n_raw_sp)

# ---- Step 1d: Intersection with environmental data ----
plots_with_env = veg_clim$plot_id_releve
plots_after_env = intersect(plots_after_coverage, plots_with_env)
sp_lost_env_intersection = character(0)  # species don't change here, only plots

# After env intersection, rebuild Y-like matrix to check for species dropout
veg_abund_env = veg_abund %>%
  filter(plot_id_releve %in% plots_after_env) %>%
  column_to_rownames("plot_id_releve")

# Check for species that became all-zero after plot removal
sp_sums_after_env = colSums(veg_abund_env > 0)
sp_lost_env = names(sp_sums_after_env[sp_sums_after_env == 0])
sp_after_env = setdiff(sp_after_coverage, sp_lost_env)

print_row("1d. After env data intersection", length(sp_after_env),
          length(sp_lost_env),
          sprintf("%d plots removed, species dropped if 0 occ",
                  length(plots_after_coverage) - length(plots_after_env)),
          n_raw_sp)

# ---- Step 1e: Final Y after nutrient/disturbance join + na.omit ----
# The final Y may lose more sites from the CWM nutrient/disturbance join + na.omit
sp_in_Y = colnames(Y)
sp_lost_final_prep = setdiff(sp_after_env, sp_in_Y)

# Also catch anything we missed in reconstruction
sp_missed = setdiff(sp_in_Y, sp_after_env)
if (length(sp_missed) > 0) {
  cat(sprintf("  NOTE: %d species in final Y not captured by reconstruction: %s\n",
              length(sp_missed), paste(head(sp_missed, 5), collapse = ", ")))
}

n_after_prep = length(sp_in_Y)
print_row("1e. Final Y (after CWM join + na.omit)", n_after_prep,
          length(sp_lost_final_prep),
          "sites lost to nutrient/disturbance NA", n_raw_sp)

# ---- Step 2: AUC >= 0.7 filter ----
sp_auc_pass = species_cv %>%
  filter(test_auc >= 0.7, !is.na(test_auc)) %>%
  pull(species)
sp_auc_fail = setdiff(sp_in_Y, sp_auc_pass)
n_after_auc = length(sp_auc_pass)

print_row("2. After AUC >= 0.7 filter", n_after_auc,
          length(sp_auc_fail), "AUC < 0.7", n_raw_sp)

# ---- Step 3: Script 09 — no additional species filter ----
print_row("3. Script 09 (env gradient)", n_after_auc, 0, "none", n_raw_sp)

# ---- Step 4a: Script 10 — remove woody species ----
actually_removed_woody = intersect(sp_auc_pass, woody_species)
sp_after_woody = setdiff(sp_auc_pass, woody_species)
n_after_woody = length(sp_after_woody)

print_row("4a. Script 10: remove woody", n_after_woody,
          length(actually_removed_woody), "woody exclusion", n_raw_sp)

# ---- Step 4b: Script 10 — trait name merge ----
trait_species = traits_medians$species_TNRS
sp_after_trait_merge = intersect(sp_after_woody, trait_species)
sp_trait_mismatch = setdiff(sp_after_woody, trait_species)
n_after_trait_merge = length(sp_after_trait_merge)

# ---- Step 4b: Script 10 — trait merge + drop_na ----
sp_med_df = traits_medians %>%
  rename(species_name = species_TNRS) %>%
  filter(species_name %in% sp_after_trait_merge) %>%
  select(species_name, all_of(sp_median_cols))

sp_kcv_df = traits_kcv %>%
  rename(species_name = species_TNRS) %>%
  select(species_name, any_of(sp_kcv_cols))

sp_dist_df = sp_distinct %>%
  rename(species_name = species) %>%
  select(species_name, any_of(sp_distinct_cols))

sp_trait_check = sp_med_df %>%
  left_join(sp_kcv_df, by = "species_name") %>%
  left_join(sp_dist_df, by = "species_name")

available_trait_cols = intersect(sp_trait_cols, names(sp_trait_check))
sp_trait_complete = sp_trait_check %>% drop_na(all_of(available_trait_cols))
sp_after_dropna = sp_trait_complete$species_name
sp_trait_na = setdiff(sp_after_trait_merge, sp_after_dropna)
n_after_dropna = length(sp_after_dropna)

n_all_trait_missing = length(sp_trait_mismatch) + length(sp_trait_na)
print_row("4b. Script 10: missing trait data (final)", n_after_dropna,
          n_all_trait_missing, "no trait data or incomplete traits", n_raw_sp)

cat(sprintf("\n  SPECIES TOTAL: %d raw -> %d final functional (%.1f%% retained)\n",
            n_raw_sp, n_after_dropna, 100 * n_after_dropna / n_raw_sp))

# ==============================================================================
# SPECIES REMOVED — DETAILED LISTS
# ==============================================================================
divider("SPECIES REMOVED — DETAILED LISTS")

cat(sprintf("\n  -- Genus-only names (%d) --\n", length(removed_genus_only)))
cat(paste0("    ", sort(removed_genus_only)), sep = "\n")

cat(sprintf("\n  -- Rare species: < 5%% prevalence (%d) --\n", length(rare_species)))
cat(paste0("    ", sort(rare_species)), sep = "\n")

cat(sprintf("\n  -- Lost to 70%% coverage plot filter (%d) --\n", length(sp_lost_coverage)))
if (length(sp_lost_coverage) > 0) { cat(paste0("    ", sort(sp_lost_coverage)), sep = "\n")
} else cat("    (none)\n")

cat(sprintf("\n  -- Lost to env data intersection (%d) --\n", length(sp_lost_env)))
if (length(sp_lost_env) > 0) { cat(paste0("    ", sort(sp_lost_env)), sep = "\n")
} else cat("    (none)\n")

cat(sprintf("\n  -- Lost to CWM nutrient/disturbance join + na.omit (%d) --\n",
            length(sp_lost_final_prep)))
if (length(sp_lost_final_prep) > 0) { cat(paste0("    ", sort(sp_lost_final_prep)), sep = "\n")
} else cat("    (none)\n")

cat(sprintf("\n  -- AUC < 0.7 (%d) --\n", length(sp_auc_fail)))
cat(paste0("    ", sort(sp_auc_fail)), sep = "\n")

cat(sprintf("\n  -- Woody species removed (%d) --\n", length(actually_removed_woody)))
cat(paste0("    ", sort(actually_removed_woody)), sep = "\n")

sp_all_trait_missing = c(sp_trait_mismatch, sp_trait_na)
cat(sprintf("\n  -- Missing trait data (%d) --\n", length(sp_all_trait_missing)))
if (length(sp_all_trait_missing) > 0) {
  for (sp in sort(sp_all_trait_missing)) {
    if (sp %in% sp_trait_mismatch) {
      cat(sprintf("    %-40s  no trait data available\n", sp))
    } else {
      row = sp_trait_check %>% filter(species_name == sp)
      missing_cols = available_trait_cols[is.na(row[1, available_trait_cols])]
      cat(sprintf("    %-40s  missing: %s\n", sp, paste(missing_cols, collapse = ", ")))
    }
  }
} else {
  cat("    (none)\n")
}


# ##############################################################################
#
#                          SITE ATTRITION
#
# ##############################################################################
divider("SITE ATTRITION")
cat(sprintf("  %-55s  %5s  %13s  %-35s  %6s\n",
            "Stage", "N", "", "Reason", "Cumul."))
cat(sprintf("  %s\n", strrep("-", 130)))

# ---- Step 1: Raw unique plots ----
raw_sites = raw_veg %>%
  mutate(plot_id_releve = paste0(plot_id, releve_id)) %>%
  pull(plot_id_releve) %>%
  unique()
n_raw_sites = length(raw_sites)

print_row("1. Raw veg: unique plots", n_raw_sites, 0,
          "starting count", n_raw_sites)

# ---- Step 2: After genus-only removal ----
sites_after_genus = unique(veg_clean$plot_id_releve)
n_removed_genus_sites = n_raw_sites - length(sites_after_genus)

print_row("2. After removing genus-only species", length(sites_after_genus),
          n_removed_genus_sites, "plots with only genus-level taxa", n_raw_sites)

# ---- Step 3: After 70% coverage filter ----
n_removed_coverage = length(sites_after_genus) - length(plots_after_coverage)

print_row("3. After 70% coverage filter", length(plots_after_coverage),
          n_removed_coverage, "< 70% of abundance from non-rare species", n_raw_sites)

# ---- Step 4: After env data intersection ----
n_removed_env = length(plots_after_coverage) - length(plots_after_env)

print_row("4. After env data intersection", length(plots_after_env),
          n_removed_env, "missing climate/topo/ET data", n_raw_sites)

# ---- Step 5: Final Y (after CWM join + na.omit) ----
sites_in_Y = rownames(Y)
n_removed_final_prep = length(plots_after_env) - length(sites_in_Y)

print_row("5. Final Y (after CWM join + na.omit)", length(sites_in_Y),
          n_removed_final_prep, "nutrient/disturbance NA", n_raw_sites)

# ---- Step 6: After logloss <= log(2) ----
sites_ll_pass = site_cv %>%
  filter(logloss <= log(2), !is.na(logloss)) %>%
  pull(site)
n_removed_ll = length(sites_in_Y) - length(sites_ll_pass)

print_row("6. After logloss <= log(2) filter", length(sites_ll_pass),
          n_removed_ll, "high logloss", n_raw_sites)

# ---- Step 7: After ci_width < 0.1 ----
site_names_all = rownames(X)
site_ci = vp_si %>%
  filter(component %in% c("env", "spa", "codist")) %>%
  group_by(site_idx) %>%
  summarise(max_ci = max(ci_width, na.rm = TRUE), .groups = "drop") %>%
  mutate(site_name = site_names_all[site_idx])

sites_ci_pass = site_ci %>%
  filter(site_name %in% sites_ll_pass, max_ci < 0.1) %>%
  pull(site_name)
n_removed_ci = length(sites_ll_pass) - length(sites_ci_pass)

print_row("7. After ci_width < 0.1 filter", length(sites_ci_pass),
          n_removed_ci, "unstable VP estimates", n_raw_sites)

# ---- Step 8: Script 10 community trait merge ----
comm_sites = comm_traits$plot_id_releve
sites_comm_pass = intersect(sites_ci_pass, comm_sites)
n_removed_comm = length(sites_ci_pass) - length(sites_comm_pass)

print_row("8. Script 10: community trait merge (final)", length(sites_comm_pass),
          n_removed_comm, "woody-dominated plots excluded", n_raw_sites)

cat(sprintf("\n  SITE TOTAL: %d raw -> %d final community (%.1f%% retained)\n",
            n_raw_sites, length(sites_comm_pass),
            100 * length(sites_comm_pass) / n_raw_sites))


# ##############################################################################
#
#                    BUILD & SAVE SPECIES ATTRITION CSV
#
# ##############################################################################
divider("SAVING SPECIES ATTRITION CSV")

# Pre-compute missing trait reason for each species (vectorized lookup)
trait_na_reason = sapply(sp_trait_na, function(sp) {
  row = sp_trait_check %>% filter(species_name == sp)
  if (nrow(row) > 0) {
    missing = available_trait_cols[is.na(row[1, available_trait_cols])]
    paste("missing:", paste(missing, collapse = ", "))
  } else {
    "missing trait data"
  }
})
names(trait_na_reason) = sp_trait_na

# One row per species from the raw binomial list
sp_attrition = tibble(species = sort(raw_species_binomial)) %>%
  mutate(
    removed_at = case_when(
      species %in% rare_species           ~ "1b_prevalence_filter",
      species %in% sp_lost_coverage       ~ "1c_coverage_filter",
      species %in% sp_lost_env            ~ "1d_env_intersection",
      species %in% sp_lost_final_prep     ~ "1e_cwm_na_omit",
      species %in% sp_auc_fail            ~ "2_auc_filter",
      species %in% actually_removed_woody ~ "4a_woody_exclusion",
      species %in% sp_trait_mismatch      ~ "4b_missing_trait_data",
      species %in% sp_trait_na            ~ "4b_missing_trait_data",
      species %in% sp_after_dropna        ~ "retained",
      TRUE                                ~ "retained_pre_script10"
    ),
    reason = case_when(
      removed_at == "1b_prevalence_filter" ~
        sprintf("< 5%% prevalence (%d/%d plots = %.1f%%)",
                sp_occurrences[species], n_plots_for_prev,
                100 * sp_occurrences[species] / n_plots_for_prev),
      removed_at == "1c_coverage_filter"     ~ "all plots lost to 70% coverage filter",
      removed_at == "1d_env_intersection"    ~ "all plots lost when intersecting with env data",
      removed_at == "1e_cwm_na_omit"         ~ "sites lost during CWM nutrient/disturbance join + na.omit",
      removed_at == "2_auc_filter" ~
        sprintf("AUC = %.3f (< 0.7 threshold)",
                species_cv$test_auc[match(species, species_cv$species)]),
      removed_at == "4a_woody_exclusion"     ~ "woody species excluded from functional analysis",
      removed_at == "4b_missing_trait_data" & species %in% sp_trait_mismatch ~ "no trait data available",
      removed_at == "4b_missing_trait_data" & species %in% sp_trait_na ~ trait_na_reason[species],
      removed_at == "retained"               ~ "retained through all filters",
      removed_at == "retained_pre_script10"  ~ "retained in env gradient (script 09); not in functional (script 10)",
      TRUE ~ ""
    ),
    in_Y_matrix = species %in% sp_in_Y,
    in_env_gradient_09 = species %in% sp_auc_pass,
    in_functional_10 = species %in% sp_after_dropna,
    prevalence_in_raw = sp_occurrences[species] / n_plots_for_prev,
    auc = species_cv$test_auc[match(species, species_cv$species)]
  ) %>%
  mutate(
    prevalence_in_raw = round(prevalence_in_raw, 4),
    auc = round(auc, 4)
  )

write_csv(sp_attrition, here("Calanda_JSDM", "output", "species_attrition.csv"))
cat(sprintf("  Saved output/species_attrition.csv (%d species)\n", nrow(sp_attrition)))

# Quick summary of the CSV
cat("\n  Species by removal stage:\n")
sp_attrition %>%
  count(removed_at) %>%
  arrange(removed_at) %>%
  mutate(label = sprintf("    %-35s  %d", removed_at, n)) %>%
  pull(label) %>%
  cat(sep = "\n")

# ##############################################################################
#
#                    BUILD & SAVE SITE ATTRITION CSV
#
# ##############################################################################
cat("\n")
divider("SAVING SITE ATTRITION CSV")

# One row per site from the raw data
site_attrition = tibble(site = sort(raw_sites)) %>%
  mutate(
    in_after_genus = site %in% sites_after_genus,
    in_after_coverage = site %in% plots_after_coverage,
    in_after_env = site %in% plots_after_env,
    in_Y_matrix = site %in% sites_in_Y,
    in_after_logloss = site %in% sites_ll_pass,
    in_after_ci_width = site %in% sites_ci_pass,
    in_community_analysis = site %in% sites_comm_pass,
    removed_at = case_when(
      !in_after_genus    ~ "2_genus_only_plots",
      !in_after_coverage ~ "3_coverage_filter",
      !in_after_env      ~ "4_env_intersection",
      !in_Y_matrix       ~ "5_cwm_na_omit",
      !in_after_logloss  ~ "6_logloss_filter",
      !in_after_ci_width ~ "7_ci_width_filter",
      !in_community_analysis ~ "8_community_trait_merge",
      TRUE               ~ "retained"
    ),
    logloss = site_cv$logloss[match(site, site_cv$site)],
    max_ci_width = site_ci$max_ci[match(site, site_ci$site_name)]
  ) %>%
  mutate(
    logloss = round(logloss, 4),
    max_ci_width = round(max_ci_width, 4)
  )

write_csv(site_attrition, here("Calanda_JSDM", "output", "site_attrition.csv"))
cat(sprintf("  Saved output/site_attrition.csv (%d sites)\n", nrow(site_attrition)))

cat("\n  Sites by removal stage:\n")
site_attrition %>%
  count(removed_at) %>%
  arrange(removed_at) %>%
  mutate(label = sprintf("    %-35s  %d", removed_at, n)) %>%
  pull(label) %>%
  cat(sep = "\n")

# ##############################################################################
#
#                         SUMMARY TABLES
#
# ##############################################################################
divider("SUMMARY TABLE — SPECIES")

sp_stages = tibble(
  Stage = c(
    "1a. Raw binomial species",
    "1b. After 5% prevalence filter",
    "1c. After 70% coverage filter",
    "1d. After env intersection",
    "1e. Final Y (CWM join + na.omit)",
    "2.  After AUC >= 0.7",
    "3.  Script 09 (env gradient)",
    "4a. Script 10: remove woody",
    "4b. Script 10: missing trait data (final)"
  ),
  N = c(n_raw_sp, length(sp_after_rare), length(sp_after_coverage),
        length(sp_after_env), n_after_prep, n_after_auc, n_after_auc,
        n_after_woody, n_after_dropna),
  Removed = c(length(removed_genus_only), length(rare_species),
              length(sp_lost_coverage), length(sp_lost_env),
              length(sp_lost_final_prep), length(sp_auc_fail), 0,
              length(actually_removed_woody), n_all_trait_missing),
  Pct_of_raw = sprintf("%.1f%%", 100 * N / n_raw_sp)
)
print(sp_stages, n = Inf)

divider("SUMMARY TABLE — SITES")

si_stages = tibble(
  Stage = c(
    "1.  Raw unique plots",
    "2.  After genus-only removal",
    "3.  After 70% coverage filter",
    "4.  After env intersection",
    "5.  Final Y (CWM join + na.omit)",
    "6.  After logloss <= log(2)",
    "7.  After ci_width < 0.1",
    "8.  Script 10 community merge (final)"
  ),
  N = c(n_raw_sites, length(sites_after_genus), length(plots_after_coverage),
        length(plots_after_env), length(sites_in_Y), length(sites_ll_pass),
        length(sites_ci_pass), length(sites_comm_pass)),
  Removed = c(0, n_removed_genus_sites, n_removed_coverage, n_removed_env,
              n_removed_final_prep, n_removed_ll, n_removed_ci, n_removed_comm),
  Pct_of_raw = sprintf("%.1f%%", 100 * N / n_raw_sites)
)
print(si_stages, n = Inf)

# ##############################################################################
#
#                    SAVE SUMMARY CSV
#
# ##############################################################################
divider("SAVING SUMMARY CSV")

summary_attrition = bind_rows(
  sp_stages %>% mutate(type = "species"),
  si_stages %>% mutate(type = "sites")
) %>%
  select(type, Stage, N, Removed, Pct_of_raw)

write_csv(summary_attrition, here("Calanda_JSDM", "output", "attrition_summary.csv"))
cat(sprintf("  Saved output/attrition_summary.csv (%d rows)\n", nrow(summary_attrition)))

cat("\n\nDone.\n")
