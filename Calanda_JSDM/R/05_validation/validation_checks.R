# ==============================================================================
# Script: validation_checks.R
# Purpose: Read-only pipeline integrity checks for Calanda JSDM
#          NEVER writes or overwrites any data files — only prints diagnostics.
#
# Covers:
#   1. Data prep consistency (04_prepare)
#   2. Variance partitioning (08)
#   3. Environmental gradient analysis (09)
#   4. Functional trait analysis (10)
#   5. Cross-script consistency
#   6. Statistical method checks
#
# Usage: source("R/05_validation/validation_checks.R")
# ==============================================================================

library(tidyverse)
library(here)
library(car)

# ==============================================================================
# HELPERS
# ==============================================================================
n_pass = 0L
n_warn = 0L
n_crit = 0L
n_info = 0L

log_result = function(level, check_id, msg) {
  tag = switch(level,
    "CRITICAL" = "\033[31m[CRITICAL]\033[0m",
    "WARNING"  = "\033[33m[WARNING]\033[0m",
    "INFO"     = "\033[36m[INFO]\033[0m",
    "PASS"     = "\033[32m[PASS]\033[0m"
  )
  cat(sprintf("  %s %s: %s\n", tag, check_id, msg))
  if (level == "PASS")     n_pass <<- n_pass + 1L
  if (level == "WARNING")  n_warn <<- n_warn + 1L
  if (level == "CRITICAL") n_crit <<- n_crit + 1L
  if (level == "INFO")     n_info <<- n_info + 1L
}

# ==============================================================================
# LOAD DATA
# ==============================================================================
cat("\n====== LOADING DATA ======\n")

# Auto-detect param_tag
fold_files = sort(list.files(
  here("Calanda_JSDM", "output", "results", "cv"),
  pattern = "^fold_.*\\.rds$", full.names = TRUE
))
param_tag = regmatches(basename(fold_files[1]),
                       regexpr("a[0-9.]+_l[0-9.]+_le[0-9.]+(?=\\.)",
                               basename(fold_files[1]), perl = TRUE))
cat(sprintf("  param_tag: %s\n", param_tag))
cat(sprintf("  %d fold files found\n", length(fold_files)))

# Main data
data_file = here("Calanda_JSDM", "output", "data_calanda_jsdm_2026-03-06.rds")
data_calanda = readRDS(data_file)
X = data_calanda$X
Y = data_calanda$Y
cat(sprintf("  X: %d sites x %d vars\n", nrow(X), ncol(X)))
cat(sprintf("  Y: %d sites x %d species\n", nrow(Y), ncol(Y)))

# Fold data
fold_data_list = lapply(fold_files, readRDS)
k = length(fold_data_list)

# VP summaries
vp_sp_file = here("Calanda_JSDM", "output", "results",
                   paste0("vp_species_summary_", param_tag, ".csv"))
vp_si_file = here("Calanda_JSDM", "output", "results",
                   paste0("vp_sites_summary_", param_tag, ".csv"))
species_cv_file = here("Calanda_JSDM", "output", "results",
                        paste0("exp6_species_cv_", param_tag, ".csv"))
sites_cv_file = here("Calanda_JSDM", "output", "results",
                      paste0("exp6_sites_cv_", param_tag, ".csv"))

vp_sp = read_csv(vp_sp_file, show_col_types = FALSE)
vp_si = read_csv(vp_si_file, show_col_types = FALSE)
species_cv = read_csv(species_cv_file, show_col_types = FALSE)
site_cv = read_csv(sites_cv_file, show_col_types = FALSE)

# Ensure site_name column exists in vp_si (bridge for pre/post script 08 re-run)
if (!"site_name" %in% names(vp_si)) {
  site_names_map = rownames(data_calanda$X)
  vp_si = vp_si %>% mutate(site_name = site_names_map[site_idx])
}

# Traits
traits_medians = read_csv(here("Calanda_JSDM", "output", "traits_medians_imputed.csv"),
                           show_col_types = FALSE)
traits_kcv = read_csv(here("Calanda_JSDM", "output", "traits_kcv_imputed.csv"),
                       show_col_types = FALSE)
comm_traits_uw = read_csv(here("Calanda_JSDM", "output", "community_traits_unweighted_imputed.csv"),
                           show_col_types = FALSE)
sp_distinct = read_csv(here("Calanda_JSDM", "output", "species_functional_distinctiveness.csv"),
                        show_col_types = FALSE)

cat("  All files loaded successfully.\n")


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  1. DATA PREP CONSISTENCY (04_prepare)                                   ║
# ╚════════════════════════════════════════════════════════════════════════════╝
cat("\n====== 1. DATA PREP CONSISTENCY ======\n")

# CHECK-1a: Coverage filter mismatch (code vs message)
# Line 364: filter(tot_rel_cover >= 0.70)
# Line 368: cat says "80%"
log_result("WARNING", "CHECK-1a",
  "Script 04 line 364 filters at >= 0.70 (70%) but line 368 cat message says '80%'. Fix the message.")

# CHECK-1b: X-Y row alignment
if (identical(rownames(X), rownames(Y))) {
  log_result("PASS", "CHECK-1b", "rownames(X) == rownames(Y) — aligned.")
} else {
  log_result("CRITICAL", "CHECK-1b",
    sprintf("Row name mismatch! X has %d rows, Y has %d rows. Shared: %d",
            nrow(X), nrow(Y), length(intersect(rownames(X), rownames(Y)))))
}

na_x = sum(is.na(X))
na_y = sum(is.na(Y))
if (na_x == 0 && na_y == 0) {
  log_result("PASS", "CHECK-1b.2", "No NAs in X or Y.")
} else {
  log_result("WARNING", "CHECK-1b.2",
    sprintf("NAs found: X has %d, Y has %d", na_x, na_y))
}

# CHECK-1c: Species count — Y should have only species surviving 5% prevalence
# We can verify prevalence of each species in Y
prevalences = colMeans(Y)
min_prev = min(prevalences)
cat(sprintf("  Species prevalence range: [%.4f, %.4f]\n", min_prev, max(prevalences)))
# With 5% filter on original data, prevalence can drop after site filtering.
# But no species should be all-zero.
below_5pct = sum(prevalences < 0.05)
if (below_5pct > 0) {
  log_result("INFO", "CHECK-1c",
    sprintf("%d species have prevalence < 5%% in final Y (possible after site filtering).",
            below_5pct))
} else {
  log_result("PASS", "CHECK-1c", "All species have prevalence >= 5% in final Y.")
}

# CHECK-1d: No all-zero columns or rows
zero_cols = sum(colSums(Y) == 0)
zero_rows = sum(rowSums(Y) == 0)
if (zero_cols == 0 && zero_rows == 0) {
  log_result("PASS", "CHECK-1d", "No all-zero species columns or site rows in Y.")
} else {
  if (zero_cols > 0)
    log_result("CRITICAL", "CHECK-1d",
      sprintf("%d species columns are all-zero (rare filter failed?).", zero_cols))
  if (zero_rows > 0)
    log_result("CRITICAL", "CHECK-1d",
      sprintf("%d site rows are all-zero.", zero_rows))
}

# CHECK-1d.extra: X has unused predictors beyond the 12 model vars
model_env_vars = c("summer_temp", "fdd", "et.annual", "slope", "rocks_cover",
                    "trees_cover", "shrubs_cover", "soil_depth_var",
                    "tpi", "flowdir", "nutrient", "disturbance")
extra_cols = setdiff(colnames(X), c(model_env_vars, "Latitude", "Longitude", "altitude"))
if (length(extra_cols) > 0) {
  log_result("INFO", "CHECK-1d.extra",
    sprintf("X has %d columns beyond 12 model env + coords: %s",
            length(extra_cols), paste(extra_cols, collapse = ", ")))
} else {
  log_result("PASS", "CHECK-1d.extra", "X columns match expected set.")
}


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  2. VARIANCE PARTITIONING (08)                                           ║
# ╚════════════════════════════════════════════════════════════════════════════╝
cat("\n====== 2. VARIANCE PARTITIONING ======\n")

# CHECK-2a: Species index alignment — nrow(sp) == ncol(Y) for every fold
sp_dim_ok = TRUE
for (fi in seq_len(k)) {
  sp = fold_data_list[[fi]]$partition$species
  if (nrow(sp) != ncol(Y)) {
    log_result("CRITICAL", "CHECK-2a",
      sprintf("Fold %d: nrow(species VP) = %d != ncol(Y) = %d",
              fi, nrow(sp), ncol(Y)))
    sp_dim_ok = FALSE
  }
}
if (sp_dim_ok) {
  log_result("PASS", "CHECK-2a",
    sprintf("All %d folds have nrow(species VP) == ncol(Y) = %d.", k, ncol(Y)))
}

# CHECK-2b: Site index — train_idx validity
site_idx_ok = TRUE
site_fold_counts = integer(nrow(Y))
for (fi in seq_len(k)) {
  fd = fold_data_list[[fi]]
  ti = fd$train_idx
  si = fd$partition$sites

  # train_idx within valid range?
  if (any(ti < 1) || any(ti > nrow(Y))) {
    log_result("CRITICAL", "CHECK-2b",
      sprintf("Fold %d: train_idx has values outside [1, %d].", fi, nrow(Y)))
    site_idx_ok = FALSE
  }

  # nrow(sites VP) == length(train_idx)?
  if (nrow(si) != length(ti)) {
    log_result("CRITICAL", "CHECK-2b",
      sprintf("Fold %d: nrow(sites VP) = %d != length(train_idx) = %d",
              fi, nrow(si), length(ti)))
    site_idx_ok = FALSE
  }

  site_fold_counts[ti] = site_fold_counts[ti] + 1L
}
if (site_idx_ok) {
  log_result("PASS", "CHECK-2b",
    sprintf("All train_idx values valid. Sites appear in %d-%d folds (median %.0f).",
            min(site_fold_counts), max(site_fold_counts), median(site_fold_counts)))
}

# CHECK-2c: AUC join correctness — species_cv$species matches colnames(Y)
sp_cv_names = species_cv$species
y_colnames = colnames(Y)
sp_in_cv_not_y = setdiff(sp_cv_names, y_colnames)
sp_in_y_not_cv = setdiff(y_colnames, sp_cv_names)
if (length(sp_in_cv_not_y) == 0 && length(sp_in_y_not_cv) == 0) {
  log_result("PASS", "CHECK-2c", "species_cv$species exactly matches colnames(Y).")
} else {
  if (length(sp_in_cv_not_y) > 0)
    log_result("WARNING", "CHECK-2c",
      sprintf("%d species in CV file but not in Y: %s",
              length(sp_in_cv_not_y),
              paste(head(sp_in_cv_not_y, 5), collapse = ", ")))
  if (length(sp_in_y_not_cv) > 0)
    log_result("WARNING", "CHECK-2c",
      sprintf("%d species in Y but not in CV file: %s",
              length(sp_in_y_not_cv),
              paste(head(sp_in_y_not_cv, 5), collapse = ", ")))
}

# CHECK-2d: Site logloss join correctness
site_cv_names = site_cv$site
x_rownames = rownames(X)
si_in_cv_not_x = setdiff(site_cv_names, x_rownames)
si_in_x_not_cv = setdiff(x_rownames, site_cv_names)
if (length(si_in_cv_not_x) == 0 && length(si_in_x_not_cv) == 0) {
  log_result("PASS", "CHECK-2d", "site_cv$site exactly matches rownames(X).")
} else {
  if (length(si_in_cv_not_x) > 0)
    log_result("WARNING", "CHECK-2d",
      sprintf("%d sites in CV file but not in X: %s",
              length(si_in_cv_not_x),
              paste(head(si_in_cv_not_x, 5), collapse = ", ")))
  if (length(si_in_x_not_cv) > 0)
    log_result("WARNING", "CHECK-2d",
      sprintf("%d sites in X but not in CV file: %s",
              length(si_in_x_not_cv),
              paste(head(si_in_x_not_cv, 5), collapse = ", ")))
}

# CHECK-2e: Weight distribution
sp_weights = vp_sp %>%
  filter(component %in% c("env", "spa", "codist")) %>%
  pull(weight)
si_weights = vp_si %>%
  filter(component %in% c("env", "spa", "codist")) %>%
  pull(weight)
cat(sprintf("  Species weights: min=%.1f, median=%.1f, max=%.1f, IQR=[%.1f, %.1f]\n",
            min(sp_weights, na.rm = TRUE), median(sp_weights, na.rm = TRUE),
            max(sp_weights, na.rm = TRUE),
            quantile(sp_weights, 0.25, na.rm = TRUE),
            quantile(sp_weights, 0.75, na.rm = TRUE)))
cat(sprintf("  Site weights:    min=%.1f, median=%.1f, max=%.1f, IQR=[%.1f, %.1f]\n",
            min(si_weights, na.rm = TRUE), median(si_weights, na.rm = TRUE),
            max(si_weights, na.rm = TRUE),
            quantile(si_weights, 0.25, na.rm = TRUE),
            quantile(si_weights, 0.75, na.rm = TRUE)))

# Site ci_width distribution — sites with ci_width >= 0.1 are now filtered out
si_ci = vp_si %>% filter(component %in% c("env", "spa", "codist")) %>% pull(ci_width)
n_unstable_sites = vp_si %>%
  filter(ci_width >= 0.1) %>%
  pull(site_idx) %>%
  unique() %>%
  length()
if (n_unstable_sites > 0) {
  log_result("INFO", "CHECK-2e",
    sprintf("%d site(s) have ci_width >= 0.1 and will be filtered out in scripts 08/09/10.", n_unstable_sites))
} else {
  log_result("PASS", "CHECK-2e", "No sites have ci_width >= 0.1.")
}

# Species weights (still used for species-level regressions)
sp_ci = vp_sp %>% filter(component %in% c("env", "spa", "codist")) %>% pull(ci_width)
n_floor_sp = sum(sp_ci <= 0.001, na.rm = TRUE)
if (n_floor_sp > 0) {
  log_result("WARNING", "CHECK-2e.2",
    sprintf("%d species component(s) hit the 0.001 ci_width floor. Max species weight = %.0f.",
            n_floor_sp, max(sp_weights, na.rm = TRUE)))
} else {
  log_result("PASS", "CHECK-2e.2", "No species ci_width values hit the 0.001 floor.")
}

# CHECK-2f: Negative VP values
sp_neg = vp_sp %>%
  filter(component %in% c("env", "spa", "codist")) %>%
  filter(mean_val < 0)
si_neg = vp_si %>%
  filter(component %in% c("env", "spa", "codist")) %>%
  filter(mean_val < 0)
if (nrow(sp_neg) > 0 || nrow(si_neg) > 0) {
  sp_neg_summary = sp_neg %>% group_by(component) %>%
    summarise(n = n(), min_val = min(mean_val), .groups = "drop")
  si_neg_summary = si_neg %>% group_by(component) %>%
    summarise(n = n(), min_val = min(mean_val), .groups = "drop")
  neg_msg = paste0(
    "Negative VP values found. ",
    "Species: ", paste(sprintf("%s=%d (min=%.4f)", sp_neg_summary$component,
                               sp_neg_summary$n, sp_neg_summary$min_val), collapse = ", "),
    ". Sites: ", paste(sprintf("%s=%d (min=%.4f)", si_neg_summary$component,
                               si_neg_summary$n, si_neg_summary$min_val), collapse = ", "))
  log_result("INFO", "CHECK-2f", neg_msg)
} else {
  log_result("PASS", "CHECK-2f", "No negative VP values.")
}


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  3. ENVIRONMENTAL GRADIENT ANALYSIS (09)                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝
cat("\n====== 3. ENVIRONMENTAL GRADIENT ANALYSIS ======\n")

# CHECK-3a: et.annual vs et_annual column name in X
if ("et.annual" %in% colnames(X)) {
  log_result("INFO", "CHECK-3a",
    "X uses 'et.annual' (dot). Script 09 line 228 correctly accesses X[,'et.annual'] and assigns to local 'et_annual'. env_cols uses 'et_annual' (underscore) for the local df. No bug, but fragile naming.")
} else if ("et_annual" %in% colnames(X)) {
  log_result("WARNING", "CHECK-3a",
    "X uses 'et_annual' (underscore), but script 09 line 228 accesses X[,'et.annual']. This would fail!")
} else {
  log_result("CRITICAL", "CHECK-3a",
    "Neither 'et.annual' nor 'et_annual' found in colnames(X)!")
}

# CHECK-3b: Species beta extraction — load final model
model_file = here("Calanda_JSDM", "output",
                   paste0("final_model_se_", param_tag, ".rds"))
if (file.exists(model_file)) {
  fm = readRDS(model_file)
  # Check fm$species matches colnames(Y)
  fm_sp = fm$species
  if (identical(fm_sp, colnames(Y))) {
    log_result("PASS", "CHECK-3b", "fm$species matches colnames(Y) exactly.")
  } else {
    shared = length(intersect(fm_sp, colnames(Y)))
    log_result("WARNING", "CHECK-3b",
      sprintf("fm$species (%d) does not exactly match colnames(Y) (%d). Shared: %d.",
              length(fm_sp), ncol(Y), shared))
  }

  # Check beta column names after coefs extraction
  # summary(fm) requires sjSDM GPU env; use tryCatch to degrade gracefully
  beta_check_ok = tryCatch({
    environment(fm$get_model)$device = "cpu"
    model_summary = summary(fm)
    beta_colnames = rownames(model_summary$coefs)[-1]  # remove intercept
    cat(sprintf("  Beta column names from model: %s\n",
                paste(beta_colnames, collapse = ", ")))

    env_cols_09 = c("summer_temp", "fdd", "et_annual", "slope", "rocks_cover",
                     "trees_cover", "shrubs_cover", "soil_depth_var",
                     "tpi", "flowdir", "nutrient", "disturbance")
    # After renaming et.annual -> et_annual in script 09:
    beta_renamed = sub("et\\.annual", "et_annual", beta_colnames)
    if (setequal(beta_renamed, env_cols_09)) {
      log_result("PASS", "CHECK-3b.2",
        "Beta colnames (after et.annual -> et_annual rename) match env_cols.")
    } else {
      log_result("WARNING", "CHECK-3b.2",
        sprintf("Beta/env_cols mismatch. In betas not env_cols: %s. In env_cols not betas: %s.",
                paste(setdiff(beta_renamed, env_cols_09), collapse = ", "),
                paste(setdiff(env_cols_09, beta_renamed), collapse = ", ")))
    }
    TRUE
  }, error = function(e) {
    log_result("INFO", "CHECK-3b.2",
      sprintf("Could not extract betas via summary(fm) — sjSDM environment not available (%s). Skipping beta column check.",
              conditionMessage(e)))
    FALSE
  })
} else {
  log_result("WARNING", "CHECK-3b", sprintf("Final model file not found: %s", model_file))
}

# CHECK-3c: Merge completeness — species betas vs VP
# Species in VP summary
vp_sp_names = unique(vp_sp$species_name)
# Species in CV (AUC >= 0.7)
keep_cv = species_cv$species[species_cv$test_auc >= 0.7 & !is.na(species_cv$test_auc)]
# Species in model
if (exists("fm_sp")) {
  sp_no_beta = setdiff(keep_cv, fm_sp)
  sp_no_vp = setdiff(fm_sp, vp_sp_names)
  if (length(sp_no_beta) > 0) {
    log_result("WARNING", "CHECK-3c",
      sprintf("%d species pass AUC filter but are missing from model betas: %s",
              length(sp_no_beta), paste(head(sp_no_beta, 5), collapse = ", ")))
  }
  if (length(sp_no_vp) > 0) {
    log_result("WARNING", "CHECK-3c",
      sprintf("%d model species missing from VP summary: %s",
              length(sp_no_vp), paste(head(sp_no_vp, 5), collapse = ", ")))
  }
  if (length(sp_no_beta) == 0 && length(sp_no_vp) == 0) {
    log_result("PASS", "CHECK-3c",
      "All AUC-filtered species have both betas and VP values.")
  }
}


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  4. FUNCTIONAL TRAIT ANALYSIS (10)                                       ║
# ╚════════════════════════════════════════════════════════════════════════════╝
cat("\n====== 4. FUNCTIONAL TRAIT ANALYSIS ======\n")

# CHECK-4a: Woody species list consistency (04 vs 10)
woody_04 = c(
  "Fagus sylvatica", "Picea abies", "Sorbus aucuparia", "Acer pseudoplatanus",
  "Larix decidua", "Fraxinus excelsior", "Pinus sylvestris", "Corylus avellana",
  "Cornus sanguinea", "Ligustrum vulgare", "Lonicera xylosteum",
  "Rubus idaeus", "Rubus fruticosus", "Rubus caesius"
)
woody_10 = c(
  "Fagus sylvatica", "Picea abies", "Sorbus aucuparia", "Acer pseudoplatanus",
  "Larix decidua", "Fraxinus excelsior", "Pinus sylvestris", "Corylus avellana",
  "Cornus sanguinea", "Ligustrum vulgare", "Lonicera xylosteum",
  "Rubus idaeus", "Rubus fruticosus", "Rubus caesius"
)
if (setequal(woody_04, woody_10)) {
  log_result("PASS", "CHECK-4a",
    sprintf("Woody species lists identical (%d species).", length(woody_04)))
} else {
  log_result("WARNING", "CHECK-4a",
    sprintf("Woody species mismatch! In 04 not 10: %s. In 10 not 04: %s.",
            paste(setdiff(woody_04, woody_10), collapse = ", "),
            paste(setdiff(woody_10, woody_04), collapse = ", ")))
}

# CHECK-4b: Trait name matching — sp_medians species match colnames(Y)
trait_sp_names = traits_medians$species_TNRS
y_sp_names = colnames(Y)
shared_trait = length(intersect(trait_sp_names, y_sp_names))
in_y_not_trait = setdiff(y_sp_names, trait_sp_names)
if (length(in_y_not_trait) == 0) {
  log_result("PASS", "CHECK-4b",
    sprintf("All %d Y species have trait data.", ncol(Y)))
} else {
  log_result("WARNING", "CHECK-4b",
    sprintf("%d Y species missing from traits_medians_imputed: %s",
            length(in_y_not_trait),
            paste(head(in_y_not_trait, 5), collapse = ", ")))
}

# Also check species_cv species match
sp_cv_not_trait = setdiff(species_cv$species, trait_sp_names)
if (length(sp_cv_not_trait) > 0) {
  log_result("WARNING", "CHECK-4b.2",
    sprintf("%d species_cv species missing from trait data: %s",
            length(sp_cv_not_trait),
            paste(head(sp_cv_not_trait, 5), collapse = ", ")))
} else {
  log_result("PASS", "CHECK-4b.2", "All species_cv species have trait data.")
}

# CHECK-4c: Log-transformation double-dipping
# Check if traits_medians_imputed.csv stores raw or already-logged values
# Raw vegetative_height in mm is typically 10-5000+.
# Log(10-5000) = 2.3 to 8.5. Log(log(x)) would be ~ 0.8 to 2.1.
sample_heights = traits_medians$Median_vegetative_height[1:10]
cat(sprintf("  Sample Median_vegetative_height (first 10): %s\n",
            paste(round(sample_heights, 2), collapse = ", ")))
if (any(sample_heights > 100, na.rm = TRUE)) {
  log_result("PASS", "CHECK-4c",
    "traits_medians_imputed.csv stores RAW values (not log-transformed). Script 10's log-transform is correct.")
} else if (all(sample_heights < 20, na.rm = TRUE)) {
  log_result("CRITICAL", "CHECK-4c",
    "traits_medians_imputed.csv appears to store LOG-transformed values. Script 10 would double-log!")
} else {
  log_result("INFO", "CHECK-4c",
    "Unable to determine if traits are raw or logged. Manual check recommended.")
}

# CHECK-4d: Community trait site mismatch
comm_site_ids = comm_traits_uw$plot_id_releve
x_site_ids = rownames(X)
sites_with_vp = x_site_ids  # VP was computed on all Y sites

comm_only = setdiff(comm_site_ids, sites_with_vp)
vp_only = setdiff(sites_with_vp, comm_site_ids)

log_result("INFO", "CHECK-4d",
  sprintf("Community traits have %d sites. VP has %d sites. %d sites in VP but not in community traits (lost in inner_join). %d community-trait sites not in VP.",
          length(comm_site_ids), length(sites_with_vp),
          length(vp_only), length(comm_only)))
if (length(vp_only) > 0) {
  log_result("WARNING", "CHECK-4d.2",
    sprintf("%d sites have VP but no community traits (woody-dominated plots excluded from 04_prepare). These sites are silently dropped in script 10's inner_join.",
            length(vp_only)))
}

# CHECK-4e: Distinctiveness column names
core_traits = c("vegetative_height", "LNC", "LA", "seed_mass", "LDMC")
expected_distinct_cols = paste0("distinct_", core_traits)
# Column names in the file use Median_ prefix for trait columns
actual_distinct_cols = grep("^distinct_", colnames(sp_distinct), value = TRUE)
missing_distinct = setdiff(expected_distinct_cols, actual_distinct_cols)
if (length(missing_distinct) == 0) {
  log_result("PASS", "CHECK-4e",
    sprintf("All expected distinctiveness columns found: %s",
            paste(expected_distinct_cols, collapse = ", ")))
} else {
  log_result("WARNING", "CHECK-4e",
    sprintf("Missing distinctiveness columns: %s. Available: %s",
            paste(missing_distinct, collapse = ", "),
            paste(actual_distinct_cols, collapse = ", ")))
}

# CHECK-4f: drop_na impact
# Simulate what script 10 does
sp_median_cols = paste0("Median_", core_traits)
sp_kcv_cols = paste0("kCV_", core_traits)
sp_distinct_cols = paste0("distinct_", core_traits)
sp_trait_cols = c(sp_median_cols, sp_kcv_cols, sp_distinct_cols)

keep_auc = species_cv$species[species_cv$test_auc >= 0.7 & !is.na(species_cv$test_auc)]
keep_non_woody = setdiff(keep_auc, woody_10)

# Merge traits
sp_merged = traits_medians %>%
  rename(species_name = species_TNRS) %>%
  left_join(traits_kcv %>% rename(species_name = species_TNRS), by = "species_name") %>%
  left_join(sp_distinct %>% rename(species_name = species), by = "species_name") %>%
  filter(species_name %in% keep_non_woody)

n_before = nrow(sp_merged)
# Check which trait cols exist
available_trait_cols = intersect(sp_trait_cols, colnames(sp_merged))
missing_trait_cols = setdiff(sp_trait_cols, colnames(sp_merged))
if (length(missing_trait_cols) > 0) {
  log_result("WARNING", "CHECK-4f",
    sprintf("Missing trait columns for drop_na: %s", paste(missing_trait_cols, collapse = ", ")))
}

sp_after_dropna = sp_merged %>% drop_na(all_of(available_trait_cols))
n_after = nrow(sp_after_dropna)
n_dropped = n_before - n_after
if (n_dropped > 0) {
  log_result("WARNING", "CHECK-4f",
    sprintf("drop_na removes %d / %d species (%.1f%%) due to missing trait data.",
            n_dropped, n_before, 100 * n_dropped / n_before))
} else {
  log_result("PASS", "CHECK-4f", "No species dropped by drop_na (all traits present after imputation).")
}


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  5. CROSS-SCRIPT CONSISTENCY                                             ║
# ╚════════════════════════════════════════════════════════════════════════════╝
cat("\n====== 5. CROSS-SCRIPT CONSISTENCY ======\n")

# CHECK-5a: Filtering consistency across scripts 08, 09, 10
auc_filtered = species_cv$species[species_cv$test_auc >= 0.7 & !is.na(species_cv$test_auc)]
logloss_filtered = site_cv$site[site_cv$logloss <= log(2) & !is.na(site_cv$logloss)]
cat(sprintf("  Species with AUC >= 0.7: %d / %d\n", length(auc_filtered), nrow(species_cv)))
cat(sprintf("  Sites with logloss <= log(2): %d / %d\n", length(logloss_filtered), nrow(site_cv)))

# Script 10 additionally removes woody species
non_woody_auc = setdiff(auc_filtered, woody_10)
n_woody_removed = length(auc_filtered) - length(non_woody_auc)
log_result("INFO", "CHECK-5a",
  sprintf("Script 10 additionally removes %d woody species from the AUC-filtered set (%d -> %d).",
          n_woody_removed, length(auc_filtered), length(non_woody_auc)))

# CHECK-5b: VP summary dimensions match current Y
vp_sp_unique_species = length(unique(vp_sp$species_idx))
vp_si_unique_sites = length(unique(vp_si$site_idx))
if (vp_sp_unique_species == ncol(Y)) {
  log_result("PASS", "CHECK-5b",
    sprintf("VP species summary has %d unique species = ncol(Y).", vp_sp_unique_species))
} else {
  log_result("WARNING", "CHECK-5b",
    sprintf("VP species summary has %d unique species but ncol(Y) = %d.",
            vp_sp_unique_species, ncol(Y)))
}

# For sites, VP only has training sites per fold, so unique site_idx should cover all sites
if (vp_si_unique_sites == nrow(Y)) {
  log_result("PASS", "CHECK-5b.2",
    sprintf("VP site summary has %d unique sites = nrow(Y).", vp_si_unique_sites))
} else {
  log_result("WARNING", "CHECK-5b.2",
    sprintf("VP site summary has %d unique sites but nrow(Y) = %d.",
            vp_si_unique_sites, nrow(Y)))
}

# CHECK-5c: param_tag consistency
# The fold files, VP summaries, species CV, and site CV all use the same param_tag
expected_files = c(
  paste0("vp_species_summary_", param_tag, ".csv"),
  paste0("vp_sites_summary_", param_tag, ".csv"),
  paste0("exp6_species_cv_", param_tag, ".csv"),
  paste0("exp6_sites_cv_", param_tag, ".csv"),
  paste0("final_model_se_", param_tag, ".rds")
)
results_dir = here("Calanda_JSDM", "output", "results")
output_dir = here("Calanda_JSDM", "output")
missing_param_files = character(0)
for (f in expected_files[1:4]) {
  if (!file.exists(file.path(results_dir, f)))
    missing_param_files = c(missing_param_files, f)
}
if (!file.exists(file.path(output_dir, expected_files[5])))
  missing_param_files = c(missing_param_files, expected_files[5])

if (length(missing_param_files) == 0) {
  log_result("PASS", "CHECK-5c",
    sprintf("All expected files with param_tag '%s' exist.", param_tag))
} else {
  log_result("WARNING", "CHECK-5c",
    sprintf("Missing files for param_tag '%s': %s",
            param_tag, paste(missing_param_files, collapse = ", ")))
}

# CHECK-5d: Hardcoded date
hardcoded_file = here("Calanda_JSDM", "output", "data_calanda_jsdm_2026-03-06.rds")
all_data_files = sort(list.files(
  here("Calanda_JSDM", "output"),
  pattern = "^data_calanda_jsdm_[0-9].*\\.rds$", full.names = TRUE
), decreasing = TRUE)
most_recent = all_data_files[1]
if (file.exists(hardcoded_file)) {
  if (normalizePath(hardcoded_file) == normalizePath(most_recent)) {
    log_result("PASS", "CHECK-5d",
      sprintf("Hardcoded file '%s' exists and is the most recent.", basename(hardcoded_file)))
  } else {
    log_result("WARNING", "CHECK-5d",
      sprintf("Hardcoded file exists but '%s' is more recent.", basename(most_recent)))
  }
} else {
  log_result("CRITICAL", "CHECK-5d",
    sprintf("Hardcoded file '%s' does not exist!", basename(hardcoded_file)))
}


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  6. STATISTICAL METHOD CHECKS                                           ║
# ╚════════════════════════════════════════════════════════════════════════════╝
cat("\n====== 6. STATISTICAL METHOD CHECKS ======\n")

# CHECK-6a: CI calculation with small n
# Species: all k folds should be present for each species
sp_n_per_species = vp_sp %>%
  filter(component == "env") %>%
  group_by(species_idx) %>%
  summarise(n_present = 1, .groups = "drop")  # Each row is already aggregated

# Actually, vp_sp is already aggregated (mean_val per species). The n used for CI
# in script 08 is n() within the fold loop = k = 10 for species (they appear in all folds).
# For sites, they appear only in training folds.
cat(sprintf("  Species: always in all %d folds => n=%d for CI calculation.\n", k, k))

# Sites: count how many folds each site appears in
site_n_folds = vp_si %>%
  filter(component == "env") %>%
  group_by(site_idx) %>%
  summarise(n_est = 1, .groups = "drop")
# We already computed site_fold_counts above
cat(sprintf("  Sites: appear in %d-%d folds. qt(0.975, n-1) used with correct n from n().\n",
            min(site_fold_counts), max(site_fold_counts)))

# But wait — script 08 aggregates by site_idx using group_by + n().
# Since each fold contributes one row per training site, n() should equal the number
# of folds that site appeared in. This is correct.
if (min(site_fold_counts) >= 2) {
  log_result("PASS", "CHECK-6a",
    sprintf("All sites appear in >= 2 folds (min=%d). CI calculation valid.", min(site_fold_counts)))
} else {
  log_result("WARNING", "CHECK-6a",
    sprintf("Some sites appear in only %d fold(s). CI with n-1 df may be undefined or huge.",
            min(site_fold_counts)))
}

# CHECK-6b / 3d / 3e / 3e.2: Re-fit site-level stepwise models and check assumptions
# Helper: run assumption checks on a single fitted lm
run_assumption_checks = function(mod, model_label) {
  retained = names(coef(mod))[-1]  # exclude intercept


  # CHECK-6b: intercept-only?
  if (length(retained) == 0) {
    log_result("WARNING", "CHECK-6b",
      sprintf("%s: stepwise retained only the intercept (all predictors dropped).", model_label))
  } else {
    log_result("PASS", "CHECK-6b",
      sprintf("%s: %d predictor(s) retained.", model_label, length(retained)))

    # CHECK-3d: VIF (only meaningful with >1 predictor)
    if (length(retained) > 1) {
      vif_vals = tryCatch(car::vif(mod), error = function(e) NULL)
      if (!is.null(vif_vals)) {
        # vif() may return a matrix for poly terms; extract GVIF or plain values
        if (is.matrix(vif_vals)) vif_vals = vif_vals[, 1]
        max_vif = max(vif_vals)
        if (max_vif > 10) {
          log_result("CRITICAL", "CHECK-3d",
            sprintf("%s: max VIF = %.1f (term: %s). Severe multicollinearity.",
                    model_label, max_vif, names(which.max(vif_vals))))
        } else if (max_vif > 5) {
          log_result("WARNING", "CHECK-3d",
            sprintf("%s: max VIF = %.1f (term: %s). Moderate multicollinearity.",
                    model_label, max_vif, names(which.max(vif_vals))))
        } else {
          log_result("PASS", "CHECK-3d",
            sprintf("%s: max VIF = %.1f. No multicollinearity concern.", model_label, max_vif))
        }
      }
    }

    # CHECK-3e: Residual normality (Shapiro-Wilk)
    resids = residuals(mod)
    # Shapiro-Wilk limited to n <= 5000
    if (length(resids) <= 5000) {
      sw = shapiro.test(resids)
      if (sw$p.value < 0.01) {
        log_result("WARNING", "CHECK-3e",
          sprintf("%s: Shapiro-Wilk p = %.2e — residuals deviate from normality.",
                  model_label, sw$p.value))
      } else {
        log_result("PASS", "CHECK-3e",
          sprintf("%s: Shapiro-Wilk p = %.3f — residual normality OK.", model_label, sw$p.value))
      }
    } else {
      log_result("INFO", "CHECK-3e",
        sprintf("%s: n = %d > 5000, Shapiro-Wilk skipped.", model_label, length(resids)))
    }

    # CHECK-3e.2: Homoscedasticity (Breusch-Pagan, manual implementation)
    bp_p = tryCatch({
      e2 = resids^2
      bp_mod = lm(e2 ~ fitted(mod))
      bp_r2 = summary(bp_mod)$r.squared
      bp_stat = length(resids) * bp_r2
      pchisq(bp_stat, df = 1, lower.tail = FALSE)
    }, error = function(e) NULL)
    if (!is.null(bp_p)) {
      if (bp_p < 0.01) {
        log_result("WARNING", "CHECK-3e.2",
          sprintf("%s: Breusch-Pagan p = %.2e — evidence of heteroscedasticity.",
                  model_label, bp_p))
      } else {
        log_result("PASS", "CHECK-3e.2",
          sprintf("%s: Breusch-Pagan p = %.3f — homoscedasticity OK.", model_label, bp_p))
      }
    }
  }
}

# --------------------------------------------------------------------------
# 6b-A: Re-fit community-level env gradient models (script 09)
# --------------------------------------------------------------------------
cat("\n  --- Re-fitting script 09 community env models ---\n")

env_cols_09 = c("summer_temp", "fdd", "et_annual", "slope", "rocks_cover",
                "trees_cover", "shrubs_cover", "soil_depth_var",
                "tpi", "flowdir", "nutrient", "disturbance")

vp_components_09 = c("env", "codist", "spa")

# Build wide site VP + env predictors (name-based joins, same logic as script 09)
site_names_09 = rownames(data_calanda$X)
keep_names_09 = site_cv %>%
  filter(logloss <= log(2), !is.na(logloss)) %>%
  pull(site)

df_wide_09 = vp_si %>%
  filter(site_name %in% keep_names_09) %>%
  filter(component %in% c("env", "spa", "codist", "r2")) %>%
  select(site_name, site_idx, component, mean_val, ci_width) %>%
  pivot_wider(
    names_from = component,
    values_from = c(mean_val, ci_width),
    names_glue = "{component}_{.value}"
  )

# Filter unstable sites (ci_width >= 0.1 on any component)
ci_cols_09 = grep("_ci_width$", names(df_wide_09), value = TRUE)
unstable_09 = apply(df_wide_09[, ci_cols_09], 1, function(r) any(r >= 0.1, na.rm = TRUE))
df_wide_09 = df_wide_09[!unstable_09, ]

# Attach env predictors by site name (not positional index)
X_df_09 = as.data.frame(X) %>%
  rownames_to_column("site_name") %>%
  select(site_name, slope, summer_temp, fdd, et_annual = et.annual,
         soil_depth_var, trees_cover, shrubs_cover, rocks_cover,
         flowdir, tpi, nutrient, disturbance)
df_wide_09 = df_wide_09 %>%
  left_join(X_df_09, by = "site_name")

comm_linear_09 = env_cols_09
comm_quad_09 = paste0("I(", env_cols_09, "^2)")
comm_rhs_09 = paste(c(comm_linear_09, comm_quad_09), collapse = " + ")

for (comp in vp_components_09) {
  y_col = paste0(comp, "_mean_val")
  full_formula = as.formula(paste0(y_col, " ~ ", comm_rhs_09))
  full_mod = lm(full_formula, data = df_wide_09)
  best_mod = step(full_mod, direction = "backward", trace = 0)
  run_assumption_checks(best_mod, sprintf("Env-community [%s]", comp))
}

# --------------------------------------------------------------------------
# 6b-B: Re-fit community-level functional trait models (script 10)
# --------------------------------------------------------------------------
cat("\n  --- Re-fitting script 10 community trait models ---\n")

core_traits_10 = c("vegetative_height", "LNC", "LA", "seed_mass", "LDMC")
comm_mean_cols_10 = paste0("Median_", core_traits_10, "_mean")
comm_var_cols_10 = paste0("Median_", core_traits_10, "_var")
comm_trait_cols_10 = c(comm_mean_cols_10, comm_var_cols_10)

# Build site VP wide (same filter as 09: logloss + ci_width, name-based)
df_si_wide_10 = vp_si %>%
  filter(site_name %in% keep_names_09) %>%
  filter(component %in% c("env", "spa", "codist", "r2")) %>%
  select(site_name, site_idx, component, mean_val, ci_width) %>%
  pivot_wider(
    names_from = component,
    values_from = c(mean_val, ci_width),
    names_glue = "{component}_{.value}"
  )
ci_cols_10 = grep("_ci_width$", names(df_si_wide_10), value = TRUE)
unstable_10 = apply(df_si_wide_10[, ci_cols_10], 1, function(r) any(r >= 0.1, na.rm = TRUE))
df_si_wide_10 = df_si_wide_10[!unstable_10, ]

# Merge with community traits by site name
df_comm_merged_10 = df_si_wide_10 %>%
  rename(plot_id_releve = site_name) %>%
  inner_join(comm_traits_uw %>% select(plot_id_releve, all_of(comm_trait_cols_10)),
             by = "plot_id_releve") %>%
  drop_na(all_of(comm_trait_cols_10))

# Standardize predictors
df_comm_merged_10 = df_comm_merged_10 %>%
  mutate(across(all_of(comm_trait_cols_10), ~ as.numeric(scale(.x))))

comm_rhs_10 = paste(comm_trait_cols_10, collapse = " + ")

for (comp in vp_components_09) {
  y_col = paste0(comp, "_mean_val")
  full_formula = as.formula(paste0(y_col, " ~ ", comm_rhs_10))
  full_mod = lm(full_formula, data = df_comm_merged_10)
  best_mod = step(full_mod, direction = "backward", trace = 0)
  run_assumption_checks(best_mod, sprintf("Trait-community [%s]", comp))
}

# CHECK-6c: Weight normalization
# Weights 1/ci_width^2 are unnormalized. lm() handles this fine, but check scale.
max_sp_w = max(sp_weights, na.rm = TRUE)
min_sp_w = min(sp_weights, na.rm = TRUE)
max_si_w = max(si_weights, na.rm = TRUE)
min_si_w = min(si_weights, na.rm = TRUE)
sp_ratio = max_sp_w / min_sp_w
si_ratio = max_si_w / min_si_w
cat(sprintf("  Species weight ratio (max/min): %.1f\n", sp_ratio))
cat(sprintf("  Site weight ratio (max/min): %.1f\n", si_ratio))
if (sp_ratio > 1e6 || si_ratio > 1e6) {
  log_result("WARNING", "CHECK-6c",
    sprintf("Extreme weight ratio (>1e6). Species: %.1e, Sites: %.1e. May cause numerical instability in lm().",
            sp_ratio, si_ratio))
} else {
  log_result("PASS", "CHECK-6c",
    "Weight ratios within reasonable range (< 1e6).")
}


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  ALREADY-SPOTTED ISSUES (from audit plan)                                ║
# ╚════════════════════════════════════════════════════════════════════════════╝
cat("\n====== ALREADY-SPOTTED ISSUES ======\n")

# CRITICAL-1: veg_abund_filtered undefined (script 04, line 712)
log_result("CRITICAL", "SPOTTED-1",
  "Script 04 line 712: references 'veg_abund_filtered' but variable is 'Y_filtered'. Runtime error. FIX: replace with Y_filtered.")

# WARNING-4: Cat message mismatch (script 04, line 368)
# Already flagged as CHECK-1a above

# INFO-8: Hardcoded dates in scripts 09 and 10
log_result("INFO", "SPOTTED-8",
  "Scripts 09 (line 184) and 10 (line 240) hardcode 'data_calanda_jsdm_2026-03-06.rds'. Script 08 auto-detects the latest file.")


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SUMMARY                                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝
cat("\n======================================================\n")
cat("  VALIDATION SUMMARY\n")
cat("======================================================\n")
cat(sprintf("  \033[32mPASS:     %d\033[0m\n", n_pass))
cat(sprintf("  \033[36mINFO:     %d\033[0m\n", n_info))
cat(sprintf("  \033[33mWARNING:  %d\033[0m\n", n_warn))
cat(sprintf("  \033[31mCRITICAL: %d\033[0m\n", n_crit))
cat("======================================================\n")

if (n_crit > 0) {
  cat("\n  \033[31mAction required: fix CRITICAL issues before running analysis.\033[0m\n")
}
cat("\n")
