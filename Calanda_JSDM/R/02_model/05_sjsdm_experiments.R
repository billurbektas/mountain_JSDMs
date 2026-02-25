# ==============================================================================
# Script: 05_sjsdm_experiments.R
# Purpose: Run sjSDM model experiments with variance partitioning
#
# After each run, extract:
#   - Overall pseudo-R² (McFadden)
#   - Full model Venn components (E, S, C from ANOVA)
#   - Species-level: R², env, spatial, codist
#   - Site-level: R², env, spatial, codist
#
# Experiments:
#   1. Full grid search (spatial_form x alpha x lambda)
#   2. Decoupled lambda sensitivity
#   3. Drop-one environmental variable
#   4. k-fold CV diagnostics (species AUC, train-test gap, site log-loss)
#   5. ANOVA sampling saturation
#
# Inputs:
#   - output/data_calanda_jsdm_<date>.rds
#
# Outputs:
#   - results/runs/*.rds
#   - results/decoupled/*.rds
#   - results/dropone/*.rds
#   - results/cv/*.rds
#   - results/summary_exp*.csv
#   - results/plots/*.pdf
#
# Requires: sjSDM (GPU), tidyverse, here, pROC
# ==============================================================================

library(sjSDM)
library(tidyverse)
library(conflicted)
library(here)
library(pROC)

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("intersect", "base")

set.seed(42)
session_info = sessionInfo()

# ==============================================================================
# USER SETTINGS
# ==============================================================================

date_data = "2026-02-25"
device = "gpu"
iterations = 750L
learning_rate = 0.01
sampling_fit = 1000L
anova_samples_fast = 5000L
woody_threshold = 10

env_formula = ~summer_temp + fdd + et.annual + slope + rocks_cover +
  trees_cover + shrubs_cover + soil_depth_var +
  tpi + flowdir + nutrient + disturbance

# Toggle experiments
run_exp1 = TRUE   # grid search
run_exp2 = TRUE   # decoupled lambdas
run_exp3 = TRUE   # drop-one predictor
run_exp4 = TRUE   # k-fold CV
run_exp5 = TRUE   # ANOVA saturation

# Representative config for Exp 3-5 (set after Exp 1, or manually)
rep_spatial = "DNN"
rep_alpha = 0.5
rep_lambda = 0.1

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("\n=== Loading data ===\n")
data_calanda_jsdm = readRDS(
  here("Calanda_JSDM", "output", paste0("data_calanda_jsdm_", date_data, ".rds"))
)
X = data_calanda_jsdm$X
Y = data_calanda_jsdm$Y

cat("X:", nrow(X), "sites x", ncol(X), "predictors\n")
cat("Y:", nrow(Y), "sites x", ncol(Y), "species\n")

altitude = X[, "altitude"]
XY = X[, c("Latitude", "Longitude")]
richness = rowSums(Y)
env_cols = all.vars(env_formula)
Xenv = X[, env_cols, drop = FALSE]

dirs = c("results", "results/runs", "results/dropone", "results/cv", "results/plots")
for (d in dirs) {
  dir.create(here("Calanda_JSDM", d), showWarnings = FALSE, recursive = TRUE)
}

# ==============================================================================
# CORE FUNCTIONS
# ==============================================================================

fit_sjsdm = function(Y_train, Xenv_train, XY_train,
                     env_form = env_formula,
                     spatial_form = "DNN",
                     lambda_env = 0.1, alpha_env = 0.5,
                     lambda_sp = 0.1, alpha_sp = 0.5,
                     lambda_bio = 0.1, alpha_bio = 0.5) {

  if (spatial_form == "DNN") {
    sp_component = DNN(
      XY_train,
      formula = ~0 + .,
      hidden = c(30L, 30L),
      activation = "relu",
      bias = FALSE,
      lambda = lambda_sp, alpha = alpha_sp
    )
  } else if (spatial_form == "LIN_XY_XY") {
    sp_component = linear(
      XY_train,
      formula = ~0 + Latitude + Longitude + I(Latitude * Longitude),
      lambda = lambda_sp, alpha = alpha_sp
    )
  } else {
    stop("Unknown spatial_form: ", spatial_form)
  }

  sjSDM(
    Y = Y_train,
    env = linear(Xenv_train, formula = env_form,
                 lambda = lambda_env, alpha = alpha_env),
    spatial = sp_component,
    biotic = bioticStruct(lambda = lambda_bio, alpha = alpha_bio,
                          df = ncol(Y_train), reg_on_Cov = FALSE),
    iter = iterations,
    device = device,
    learning_rate = learning_rate,
    sampling = sampling_fit,
    control = sjSDMControl(
      RMSprop(weight_decay = 0.0),
      scheduler = 5L,
      early_stopping_training = 25L,
      lr_reduce_factor = 0.9
    ),
    se = FALSE
  )
}


#' Extract variance partitioning: overall R², Venn components, species & site E/S/C
extract_partition = function(model, anova_samples = anova_samples_fast) {

  an = anova(model, verbose = FALSE, samples = anova_samples)
  res = internalStructure(an, fractions = "proportional",
                          Rsquared = "McFadden", plot = FALSE)

  # --- Full model Venn components (from ANOVA results table) ---
  venn = an$results
  mcf_col = which(colnames(venn$sites) == "McFadden" |
                    grepl("McFadden", colnames(venn$sites), ignore.case = TRUE))
  # Extract the 7 Venn fractions: F_A (env), F_B (biotic), F_AB, F_S, F_AS, F_BS, F_ABS
  venn_labels = c("F_A", "F_B", "F_AB", "F_S", "F_AS", "F_BS", "F_ABS")

  # --- Species-level ---
  sp = res$internals$Species
  sp$R2_total = rowSums(sp)

  # --- Site-level ---
  si = res$internals$Sites
  si$R2_total = rowSums(si)

  # --- Overall pseudo-R² ---
  overall_R2 = mean(sp$R2_total, na.rm = TRUE)

  # --- Summaries ---
  summary_row = tibble(
    overall_R2 = overall_R2,
    # Species means
    species_R2_mean = mean(sp$R2_total, na.rm = TRUE),
    species_E_mean = mean(sp$env, na.rm = TRUE),
    species_S_mean = mean(sp$spa, na.rm = TRUE),
    species_C_mean = mean(sp$codist, na.rm = TRUE),
    # Species medians
    species_R2_median = median(sp$R2_total, na.rm = TRUE),
    species_E_median = median(sp$env, na.rm = TRUE),
    species_S_median = median(sp$spa, na.rm = TRUE),
    species_C_median = median(sp$codist, na.rm = TRUE),
    # Site means
    site_R2_mean = mean(si$R2_total, na.rm = TRUE),
    site_E_mean = mean(si$env, na.rm = TRUE),
    site_S_mean = mean(si$spa, na.rm = TRUE),
    site_C_mean = mean(si$codist, na.rm = TRUE),
    # Site medians
    site_R2_median = median(si$R2_total, na.rm = TRUE),
    site_E_median = median(si$env, na.rm = TRUE),
    site_S_median = median(si$spa, na.rm = TRUE),
    site_C_median = median(si$codist, na.rm = TRUE)
  )

  list(
    anova = an,
    internal_structure = res,
    species = sp,
    sites = si,
    summary = summary_row
  )
}


build_stratified_folds = function(k = 10, alt, tree_cover, shrub_cover,
                                  disturbance, nutrient, richness,
                                  woody_thresh = woody_threshold) {
  n = length(alt)

  elev_bin = as.integer(cut(alt,
    breaks = quantile(alt, probs = seq(0, 1, length.out = 11)),
    include.lowest = TRUE))

  woody_class = ifelse((tree_cover + shrub_cover) > woody_thresh, 1L, 0L)

  landuse_score = scale(disturbance)[, 1] + scale(nutrient)[, 1]
  landuse_bin = as.integer(cut(landuse_score,
    breaks = quantile(landuse_score, probs = seq(0, 1, length.out = 6), na.rm = TRUE),
    include.lowest = TRUE))
  landuse_bin[is.na(landuse_bin)] = 1L

  richness_bin = as.integer(cut(richness,
    breaks = quantile(richness, probs = seq(0, 1, length.out = 6)),
    include.lowest = TRUE))

  strata = paste(elev_bin, woody_class, landuse_bin, richness_bin, sep = "_")
  strata_factor = as.factor(strata)
  fold_ids = rep(NA_integer_, n)

  for (s in levels(strata_factor)) {
    idx = which(strata_factor == s)
    fold_ids[idx] = sample(rep(seq_len(k), length.out = length(idx)))
  }

  lapply(seq_len(k), function(i) which(fold_ids == i))
}


# ==============================================================================
# EXPERIMENT 1: FULL GRID SEARCH
# ==============================================================================

if (run_exp1) {

  cat("\n", strrep("=", 60), "\n")
  cat("EXPERIMENT 1: Full grid search\n")
  cat(strrep("=", 60), "\n")

  grid = expand.grid(
    spatial_form = c("DNN", "LIN_XY_XY"),
    alpha_common = c(0, 0.5, 1),
    lambda_common = c(0.01, 0.1, 0.5, 1, 2),
    stringsAsFactors = FALSE
  ) %>%
    mutate(run_id = paste0(spatial_form, "_a", alpha_common, "_l", lambda_common))

  cat("Total runs:", nrow(grid), "\n\n")

  exp1_summaries = list()

  for (i in seq_len(nrow(grid))) {

    cfg = grid[i, ]
    run_file = here("Calanda_JSDM", "results", "runs", paste0(cfg$run_id, ".rds"))

    if (file.exists(run_file)) {
      cat("[", i, "/", nrow(grid), "] Cached:", cfg$run_id, "\n")
      run_data = readRDS(run_file)
    } else {
      cat("[", i, "/", nrow(grid), "] Fitting:", cfg$run_id, "\n")

      model = fit_sjsdm(
        Y_train = Y, Xenv_train = Xenv, XY_train = XY,
        spatial_form = cfg$spatial_form,
        lambda_env = cfg$lambda_common, alpha_env = cfg$alpha_common,
        lambda_sp = cfg$lambda_common, alpha_sp = cfg$alpha_common,
        lambda_bio = cfg$lambda_common, alpha_bio = cfg$alpha_common
      )

      partition = extract_partition(model)

      run_data = list(config = cfg, model = model, partition = partition)
      saveRDS(run_data, run_file)
      cat("  -> Saved\n")
    }

    exp1_summaries[[i]] = bind_cols(cfg, run_data$partition$summary)
  }

  summary_exp1 = bind_rows(exp1_summaries)
  write_csv(summary_exp1, here("Calanda_JSDM", "results", "summary_exp1_grid.csv"))
  cat("\nSaved results/summary_exp1_grid.csv\n")
}

# ==============================================================================
# EXPERIMENT 2: DECOUPLED LAMBDAS
# ==============================================================================

if (run_exp2) {

  cat("\n", strrep("=", 60), "\n")
  cat("EXPERIMENT 2: Decoupled lambda sensitivity\n")
  cat(strrep("=", 60), "\n")

  alpha_fixed = 0.5
  lambda_anchor = 0.1
  lambda_values = c(0.01, 0.1, 0.5, 1, 2)

  grid_decoupled = bind_rows(
    expand.grid(spatial_form = c("DNN", "LIN_XY_XY"),
                vary = "lambda_env", lambda_vary = lambda_values,
                stringsAsFactors = FALSE),
    expand.grid(spatial_form = c("DNN", "LIN_XY_XY"),
                vary = "lambda_sp", lambda_vary = lambda_values,
                stringsAsFactors = FALSE),
    expand.grid(spatial_form = c("DNN", "LIN_XY_XY"),
                vary = "lambda_bio", lambda_vary = lambda_values,
                stringsAsFactors = FALSE)
  ) %>%
    mutate(
      lambda_env = ifelse(vary == "lambda_env", lambda_vary, lambda_anchor),
      lambda_sp = ifelse(vary == "lambda_sp", lambda_vary, lambda_anchor),
      lambda_bio = ifelse(vary == "lambda_bio", lambda_vary, lambda_anchor),
      run_id = paste0("decoupled_", spatial_form, "_", vary, "_", lambda_vary)
    )

  dir.create(here("Calanda_JSDM", "results", "decoupled"),
             showWarnings = FALSE, recursive = TRUE)

  cat("Total runs:", nrow(grid_decoupled), "\n\n")

  exp2_summaries = list()

  for (i in seq_len(nrow(grid_decoupled))) {

    cfg = grid_decoupled[i, ]
    run_file = here("Calanda_JSDM", "results", "decoupled", paste0(cfg$run_id, ".rds"))

    if (file.exists(run_file)) {
      cat("[", i, "/", nrow(grid_decoupled), "] Cached:", cfg$run_id, "\n")
      run_data = readRDS(run_file)
    } else {
      cat("[", i, "/", nrow(grid_decoupled), "] Fitting:", cfg$run_id, "\n")

      model = fit_sjsdm(
        Y_train = Y, Xenv_train = Xenv, XY_train = XY,
        spatial_form = cfg$spatial_form,
        lambda_env = cfg$lambda_env, alpha_env = alpha_fixed,
        lambda_sp = cfg$lambda_sp, alpha_sp = alpha_fixed,
        lambda_bio = cfg$lambda_bio, alpha_bio = alpha_fixed
      )

      partition = extract_partition(model)

      run_data = list(config = cfg, model = model, partition = partition)
      saveRDS(run_data, run_file)
      cat("  -> Saved\n")
    }

    exp2_summaries[[i]] = bind_cols(
      cfg %>% select(run_id, spatial_form, vary, lambda_vary, lambda_env, lambda_sp, lambda_bio),
      run_data$partition$summary
    )
  }

  summary_exp2 = bind_rows(exp2_summaries)
  write_csv(summary_exp2, here("Calanda_JSDM", "results", "summary_exp2_decoupled.csv"))
  cat("\nSaved results/summary_exp2_decoupled.csv\n")
}

# ==============================================================================
# EXPERIMENT 3: DROP-ONE ENVIRONMENTAL VARIABLE
# ==============================================================================

if (run_exp3) {

  cat("\n", strrep("=", 60), "\n")
  cat("EXPERIMENT 3: Drop-one environmental variable\n")
  cat(strrep("=", 60), "\n")

  cat("Config:", rep_spatial, "alpha=", rep_alpha, "lambda=", rep_lambda, "\n\n")

  exp3_summaries = list()

  for (j in seq_along(env_cols)) {

    pred = env_cols[j]
    run_id = paste0("dropone_", pred)
    run_file = here("Calanda_JSDM", "results", "dropone", paste0(run_id, ".rds"))

    if (file.exists(run_file)) {
      cat("[", j, "/", length(env_cols), "] Cached:", pred, "\n")
      run_data = readRDS(run_file)
    } else {
      cat("[", j, "/", length(env_cols), "] Dropping:", pred, "\n")

      reduced_cols = setdiff(env_cols, pred)
      reduced_formula = as.formula(paste("~", paste(reduced_cols, collapse = " + ")))

      model = fit_sjsdm(
        Y_train = Y,
        Xenv_train = Xenv[, reduced_cols, drop = FALSE],
        XY_train = XY,
        env_form = reduced_formula,
        spatial_form = rep_spatial,
        lambda_env = rep_lambda, alpha_env = rep_alpha,
        lambda_sp = rep_lambda, alpha_sp = rep_alpha,
        lambda_bio = rep_lambda, alpha_bio = rep_alpha
      )

      partition = extract_partition(model)

      run_data = list(dropped = pred, model = model, partition = partition)
      saveRDS(run_data, run_file)
      cat("  -> Saved\n")
    }

    exp3_summaries[[j]] = bind_cols(
      tibble(dropped_variable = pred),
      run_data$partition$summary
    )
  }

  summary_exp3 = bind_rows(exp3_summaries)
  write_csv(summary_exp3, here("Calanda_JSDM", "results", "summary_exp3_dropone.csv"))
  cat("\nSaved results/summary_exp3_dropone.csv\n")
}

# ==============================================================================
# EXPERIMENT 4: k-FOLD CV (SPECIES AUC, TRAIN-TEST GAP, SITE LOG-LOSS)
# ==============================================================================

if (run_exp4) {

  cat("\n", strrep("=", 60), "\n")
  cat("EXPERIMENT 4: k-fold cross-validation\n")
  cat(strrep("=", 60), "\n")

  k = 10
  cat("Config:", rep_spatial, "alpha=", rep_alpha, "lambda=", rep_lambda, "\n")

  folds = build_stratified_folds(
    k = k, alt = altitude,
    tree_cover = Xenv[, "trees_cover"],
    shrub_cover = Xenv[, "shrubs_cover"],
    disturbance = Xenv[, "disturbance"],
    nutrient = Xenv[, "nutrient"],
    richness = richness
  )
  cat("Fold sizes:", paste(sapply(folds, length), collapse = ", "), "\n\n")

  # Out-of-fold predictions + per-fold train predictions
  P_oof = matrix(NA, nrow = nrow(Y), ncol = ncol(Y),
                 dimnames = list(rownames(Y), colnames(Y)))
  train_auc_per_fold = matrix(NA, nrow = k, ncol = ncol(Y))

  for (fold_i in seq_len(k)) {

    fold_file = here("Calanda_JSDM", "results", "cv", paste0("fold_", fold_i, ".rds"))
    test_idx = folds[[fold_i]]
    train_idx = setdiff(seq_len(nrow(Y)), test_idx)

    if (file.exists(fold_file)) {
      cat("[Fold", fold_i, "/", k, "] Cached\n")
      fold_data = readRDS(fold_file)
    } else {
      cat("[Fold", fold_i, "/", k, "] Fitting...\n")

      model_fold = fit_sjsdm(
        Y_train = Y[train_idx, ],
        Xenv_train = Xenv[train_idx, , drop = FALSE],
        XY_train = XY[train_idx, , drop = FALSE],
        spatial_form = rep_spatial,
        lambda_env = rep_lambda, alpha_env = rep_alpha,
        lambda_sp = rep_lambda, alpha_sp = rep_alpha,
        lambda_bio = rep_lambda, alpha_bio = rep_alpha
      )

      # Test predictions
      preds_test = predict(model_fold,
                           newdata = Xenv[test_idx, , drop = FALSE],
                           SP = XY[test_idx, , drop = FALSE])

      # Train predictions (for overfitting gap)
      preds_train = predict(model_fold,
                            newdata = Xenv[train_idx, , drop = FALSE],
                            SP = XY[train_idx, , drop = FALSE])

      fold_data = list(
        fold = fold_i,
        train_idx = train_idx,
        test_idx = test_idx,
        preds_test = preds_test,
        preds_train = preds_train
      )
      saveRDS(fold_data, fold_file)
      cat("  -> Saved\n")
    }

    P_oof[fold_data$test_idx, ] = fold_data$preds_test

    # Train AUC per species for this fold
    for (s in seq_len(ncol(Y))) {
      y_tr = Y[fold_data$train_idx, s]
      p_tr = fold_data$preds_train[, s]
      if (sum(y_tr == 1) >= 2 && sum(y_tr == 0) >= 2) {
        train_auc_per_fold[fold_i, s] = tryCatch(
          as.numeric(auc(roc(y_tr, p_tr, quiet = TRUE))),
          error = function(e) NA
        )
      }
    }
  }

  # --- Species metrics: test AUC + train-test gap ---
  species_cv = tibble(
    species = colnames(Y),
    n_presences = as.integer(colSums(Y)),
    test_auc = NA_real_,
    train_auc_mean = NA_real_,
    train_test_gap = NA_real_
  )

  for (s in seq_len(ncol(Y))) {
    y_true = Y[, s]
    y_pred = P_oof[, s]
    valid = !is.na(y_pred)

    if (sum(y_true[valid] == 1) >= 2 && sum(y_true[valid] == 0) >= 2) {
      species_cv$test_auc[s] = tryCatch(
        as.numeric(auc(roc(y_true[valid], y_pred[valid], quiet = TRUE))),
        error = function(e) NA
      )
    }

    species_cv$train_auc_mean[s] = mean(train_auc_per_fold[, s], na.rm = TRUE)
    species_cv$train_test_gap[s] = species_cv$train_auc_mean[s] - species_cv$test_auc[s]
  }

  species_cv = species_cv %>% arrange(test_auc)
  write_csv(species_cv, here("Calanda_JSDM", "results", "exp4_species_cv.csv"))
  cat("Saved results/exp4_species_cv.csv\n")

  # --- Site metrics: log-loss ---
  eps = 1e-7
  site_cv = tibble(
    site = rownames(Y),
    richness = as.integer(richness),
    logloss = NA_real_
  )

  for (i in seq_len(nrow(Y))) {
    y_true = Y[i, ]
    y_pred = pmin(pmax(P_oof[i, ], eps), 1 - eps)
    if (all(!is.na(y_pred))) {
      site_cv$logloss[i] = -mean(y_true * log(y_pred) + (1 - y_true) * log(1 - y_pred))
    }
  }

  site_cv = site_cv %>% arrange(desc(logloss))
  write_csv(site_cv, here("Calanda_JSDM", "results", "exp4_sites_cv.csv"))
  cat("Saved results/exp4_sites_cv.csv\n")

  cat("\nSpecies test AUC:\n")
  print(summary(species_cv$test_auc))
  cat("Train-test gap:\n")
  print(summary(species_cv$train_test_gap))
  cat("Site log-loss:\n")
  print(summary(site_cv$logloss))
}

# ==============================================================================
# EXPERIMENT 5: ANOVA SAMPLING SATURATION
# ==============================================================================

if (run_exp5) {

  cat("\n", strrep("=", 60), "\n")
  cat("EXPERIMENT 5: ANOVA sampling saturation\n")
  cat(strrep("=", 60), "\n")

  # Load representative model from Exp 1
  rep_run_id = paste0(rep_spatial, "_a", rep_alpha, "_l", rep_lambda)
  rep_file = here("Calanda_JSDM", "results", "runs", paste0(rep_run_id, ".rds"))

  if (file.exists(rep_file)) {
    rep_model = readRDS(rep_file)$model
  } else {
    stop("Representative run not found: ", rep_file, "\nRun Experiment 1 first.")
  }

  cat("Model:", rep_run_id, "\n\n")

  sample_sizes = c(100, 1000, 5000, 10000, 20000, 50000)
  saturation_results = list()

  for (ns in sample_sizes) {
    sat_file = here("Calanda_JSDM", "results", paste0("saturation_", ns, ".rds"))

    if (file.exists(sat_file)) {
      cat("Cached: anova_samples =", ns, "\n")
      partition = readRDS(sat_file)
    } else {
      cat("Computing: anova_samples =", ns, "...\n")
      partition = extract_partition(rep_model, anova_samples = ns)
      saveRDS(partition, sat_file)
    }

    saturation_results[[as.character(ns)]] = partition
  }

  # Compare all to 50000-sample reference
  ref = saturation_results[["50000"]]
  safe_cor = function(x, y) tryCatch(cor(x, y, use = "complete.obs"), error = function(e) NA)
  safe_mad = function(x, y) mean(abs(x - y), na.rm = TRUE)

  sat_table = tibble(anova_samples = sample_sizes)

  for (idx in seq_along(sample_sizes)) {
    p = saturation_results[[as.character(sample_sizes[idx])]]

    sat_table$cor_E_site[idx]    = safe_cor(p$sites$env, ref$sites$env)
    sat_table$cor_S_site[idx]    = safe_cor(p$sites$spa, ref$sites$spa)
    sat_table$cor_C_site[idx]    = safe_cor(p$sites$codist, ref$sites$codist)
    sat_table$cor_E_species[idx] = safe_cor(p$species$env, ref$species$env)
    sat_table$cor_S_species[idx] = safe_cor(p$species$spa, ref$species$spa)
    sat_table$cor_C_species[idx] = safe_cor(p$species$codist, ref$species$codist)

    sat_table$mad_E_site[idx]    = safe_mad(p$sites$env, ref$sites$env)
    sat_table$mad_S_site[idx]    = safe_mad(p$sites$spa, ref$sites$spa)
    sat_table$mad_C_site[idx]    = safe_mad(p$sites$codist, ref$sites$codist)
    sat_table$mad_E_species[idx] = safe_mad(p$species$env, ref$species$env)
    sat_table$mad_S_species[idx] = safe_mad(p$species$spa, ref$species$spa)
    sat_table$mad_C_species[idx] = safe_mad(p$species$codist, ref$species$codist)
  }

  write_csv(sat_table, here("Calanda_JSDM", "results", "exp5_saturation.csv"))
  cat("Saved results/exp5_saturation.csv\n")

  # Saturation plot
  p_sat = sat_table %>%
    select(anova_samples, starts_with("cor_")) %>%
    pivot_longer(-anova_samples, names_to = "component", values_to = "correlation") %>%
    ggplot(aes(x = anova_samples, y = correlation, color = component)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0.99, linetype = "dashed", color = "grey50") +
    scale_x_log10() +
    labs(title = "ANOVA sampling saturation (vs 50k reference)",
         x = "ANOVA samples", y = "Correlation with reference") +
    theme_bw() +
    theme(legend.position = "bottom")

  pdf(here("Calanda_JSDM", "results", "plots", "exp5_saturation.pdf"),
      height = 6, width = 8)
  print(p_sat)
  dev.off()
  cat("Saved results/plots/exp5_saturation.pdf\n")
}

# ==============================================================================
# SESSION INFO
# ==============================================================================

cat("\n=== All experiments complete ===\n")
print(session_info)
saveRDS(session_info, here("Calanda_JSDM", "results",
                           paste0("session_info_", Sys.Date(), ".rds")))
