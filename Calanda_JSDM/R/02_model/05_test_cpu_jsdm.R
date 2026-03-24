# ==============================================================================
# Script: 05_test_cpu_jsdm.R
# Author: Billur Bektas
# Claude (Anthropic) was used to assist with code refactoring, validation, and documentation.
#
# Purpose: Local CPU test of the sjSDM model structure (low iterations).
#          Verifies model specification before running full experiments on Colab GPU.
#          Not used for production results.
#
# Inputs:
#   - output/data_calanda_jsdm_<date>.rds (X and Y matrices from 04_prepare)
#
# Outputs:
#   - None (test only; full model fitting runs on Google Colab)
#
# Note: Full model fitting + experiments run via Colab notebook:
#       R/02_model/06_sjsdm_experiments_colab.ipynb
# ==============================================================================

library(sjSDM)
library(tidyverse)
library(conflicted)
library(here)

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# Load data ----
date.data = "2026-02-25"
data_calanda_jsdm = readRDS(here("Calanda_JSDM", "output", paste0("data_calanda_jsdm_",Sys.Date(),".rds")))
X = data_calanda_jsdm$X
Y = data_calanda_jsdm$Y

# Hyperparameters ----
lambda_env = 0.001
alpha_env = 1.0

lambda_sp = 0.002
alpha_sp = 0.2

lambda_bio = 0.001
alpha_bio = 1.0

learning_rate = 0.01
sampling = 100L
device = "cpu"
iterations = 50L
act = "selu"

# Fit model ----
cat("\n=== Fitting sjSDM model ===\n")

model = sjSDM(
  Y = Y,
  env = linear(X,
    formula = ~summer_temp + fdd + et.annual + slope + rocks_cover +
      trees_cover + shrubs_cover + soil_depth_var +
      tpi + flowdir + nutrient + disturbance,
    lambda = lambda_sp, alpha = alpha_sp),
  spatial = DNN(X %>% select(Latitude, Longitude),
    formula = ~0 + .,
    activation = act,
    hidden = rep(30, 2),
    bias = FALSE,
    lambda = lambda_sp, alpha = alpha_sp),
  biotic = bioticStruct(lambda = lambda_bio, alpha = alpha_bio,
                        df = ncol(Y), reg_on_Cov = FALSE),
  iter = iterations,
  device = device,
  learning_rate = learning_rate,
  sampling = sampling,
  control = sjSDMControl(RMSprop(weight_decay = 0.0),
                         scheduler = 5L,
                         early_stopping_training = 25L,
                         lr_reduce_factor = 0.9),
  se = FALSE
)

saveRDS(model, here("Calanda_JSDM", "output", paste0("model_sjsdm_calanda_", Sys.Date(), ".rds")))

# Variance partitioning ----
cat("Computing R², ANOVA, and internal structure...\n")

R2 = Rsquared(model, verbose = TRUE)
an = anova(model, verbose = TRUE, samples = sampling)
res = internalStructure(an, fractions = "proportional")

saveRDS(R2, here("Calanda_JSDM", "output", paste0("R2_sjsdm_calanda_", Sys.Date(), ".rds")))
saveRDS(an, here("Calanda_JSDM", "output", paste0("an_sjsdm_calanda_", Sys.Date(), ".rds")))
saveRDS(res, here("Calanda_JSDM", "output", paste0("res_sjsdm_calanda_", Sys.Date(), ".rds")))

cat("\n=== Model fitting complete ===\n")
