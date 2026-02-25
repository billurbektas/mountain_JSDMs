# ==============================================================================
# Script: 02_jsdm.R
# Purpose: Fit spatial Joint Species Distribution Model using sjSDM (GPU)
#
# Inputs:
#   - output/data_calanda_jsdm.rds (X and Y matrices from 01_prepare_data.R)
#
# Outputs:
#   - output/model_sjsdm_calanda.rds
#   - output/R2_sjsdm_calanda.rds
#   - output/an_sjsdm_calanda.rds
#   - output/res_sjsdm_calanda.rds
#
# Note: Requires GPU. Authoritative results are in results_from_Max/.
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
