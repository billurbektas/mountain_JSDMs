library(sjSDM)
library(tidyverse)
library(conflicted)

conflict_prefer("select", "dplyr")

data_calanda_jsdm = readRDS("data_calanda_jsdm.rds")
X = data_calanda_jsdm$X
Y = data_calanda_jsdm$Y

lambda.env = 0.001  
alpha.env = 1.0

lambda.sp = 0.002 
alpha.sp = 0.2

lambda.bio = 0.001
alpha.bio = 1.0

learning_rate = 0.01
sampling = 20000L 
device = "gpu"
iterations = 750L 
act = "selu"

# remove species with less than 5 occurrences (not enough information!)
Y_sub = Y[,colSums(Y) > 4]

model = 
  sjSDM(
    Y = Y_sub,
    env = linear(X,
                 formula = ~summer_temp + fdd + et.annual + slope + rocks_cover + trees_cover + shrubs_cover + soil_depth_mean + soil_depth_var + tpi + flowdir + land_use,
                 lambda = lambda.sp, alpha = alpha.sp), 
    spatial = DNN(X %>% select(Latitude, Longitude), 
                  formula = ~0+.,
                  activation = act,
                  hidden = rep(30, 2),
                  bias = FALSE,
                  lambda = lambda.sp, alpha = alpha.sp),
    #spatial = linear(data = veg.env, formula = ~0+x*y, lambda = lambda.sp, alpha = alpha.sp),  
    biotic = bioticStruct(lambda = lambda.bio, alpha = alpha.bio, df = ncol(Y), reg_on_Cov = TRUE),
    iter = iterations,
    device = device,
    learning_rate = learning_rate,
    sampling = sampling,
    step_size = 100L,
    control = sjSDMControl(RMSprop(weight_decay = 0.001), 
                           scheduler = 30L, 
                           early_stopping_training = 50L, 
                           lr_reduce_factor = 0.9),
    se=T
  )
gc()
sjSDM:::pkg.env$torch$cuda$empty_cache()

saveRDS(model, file ="model_sjsdm_calanda.rds")

R2 = Rsquared(model, verbose = TRUE) #0.3496262

gc()
sjSDM:::pkg.env$torch$cuda$empty_cache()

an = anova(model, verbose = TRUE, samples = sampling)
res = internalStructure(an, fractions = "proportional")

saveRDS(R2, file = "R2_sjsdm_calanda.rds")
saveRDS(an, file = "an_sjsdm_calanda.rds")
saveRDS(res, file = "res_sjsdm_calanda.rds")

