library(sjSDM)
library(tidyverse)
library(conflicted)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

data_calanda_jsdm = readRDS("output/data_calanda_jsdm.rds")
X = data_calanda_jsdm$X
Y = data_calanda_jsdm$Y

lambda.env = 0.001  
alpha.env = 1.0

lambda.sp = 0.002 
alpha.sp = 0.2

lambda.bio = 0.001
alpha.bio = 1.0

learning_rate = 0.01
sampling = 50000L 
device = "gpu"
iterations = 750L 
act = "selu"

model = 
  sjSDM(
  Y = Y,
  env = linear(X,
               formula = ~summer_temp + fdd + et.annual + slope + rocks_cover + trees_cover + shrubs_cover + soil_depth_mean + soil_depth_var + tpi + flowdir + roughness + land_use, 
               lambda = lambda.sp, alpha = alpha.sp), 
  spatial = DNN(X %>% select(Latitude, Longitude), 
                formula = ~0+.,
                activation = act,
                hidden = rep(30, 2),
                bias = FALSE,
                lambda = lambda.sp, alpha = alpha.sp),
  #spatial = linear(data = veg.env, formula = ~0+x*y, lambda = lambda.sp, alpha = alpha.sp),  
  biotic = bioticStruct(lambda = lambda.bio, alpha = alpha.bio, df = ncol(Y), reg_on_Cov = FALSE),
  iter = iterations,
  device = device,
  learning_rate = learning_rate,
  sampling = sampling,
  control = sjSDMControl(RMSprop(weight_decay = 0.0), 
                         scheduler = 5L, 
                         early_stopping_training = 25L, 
                         lr_reduce_factor = 0.9),
  se=T
)

saveRDS(model, file ="output/model_sjsdm_calanda_260425.RDS")

R2 = Rsquared(model, verbose = TRUE) #0.3496262
an = anova(model, verbose = TRUE, samples = sampling)
res = internalStructure(an, fractions = "proportional")

saveRDS(R2, file ="output/R2_sjsdm_calanda_260425.RDS")
saveRDS(an, file ="output/an_sjsdm_calanda_260425.RDS")
saveRDS(res, file ="output/res_sjsdm_calanda_260425.RDS")

