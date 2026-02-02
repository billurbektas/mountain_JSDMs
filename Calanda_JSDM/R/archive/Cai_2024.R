#devtools::install_github("https://github.com/TheoreticalEcology/s-jSDM", subdir = "sjSDM", ref = "master")
library(sjSDM)
library(tidyverse)
library(conflicted)
set.seed(42)
conflict_prefer("select", "dplyr")

waterspecies="FALSE"  #TRUE or FALSE
learning_rate=0.01
sampling=50000L #increase stability
act = "selu"
device = 0
#######load data#########
load("FHT_320_48sp_20211110.rdata")

if (waterspecies=="TRUE"){
  Y=data.frame(Y)%>%select(starts_with("Acti") | starts_with("Amph"))
  ex=which(rowSums(Y) == 0)
  Y=Y[-ex,]
  env=env[-ex,]
  sp=sp[-ex,]
}
Y=as.matrix(Y)
env=scale(env)
XY=scale(sp)

############### Model ################
lambda.env = 0.001
alpha.env = 1.0
lambda.sp = 0.002
alpha.sp = 0.2
lambda.bio = 0.001
alpha.bio = 1.0

model = sjSDM(
  Y = Y,
  env = linear(env,
               formula = ~HSI2_Pond_area+ HSI5_Shade+ HSI10_Macrophytes+ HSI3_Pond_drying+ HSI4_Water_quality+ agriculture.urban+ woodland+ grassland, 
               lambda = lambda.sp, alpha = alpha.sp), 
  spatial = DNN(XY, 
                formula = ~0+.,
                activation = act,
                hidden = rep(30, 2),
                bias = FALSE,
                lambda = lambda.sp, alpha = alpha.sp
  ),  
  biotic = bioticStruct(lambda = lambda.bio, alpha = alpha.bio, df = ncol(Y), reg_on_Cov = FALSE),
  iter = 500L,
  learning_rate = learning_rate,
  sampling=sampling,
  device = device,
  control = sjSDMControl(RMSprop(weight_decay = 0.0), 
                         scheduler = 5L, early_stopping_training = 25L, 
                         lr_reduce_factor = 0.9),
  se=T
)
plot(model$history)
an=anova(model, samples = sampling) 

saveRDS(model,file ="model_sjsdm_12042023.RDS")
saveRDS(an,file ="model_an_sjsdm_12042023.RDS")

R2 = Rsquared(model)
summary.p<-summary(model)
co.env.spe <- cov2cor(getCov(model))

save(R2,summary.p,co.env.spe,file = "model_48sp_result_12042023.rdata")
rm(model,an,R2,summary.p,co.env.spe)
######### CV + AUC calculation ###########
sklearn = reticulate::import("skmultilearn")
split = sklearn$model_selection$IterativeStratification(n_splits = 20L,order = 2L) 
obj = split$split(veg.env[,c("PC1", "PC2", "slope", "aspect")], Y)
splits = reticulate::iterate(obj)
CV = lapply(splits, function(sp) sp[[2]])

train_aucs = matrix(NA, length(CV), ncol(Y))
test_predictions = list()

for(i in 1:length(CV)) {
  indices_test = CV[[i]]
  model = sjSDM(
    Y = Y[-indices_test,],
    env = linear(veg.env[-indices_test,],
                 formula = ~HSI2_Pond_area+ HSI5_Shade+ HSI10_Macrophytes+ HSI3_Pond_drying+ HSI4_Water_quality+ agriculture.urban+ woodland+ grassland, 
                 lambda = lambda.sp, alpha = alpha.sp), 
    spatial = DNN(XY[-indices_test,], 
                  formula = ~0+.,
                  hidden = rep(30, 2),
                  bias = FALSE,
                  activation = act,
                  lambda = lambda.sp, alpha = alpha.sp
    ),  
    biotic = bioticStruct(lambda = lambda.bio, alpha = alpha.bio, df = ncol(Y), reg_on_Cov = FALSE),
    iter = 500L,
    learning_rate = learning_rate,
    sampling=sampling,
    device = device,
    control = sjSDMControl(RMSprop(weight_decay = 0.0), 
                           scheduler = 5L, early_stopping_training = 25L, 
                           lr_reduce_factor = 0.9)
  )
  train_pred = predict(model)
  test_predictions = append(test_predictions, list(predict(model, newdata = env[indices_test,], SP = XY[indices_test,])))
  train_aucs[i,] = sapply(1:ncol(train_pred), function(i) Metrics::auc(Y[-indices_test,i], train_pred[,i]))
}
pred = do.call(rbind, test_predictions)
test_obs = do.call(rbind, lapply(CV, function(ind) Y[ind,]))
train_aucs = apply(train_aucs, 2, mean)
test_aucs = sapply(1:ncol(test_obs), function(i) Metrics::auc(test_obs[,i], pred[,i]))
save(train_aucs, test_aucs ,file = "model_48sp_AUC_12042023.rdata")