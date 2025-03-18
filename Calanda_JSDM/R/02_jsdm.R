# Recheck the regularization stuff
# Increase the sampling
# Merge the shrub and tree cover, and integrate soil depth variation as a proxy of microclimate.
# Updated regularization parameters - maybe because now I increased I migh

lambda.env <- 0.001  
alpha.env <- 1.0

lambda.sp <- 0.002 
alpha.sp <- 0.2

lambda.bio <- 0.001
alpha.bio <- 1.0

learning_rate=0.01
sampling=500L 
device = "cpu"
iterations = 500L
act = "selu"

model = sjSDM(
  Y = Y,
  env = linear(X,
               formula = ~summer_temp + fdd + et.annual + slope + trees_cover + shrubs_cover + rocks_cover + soil_depth_mean + soil_depth_var, 
               lambda = lambda.sp, alpha = alpha.sp), 
  spatial = DNN(X %>% select(Latitude, Longitude), 
                formula = ~0+.,
                activation = act,
                hidden = rep(30, 2),
                bias = FALSE,
                lambda = lambda.sp, alpha = alpha.sp
  ),
  #spatial = linear(data = veg.env, formula = ~0+x*y, lambda = lambda.sp, alpha = alpha.sp),  
  biotic = bioticStruct(lambda = lambda.bio, alpha = alpha.bio, df = ncol(Y), reg_on_Cov = FALSE),
  iter = iterations,
  device = device,
  learning_rate = learning_rate,
  sampling = sampling,
  control = sjSDMControl(RMSprop(weight_decay = 0.0), 
                         scheduler = 5L, early_stopping_training = 25L, 
                         lr_reduce_factor = 0.9),
  se=T
)

saveRDS(model, file ="output/model_sjsdm_calanda.RDS")
plot(model$history, ylab = "Iteration losses", xlab = "Iteration") # maybe needs longer run
R2 = Rsquared(model, verbose = TRUE) #0.3496262
an = anova(model, verbose = TRUE, samples = sampling)
saveRDS(an, file ="output/an_sjsdm_calanda.RDS")
plot(an)
res = internalStructure(an, fractions = "proportional")
plot(res)
plotAssemblyEffects(res)

# AUC calculations
XY = veg.env[, c("x", "y")]
sklearn = reticulate::import("skmultilearn")
split = sklearn$model_selection$IterativeStratification(n_splits = 20L, order = 2L) 
obj = split$split(veg.env[,c("PC1", "PC2", "slope", "aspect")], Y)
splits = reticulate::iterate(obj)
CV = lapply(splits, function(sp) sp[[2]])
train_aucs = matrix(NA, length(CV), ncol(Y))
test_predictions = list()

for(i in 1:length(CV)) {
  indices_test = CV[[i]]
  modelt = sjSDM(
    Y = Y[-indices_test,],
    env = linear(veg.env[-indices_test,],
                 formula = ~PC1 + PC2 + slope + aspect, 
                 lambda = lambda.sp, alpha = alpha.sp), 
    spatial = linear(data = XY[-indices_test,], formula = ~0+x*y, lambda = lambda.sp, alpha = alpha.sp),  
    biotic = bioticStruct(lambda = lambda.bio, alpha = alpha.bio, df = ncol(Y), reg_on_Cov = FALSE),
    iter = iterations,
    learning_rate = learning_rate,
    sampling = sampling,
    device = device,
    control = sjSDMControl(RMSprop(weight_decay = 0.0), 
                           scheduler = 5L, early_stopping_training = 25L, 
                           lr_reduce_factor = 0.9)
  )
  train_pred = predict(modelt)
  test_predictions = append(test_predictions, list(predict(modelt, newdata = veg.env[indices_test,], SP = XY[indices_test,])))
  train_aucs[i,] = sapply(1:ncol(train_pred), function(i) Metrics::auc(Y[-indices_test,i], train_pred[,i]))
}
pred = do.call(rbind, test_predictions)
test_obs = do.call(rbind, lapply(CV, function(ind) Y[ind,]))
train_aucs = apply(train_aucs, 2, mean)
test_aucs = sapply(1:ncol(test_obs), function(i) Metrics::auc(test_obs[,i], pred[,i]))
save(train_aucs, test_aucs ,file = "output/model_AUC_calanda.rdata")

spnames = colnames(Y)
auc = data.frame(sp=spnames,e.aucs=train_aucs,p.aucs=test_aucs)  
aucplot(auc) 

# Get TP species colonizer/extinct
wl.determined = readRDS("~/Desktop/TransPlant_winners_losers/output/wl.determined.rds") %>%
  filter(Region == "CH_Calanda")%>%
  dplyr::select(Region, originSiteID, destSiteID, SpeciesName, wl)%>%
  distinct()%>%
  rename(taxon = SpeciesName)%>%
  left_join(cal_pools%>% dplyr::select(Region, originSiteID, destSiteID, taxon, pool))%>%
  mutate(cat = ifelse(!is.na(pool), paste0(wl, "_", pool), wl))

# Build ternary plots
x=an
type="R2_McFadden"
#x$species$R2_McFadden_shared
#modify by the original plot.anova from sjsdm package
internals = list()
df = data.frame(
  env = ifelse(x$sites[[type]]$F_A<0, 0, x$sites[[type]]$F_A),
  spa = ifelse(x$sites[[type]]$F_S<0, 0, x$sites[[type]]$F_S),
  codist = ifelse(x$sites[[type]]$F_B<0, 0, x$sites[[type]]$F_B),
  r2  = ifelse(x$sites[[type]]$Full<0, 0, x$sites[[type]]$Full)
)
internals[[1]] = df
names(internals)[1] = "Sites"

df = data.frame(
  env = ifelse(x$species[[type]]$F_A<0, 0, x$species[[type]]$F_A),
  spa = ifelse(x$species[[type]]$F_S<0, 0, x$species[[type]]$F_S),
  codist = ifelse(x$species[[type]]$F_B<0, 0, x$species[[type]]$F_B),
  r2  = ifelse(x$species[[type]]$Full<0, 0, x$species[[type]]$Full),
  group=NA,
  taxon=model$species)%>% 
  left_join(wl.determined)%>%
  mutate(cat = as.character(cat))%>%
  mutate(cat = ifelse(is.na(cat), "not_tp", cat))

internals[[2]] <- df
names(internals)[2] = "Species"

df %>%
  dplyr::select(env:r2, cat)%>%
  pivot_longer(cols = env:r2)%>%
  mutate(cat = factor(cat, levels = c("loser", "loser_extinct", "winner", "winner_colonizing", "not_tp")))%>%
  ggplot(aes(name, value, color = cat))+
  theme_bw()+
  geom_boxplot()+
  # scale_y_continuous(limits=c(0,1),
  #                    breaks=seq(0, 1,by=0.2),
  #                    labels=seq(0, 1,by=0.2))+
  scale_color_manual(values = c("winner_colonizing" = "#FE6100", "winner"="#DC267F", "loser" = "#648FFF", "loser_extinct" = "#785EF0","not_tp" = "grey70"))+
  labs(x = "", y = "R2 (Species)")

ggtern::ggtern(internals[[2]], ggplot2::aes(x = env, z = spa, y = codist, color = cat)) +
  ggtern::scale_T_continuous(limits=c(0,1),
                             breaks=seq(0, 1,by=0.2),
                             labels=seq(0,1, by= 0.2)) +
  ggtern::scale_L_continuous(limits=c(0,1),
                             breaks=seq(0, 1,by=0.2),
                             labels=seq(0, 1,by=0.2)) +
  ggtern::scale_R_continuous(limits=c(0,1),
                             breaks=seq(0, 1,by=0.2),
                             labels=seq(0, 1,by=0.2)) +
  ggplot2::labs(title = names(internals)[2],
                x = "E",
                xarrow = "Environment",
                y = "C",
                yarrow = "Co-Distribution",
                z = "S", 
                zarrow = "Space") +
  ggtern::theme_bw() +
  ggtern::theme_showarrows() +
  ggtern::theme_arrowlong() +
  ggplot2::theme(
    panel.grid = ggplot2::element_line(color = "darkgrey", linewidth = 0.3),
    plot.tag = ggplot2::element_text(size = 11),
    plot.title = ggplot2::element_text(size = 11, hjust = 0.1 , margin = ggplot2::margin(t = 10, b = -20)),
    tern.axis.arrow = ggplot2::element_line(linewidth = 1),
    tern.axis.arrow.text = ggplot2::element_text(size = 6),
    axis.text = ggplot2::element_text(size = 4),
    axis.title = ggplot2::element_text(size = 6),
    legend.text = ggplot2::element_text(size = 15),
    legend.title = ggplot2::element_text(size = 8),
    strip.text = ggplot2::element_text(size = 15),
    #plot.margin = unit(c(top,1,1,1)*0.2, "cm"),
    strip.background = ggplot2::element_rect(color = NA),
  ) +
  ggplot2::theme(tern.axis.arrow.text = element_text(size = 15),legend.position = "right", legend.margin = margin(r = 5), legend.box="vertical")+
  ggplot2::guides(size = ggplot2::guide_legend(title = expression(R^2), order = 1, nrow = 1, label.position = "bottom"))+
  scale_color_manual(values = c("winner_colonizing" = "#FE6100", "winner"="#DC267F", "loser" = "#648FFF", "loser_extinct" = "#785EF0","not_tp" = "grey70"))+
  geom_point(alpha = 0.4, size = 3, position= position_jitter_tern(x=0.01, y=0.01, z=0.00))+
  labs(color = "")


# Assuming your data frame is named 'df'
df_long <- df %>%
  dplyr::select(env:r2, cat) %>%                       # Select relevant columns
  pivot_longer(cols = env:r2,                           # Pivot from wide to long format
               names_to = "name",
               values_to = "value") %>%
  mutate(cat = factor(cat,                            # Define the order of categories
                      levels = c("loser", 
                                 "loser_extinct", 
                                 "winner", 
                                 "winner_colonizing", 
                                 "not_tp")))
# Fit the linear model
# Including 'name' allows for different effects across 'name' categories
model <- lm(value ~ cat * name, data = df_long %>% filter(cat != "not_tp"))
# Obtain estimated marginal means for 'cat' within each 'name'
emm <- emmeans(model, ~ cat | name)
# Perform pairwise comparisons with adjustment for multiple testing
pairwise_diff <- pairs(emm, adjust = "mvt")  # You can use "tukey" or other methods
# View the summary of pairwise differences
summary(pairwise_diff)



# Try gjam
modelList <- list(seed=123, nIter=2000, burnin=500, thin=10, typeNames = "FC")
f = as.formula( ~ MAT + AP + x + y )
gjam_fit <- gjam(f, xdata = as.matrix(veg.env), ydata = as.matrix(veg.abund), modelList=modelList)
heatmap(gjam_fit$parameters$sigMu)
heatmap(gjam_fit$parameters$sigMu, Colv = NA, Rowv = NA, scale="column", col = cm.colors(256))
plotPars  <- list(GRIDPLOTS=T, CLUSTERPLOTS=T,SAVEPLOTS = F)
gjamPlot(gjam_fit, plotPars)