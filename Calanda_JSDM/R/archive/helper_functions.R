bioclim_data = data.frame(
  Variable = c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", 
               "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", 
               "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", "AI", "etpY","gdd0Y", "gdd3Y", "gdd5Y"),
  Description = c("Annual Mean Temperature",
                  "Mean Diurnal Range (Mean of monthly (max temp - min temp))",
                  "Isothermality",
                  "Temperature Seasonality (standard deviation ×100)",
                  "Max Temperature of Warmest Month",
                  "Min Temperature of Coldest Month",
                  "Temperature Annual Range",
                  "Mean Temperature of Wettest Quarter",
                  "Mean Temperature of Driest Quarter",
                  "Mean Temperature of Warmest Quarter",
                  "Mean Temperature of Coldest Quarter",
                  "Annual Precipitation",
                  "Precipitation of Wettest Month",
                  "Precipitation of Driest Month",
                  "Precipitation Seasonality",
                  "Precipitation of Wettest Quarter",
                  "Precipitation of Driest Quarter",
                  "Precipitation of Warmest Quarter",
                  "Precipitation of Coldest Quarter",
                  "Aridity Index", 
                  "Evapotranspiration",
                  "Sum of growing degree days (above 0 degrees C)",
                  "Sum of growing degree days (above 3 degrees C)",
                  "Sum of growing degree days (above 5 degrees C)"))
match_clim_var = function(var){
bioclim_data = data.frame(
  Variable = c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", 
               "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", 
               "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", "AI", "etpY","gdd0Y", "gdd3Y", "gdd5Y"),
  Description = c("Annual Mean Temperature",
                  "Mean Diurnal Range (Mean of monthly (max temp - min temp))",
                  "Isothermality",
                  "Temperature Seasonality (standard deviation ×100)",
                  "Max Temperature of Warmest Month",
                  "Min Temperature of Coldest Month",
                  "Temperature Annual Range",
                  "Mean Temperature of Wettest Quarter",
                  "Mean Temperature of Driest Quarter",
                  "Mean Temperature of Warmest Quarter",
                  "Mean Temperature of Coldest Quarter",
                  "Annual Precipitation",
                  "Precipitation of Wettest Month",
                  "Precipitation of Driest Month",
                  "Precipitation Seasonality",
                  "Precipitation of Wettest Quarter",
                  "Precipitation of Driest Quarter",
                  "Precipitation of Warmest Quarter",
                  "Precipitation of Coldest Quarter",
                  "Aridity Index", 
                  "Evapotranspiration",
                  "Sum of growing degree days (above 0 degrees C)",
                  "Sum of growing degree days (above 3 degrees C)",
                  "Sum of growing degree days (above 5 degrees C)"))
  return(bioclim_data$Description[bioclim_data$Variable == var])
}

aucplot =  function(auc) {
  mean.e=signif(mean(auc$e.aucs,na.rm = T),digits = 2)  
  mean.p=signif(mean(auc$p.aucs,na.rm = T),digits = 2) 
  ggplot(auc,aes(x =e.aucs, y =p.aucs))+
    scale_color_viridis(option="D",discrete=TRUE)+
    xlab("explanatory AUC") + 
    ylab("predictive AUC") + 
    xlim(0.3,1)+
    ylim(0.3,1)+
    geom_abline(intercept=0,slope=1,colour='azure3')+
    geom_hline(yintercept = 0.5,colour='azure3')+
    geom_vline(xintercept = 0.5,colour='azure3')+
    geom_hline(yintercept = mean.p,linetype="dotted")+
    geom_vline(xintercept = mean.e,linetype="dotted")+
    labs(fill="Species Index") + 
    geom_point()+ 
    ggplot2::annotate("text",label=mean.p, x = 0.45, y = mean.p+0.02)+
    ggplot2::annotate("text",label=mean.e, x = 0.75, y = 0.4)+
    theme(panel.grid.major =element_blank(), panel.grid.minor  =element_blank(),panel.background = element_blank(),panel.border = element_blank())+
    geom_text_repel(aes(e.aucs, p.aucs, label=sp,size=20,fontface = "italic"),show.legend = FALSE)
}

