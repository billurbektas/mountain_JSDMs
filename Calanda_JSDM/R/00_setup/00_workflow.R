library(sjSDM)
library(tidyverse)
library(terra)
library(TNRS)
library(stringi)
library(sf)
library(ggpubr)
library(FactoMineR)
library(factoextra)
library(corrplot)
library(gt)
library(gridExtra)
library(grid)
library(conflicted)
library(miceRanger)
library(patchwork)
library(ggrepel)
library(spatialEco)
library(randomForest)
library(caret)
library(mvabund)
library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
#library(ggtern) ## doesnt work. 
library(patchwork)
library(factoextra)
library(FactoMineR)
library(qgam)
library(mvabund)
library(here)
library(ggforce)
library(ggrepel)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("aes", "ggplot2")
conflict_prefer("extract", "terra")
conflict_prefer("year", "lubridate")
conflict_prefer("first", "dplyr")
conflict_prefer("intersect", "base")
conflict_prefer("theme_minimal", "ggplot2")
conflict_prefer("theme_bw", "ggplot2")
conflict_prefer("shift", "terra")
conflict_prefer("margin", "ggplot2")
conflicts_prefer(ggplot2::theme_classic)

setwd("Calanda_JSDM/")
source("R/functions_calanda.R")

# Load data
load("data/vegetation/TransPlantNetwork_101024.RData")
load("output/starter_data_25.04.25.RData")
#https://github.com/TheoreticalEcology/s-jSDM
#https://cran.r-project.org/web/packages/sjSDM/vignettes/sjSDM_Introduction.html

# Unset any existing Python configurations
# Sys.unsetenv("RETICULATE_PYTHON")
# 
# # Load libraries in order
# library(reticulate)
# use_condaenv("r-sjsdm", required = TRUE)
# library(sjSDM)

# Protocol of the CapHE data: "/Volumes/green_groups_PLEC_public/Research/Calanda/2024_CAPHE/CAPHE_SpeciesDistribution/CAPHE_SpeDis_Documents/CAPHE_SpeDis_Methods/2023_CAPHE_Methods_20240207.docx"

# reticulate::py_config()
# reticulate::py_install("scikit-multilearn", pip = TRUE)
# sklearn = reticulate::import("skmultilearn")


veg_coord = read_csv("data/vegetation/2024_CAPHE_SpeDis_CleanData_20240214.csv")%>%
  select(plot_id, releve_id, x, y )%>%
  distinct()
write_csv(veg_coord, file = "data/vegetation/veg.coord.csv")
# Get coordinates data for downloads
veg_coord = read_csv("data/vegetation/veg.coord.csv")[,-1]

veg_coord_ecostress = 
  veg_coord %>%
  mutate(ID = paste0(plot_id, "_", releve_id))%>%
  rename(Longitude = x, Latitude = y)%>%
  dplyr::select(ID, Latitude, Longitude)
write.csv(veg_coord_ecostress, file = "data/veg_coord_ecostress.csv")

# Get functions
source("R/functions_calanda.R")


# Workflow
if(file.exists("output/starter_data_25.04.25.RData")){
  load("output/starter_data_25.04.25.RData")
}else{source("R/01_prepare_data.R")}

