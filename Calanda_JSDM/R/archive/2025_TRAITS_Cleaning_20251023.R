################################
## C&N ANALYSIS DATA CLEANING ##
################################

## 1. Preparation ====
setwd("/Volumes/green_groups_plec_public/Research/Calanda/Calanda_Projects/2024_CAPHE/CAPHE_Traits")
## set working directory
library(tidyr)
library(dplyr)
library(lubridate)
library(tidyverse)
library(sf)

  ## 1.1. Read in file
  CN_total <- read.csv("/2024_TRAITS_OriginalData/2025_CN_analysis/CorrectedData/CSV/2025_TRAITS_CN_total_20250908.csv", header = TRUE, sep = ";")
    ## !! IMPORTANT: combined CN data from files Trial 1, Trial 2, Sample 1, Sample 2, Sample 3 - NEEDS TO BE UPDATED WITH ADDITIONAL SAMPLES!

  ## 1.2. Change variable class
  CN_total$N_content <- as.numeric(CN_total$N_content)
  CN_total$N_corr_factor <- as.numeric(CN_total$N_corr_factor)
  CN_total$date <- as.Date(CN_total$date, "%d.%m.%y")
  
## 2. Delete unnecessary units ====
  
  ## 2.1. Remove reference rows (Ali, bl, Caf, Tyr)
  CN_total <- CN_total[!grepl("Ali|bl|Caf|Tyr", CN_total$sample_id), ]
  
  ## 2.2. Remove unnecessary rows
  CN_colnames <- colnames(CN_total)
    ## get vector with colnames
  CN_red_col <- CN_total[, CN_colnames[c(1, 3:5, 7:8, 10, 12:13, 15, 17:18, 20)]] 
    ## select relevant columns
  
## 3. Long to wide data format ====
  
  ## 3.1. Split dataframe according to variable in "sheet"
  C_CN <- CN_red_col %>% filter(sheet == "C") %>% select(-N_content, -N_corr_factor, -N_content_corr)
    ## all columns apart from the N ones
  N_CN <- CN_red_col %>% filter(sheet == "N") %>% select(analysis_id, N_content, N_corr_factor, N_content_corr)
    ## just N variables + the one to merge the sheet
  
  ## 3.2. Merge files again
  CN_joint <- left_join(C_CN, N_CN, by = "analysis_id")
  
  ## 3.3. Adjust columns
  CN_joint <- CN_joint[, c(1:5, 7:9, 11:13)]
    ## remove "sheet" and comment
  
## 4. Delete duplicates ====
  
  ## 4.1. Delete duplicates
  CN_clean <- CN_joint %>% 
    group_by(sample_id) %>%         # Group by CustomerID
    filter(date == max(date)) %>%    # Filter date to return max date of each group
    ungroup()                        # Ungroup
  ## keep most recent date
  
  ## 4.2. Order column sample id
  CN_clean$sample_id <- as.numeric(CN_clean$sample_id)
  CN_clean <- CN_clean[order(CN_clean$sample_id),]
  CN_clean <- CN_clean %>% select(-file)
    ## delete column "file" (not used anymore)
  
  ## 4.3. Add comments
  comments <- CN_total %>% select(sample_id, comment)
    ## get sample ids and comments
  comments$sample_id <- as.numeric(comments$sample_id)
    ## same class as in overall df
  comments[comments == ""] <- NA
    ## set blank cells to NA
  comments_index <- which(!is.na(comments$comment))
    ## get indices from the cells with content
  comments <- comments[comments_index,]
    ## get a small df with sample_ids and comments
  CN_clean_comm <- left_join(CN_clean, comments, by = "sample_id")
  
## 5. Match sample ids with plant id (scans) ====
  
  ## 5.1. Import file with names
  sample_id <- read.csv("2025_CAPHE_traits_sample_id.csv", header = TRUE, sep = ",")
    ## import document that matches sample ids with plant ids
  sample_id <- sample_id[, c(1:2)]
    ## delete unnecessary columns
  colnames(sample_id) <- c("sample_id", "plant_id")
    ## rename columns for merging
  sample_id$sample_id <- as.numeric(sample_id$sample_id)
  CN_id <- left_join(CN_clean_comm, sample_id, by = "sample_id")
    ## merge files
  CN_final <- CN_id[, c(12, 1, 4:11)]
    ## rearrange colum order
  
## 6. Export file ====
write.csv(CN_final, file = "/2024_TRAITS_CleanData/2025_TRAITS_CN_clean_20251023.csv", row.names = FALSE)
  
  
##############################
## MERGE C&N AND TRAIT FILE ##
##############################
  
## 1. Preparations====
  
  setwd("~/Desktop/CN_Data_Cleaning/OriginalDocuments")
  ## set working directory  
  
  ## 1.1. Load packages
  library(tidyr)
  library(dplyr)
  library(lubridate)
  library(tidyverse)
  library(sf)
  library(rstatix)
  
  ## 1.2. Read in files
  data_CN <- read.csv("/2024_TRAITS_CleanData/2025_TRAITS_CN_clean_20251023.csv", header = TRUE, sep = ",")
  data_traits <-  read.csv("/2024_TRAITS_OriginalData/2024_Traits/CorrectedData/CSV/2024_TRAITS_Data_total.csv", header = TRUE, sep = ";")
  
## 2. Cleaning====
  
  ## 2.1. Delete redundant columns
  data_traits <- data_traits[, -c(8, 15:17)]
  
  ## 2.2. Change class of columns
  data_traits$plant_species <- as.factor(data_traits$plant_species)
  data_traits$site <- as.factor(data_traits$site) 
  data_traits$collector <- as.factor(data_traits$collector)
  data_traits$site <- as.factor(data_traits$site) 
  data_traits$collector <- as.factor(data_traits$collector) 
  data_traits$vegetative_height <- as.integer(data_traits$vegetative_height)
  data_traits$reproductive_height <- as.integer(data_traits$reproductive_height)
  data_traits$diameter_long <- as.integer(data_traits$diameter_long)
  data_traits$diameter_90 <- as.integer(data_traits$diameter_90)
  data_traits$vegetative_height_stretched <- as.integer(data_traits$vegetative_height_stretched)
  data_traits$reproductive_height_stretched <- as.integer(data_traits$reproductive_height_stretched)
  data_traits$flowering_unit <- as.factor(data_traits$flowering_unit)
  data_traits$phenology_state <- as.factor(data_traits$phenology_state)
  data_traits$flowering_unit_count <- as.integer(data_traits$flowering_unit_count)
  data_traits$inflorescence_count <- as.integer(data_traits$inflorescence_count)
  data_traits$grass_spikelet_count <- as.integer(data_traits$grass_spikelet_count)
  data_traits$leaves_sampled <- as.integer(data_traits$leaves_sampled)
  data_traits$leaves_other <- as.factor(data_traits$leaves_other)
  data_traits$flowering_unit_fresh_weight <- as.numeric(data_traits$flowering_unit_fresh_weight)
  data_traits$leaves_scans <- as.factor(data_traits$leaves_scans)
  
  ## 2.3. Clean column by column
  
    ## 2.3.1. Plant species
    levels(data_traits$plant_species)
      ## displays categories in column plant species
    data_traits[(data_traits$plant_species == ""), ] 
    which(data_traits$plant_species == "") 
      ## empty rows (delete them)
    data_traits <- data_traits[-c(1821:1825),]
      ## delete empty rows
    data_traits[(data_traits$plant_species == "Agrostis capillaris "), "plant_species"] <- "Agrostis capillaris"
    ## clear spelling mistake 
    data_traits[(data_traits$plant_species == "Anthoxanthum odoratum "), "plant_species"] <- "Anthoxanthum odoratum"
      ## clear spelling mistake
    data_traits[(data_traits$plant_species == "Brachipodium pinatum"), "plant_species"] <- "Brachypodium pinnatum"
      ## clear spelling mistake 
    data_traits[(data_traits$plant_species == "Carex caryophyllea "), "plant_species"] <- "Carex caryophyllea"
      ## clear spelling mistake
    data_traits[(data_traits$plant_species == "Cirsium acaule "), "plant_species"] <- "Cirsium acaule"
      ## clear spelling mistake
    data_traits[(data_traits$plant_species == "Dactylis glomerata "), "plant_species"] <- "Dactylis glomerata"
      ## clear spelling mistake
    data_traits[(data_traits$plant_species == "Euphrasia rostkoviana "), "plant_species"] <- "Euphrasia rostkoviana"
      ## clear spelling mistake
    data_traits[(data_traits$plant_species == "Festuca ovina agg."), "plant_species"] <- "Festuca ovina"
      ## clear spelling mistake
    data_traits[(data_traits$plant_species == "Geranium Sylvaticum"), "plant_species"] <- "Geranium sylvaticum"
      ## clear spelling mistake
    data_traits[(data_traits$plant_species == "Geranium Sylvaticum"), "plant_species"] <- "Geranium sylvaticum"
      ## clear spelling mistake
    data_traits[(data_traits$plant_species == "Helianthemum numularium"), "plant_species"] <- "Helianthemum nummularium"
      ## clear spelling mistake
    data_traits[(data_traits$plant_species == "Helictrotrichon pubescens"), "plant_species"] <- "Helictotrichon pubescens"
      ## clear spelling mistake
    data_traits[(data_traits$plant_species == "Leontodon hispidus "), "plant_species"] <- "Leontodon hispidus"
      ## clear spelling mistake
    data_traits[(data_traits$plant_species == "Leucanthemum Vulgare"), "plant_species"] <- "Leucanthemum vulgare"
      ## clear spelling mistake
    data_traits[(data_traits$plant_species == "Lotus corniculatus "), "plant_species"] <- "Lotus corniculatus"
      ## clear spelling mistake
    data_traits[(data_traits$plant_species == "Phleum alpinum "), "plant_species"] <- "Phleum alpinum"
      ## clear spelling mistake
    data_traits[(data_traits$plant_species == "Phleum alpinum agg."), "plant_species"] <- "Phleum alpinum"
      ## clear spelling mistake
    data_traits[(data_traits$plant_species == "Plantago atrata "), "plant_species"] <- "Plantago atrata"
      ## clear spelling mistake
    data_traits[(data_traits$plant_species == "Plantago lanceolota"), "plant_species"] <- "Plantago lanceolata"
      ## clear spelling mistake
    data_traits[(data_traits$plant_species == "Polygala chamaebuxud"), "plant_species"] <- "Polygala chamaebuxus"
      ## clear spelling mistake
    data_traits[(data_traits$plant_species == "Potentilla erecta "), "plant_species"] <- "Potentilla erecta"
      ## clear spelling mistake
    data_traits[(data_traits$plant_species == "Potentilla grandiflora "), "plant_species"] <- "Potentilla grandiflora"
      ## clear spelling mistake
    data_traits[(data_traits$plant_species == "Prunella grandiflora "), "plant_species"] <- "Prunella grandiflora"
    ## clear spelling mistake
    data_traits[(data_traits$plant_species == "Ranunculus acris "), "plant_species"] <- "Ranunculus acris"
      ## clear spelling mistake
    data_traits[(data_traits$plant_species == "Salvia pratensis "), "plant_species"] <- "Salvia pratensis"
      ## clear spelling mistake
    data_traits[(data_traits$plant_species == "Sesleria cearula"), "plant_species"] <- "Sesleria caerulea"
      ## clear spelling mistake
    data_traits[(data_traits$plant_species == "Sesleria cearulea "), "plant_species"] <- "Sesleria caerulea"
      ## clear spelling mistake
    data_traits[(data_traits$plant_species == "Sesleria caerula"), "plant_species"] <- "Sesleria caerulea"
      ## clear spelling mistake
    data_traits[(data_traits$plant_species == "Sesleria caerulea "), "plant_species"] <- "Sesleria caerulea"
      ## clear spelling mistake
    data_traits[(data_traits$plant_species == "Taraxacum officinale "), "plant_species"] <- "Taraxacum officinale"
      ## clear spelling mistake
    data_traits[(data_traits$plant_species == "Teucreum chamaedris"), "plant_species"] <- "Teucrium chamaedrys"
      ## clear spelling mistake
    data_traits[(data_traits$plant_species == "Teucreum montanum"), "plant_species"] <- "Teucrium montanum"
      ## clear spelling mistake
    data_traits[(data_traits$plant_species == "Teucrium chamaedris"), "plant_species"] <- "Teucrium chamaedrys"
      ## clear spelling mistake
    data_traits[(data_traits$plant_species == "Thymus puleginoides"), "plant_species"] <- "Thymus pulegioides"
      ## clear spelling mistake
    data_traits[(data_traits$plant_species == "Trifolium montanum "), "plant_species"] <- "Trifolium montanum"
      ## clear spelling mistake
    data_traits[(data_traits$plant_species == "Veronica chamaedris"), "plant_species"] <- "Veronica chamaedrys"
      ## clear spelling mistake
    data_traits[(data_traits$plant_species == "Viola reichenbachiana"), "plant_species"] <- "Viola riviana"
      ## species misidentification
    data_traits$plant_species <- factor(data_traits$plant_species)
      ## get rid of levels that are not used
    
    ## 2.3.2. Individual number
    table(data_traits$individual_nr)
      ## nothing special
    
    ## 2.3.3. Site
    levels(data_traits$site)
      ## displays categories in column plant species
    which(data_traits$site == "") 
    data_traits[(data_traits$site == "Arella "), "site"] <- "Arella"
      ## clear spelling mistake
    data_traits[(data_traits$site == "Bofel "), "site"] <- "Bofel"
      ## clear spelling mistake
    data_traits[(data_traits$site == "Calanda "), "site"] <- "Calanda"
      ## clear spelling mistake
    data_traits[(data_traits$site == "Nesselboden "), "site"] <- "Nesselboden"
      ## clear spelling mistake
    data_traits[data_traits$site == "Lower site", ]
      ## check coordinates of "lower site" - separate category at Bofel elevation (but not on Bofel)
    data_traits$site <- factor(data_traits$site)
      ## get rid of levels that are not used
    
    ## 2.3.4. Date
    table(data_traits$date)
      ## different date formats
    data_traits$date <- parse_date_time(data_traits$date, orders = c("d.m.y", "d.m.Y"))
      ## different date formats are parsed
    year(data_traits$date) <- 2024
      ## set all the years to 2024
    
    ## 2.3.5. Collector
    table(data_traits$collector)
    data_traits$collector <- gsub("CB1", "CB", data_traits$collector)
      ## change CB1 to CB
    
    ## 2.3.6. Coordinates
    colnames(data_traits)[7] <- "coordinates"
      ## change column name
    data_traits <- data_traits %>%
      separate(col = coordinates, into = c("x", "y", "z"), sep = " ", fill = "right", remove = TRUE)
      ## split coordinates into different columns
    data_traits$x <- gsub("[',]", "", data_traits$x)
    data_traits$y <- gsub("[',]", "", data_traits$y)
      ## remove ' and ,
    data_traits$coord_type <- NA
      ## make a new column for coordinate type
    data_traits$coord_type[is.na(data_traits$z)] <- "LV95"
    data_traits$coord_type[!is.na(data_traits$z)] <- "LV03"
      ## add coordinate type based on column z
    data_traits[data_traits$plant_id == 1989,"coord_type"] <- "LV03"
      ## change manual mistake in the coordinate format
    data_traits$x <- ifelse(
      data_traits$coord_type == "LV95",
      ifelse(
        grepl("\\.", as.character(data_traits$x)),
        as.character(data_traits$x),
        sub("^(\\d{7})(\\d+)$", "\\1.\\2", as.character(data_traits$x))
      ),
      data_traits$x
    )
      ## some x coordinate values had a "," instead of a "." as decimal sign (and it therefore got deleted) - now it's added again
    data_traits$y <- ifelse(
      data_traits$coord_type == "LV95",
      ifelse(
        grepl("\\.", as.character(data_traits$y)),
        as.character(data_traits$y),
        sub("^(\\d{7})(\\d+)$", "\\1.\\2", as.character(data_traits$y))
      ),
      data_traits$y
    )
      ## some y coordinate values had a "," instead of a "." as decimal sign (and it therefore got deleted) - now it's added again
    data_traits <- data_traits[, -9]
      ## delete column z
    data_traits$x <- as.numeric(data_traits$x)
    data_traits$y <- as.numeric(data_traits$y)
      ## change class of columns to numeric - seems like precision is missing when converting, but just displayed like that 
    convert_row <- function(x, y, type) {
      if (type == "LV95") {
        crs_code <- 2056
        x_in <- x
        y_in <- y
      } else { # LV03
        crs_code <- 21781
        x_in <- x
        y_in <- y
      }
      sf::st_transform(
        sf::st_sfc(
          sf::st_point(c(x_in, y_in)),
          crs = crs_code
        ),
        4326
      )[[1]] 
    }
      ## support function  
    wgs84_coords <- t(mapply(
      convert_row,
      data_traits$x,
      data_traits$y,
      data_traits$coord_type
    ))
      ## transformation for all cells
    data_traits$lon <- wgs84_coords[,1]
    data_traits$lat <- wgs84_coords[,2]
      ## add to data frame
    data_traits <- data_traits %>% select(-c(x,y,coord_type))
      ## remove old coordinate colums
    data_traits <- data_traits[, c(1:6, 33:32, 7:31)]
      ## rearrange column
    
    ## 2.3.7. Field measurements
    data_traits <- data_traits[!(is.na(data_traits$vegetative_height) & is.na(data_traits$reproductive_height)), ]
      ## delete all rows that have missing values in vegetative height & reproductive height
    
    ## 2.3.8. Flowering unit
    levels(data_traits$flowering_unit)
    table(data_traits$flowering_unit)
    data_traits[(!is.na(data_traits$flowering_unit) & data_traits$flowering_unit == ""), "flowering_unit"] <- NA
      ## add NA
    data_traits[(!is.na(data_traits$flowering_unit) & data_traits$flowering_unit == "capitulum"), "flowering_unit"] <- "Capitulum"
    data_traits[(!is.na(data_traits$flowering_unit) & data_traits$flowering_unit == " capitulum"), "flowering_unit"] <- "Capitulum"
    data_traits[(!is.na(data_traits$flowering_unit) & data_traits$flowering_unit == "\"Körbchen\""), "flowering_unit"] <- "Capitulum"
    data_traits[(!is.na(data_traits$flowering_unit) & data_traits$flowering_unit == "capitulum "), "flowering_unit"] <- "Capitulum"
    data_traits[(!is.na(data_traits$flowering_unit) & data_traits$flowering_unit == "Capitulum "), "flowering_unit"] <- "Capitulum"
    data_traits[(!is.na(data_traits$flowering_unit) & data_traits$flowering_unit == "flower"), "flowering_unit"] <- "Flower"
    data_traits[(!is.na(data_traits$flowering_unit) & data_traits$flowering_unit == "flower "), "flowering_unit"] <- "Flower"
    data_traits[(!is.na(data_traits$flowering_unit) & data_traits$flowering_unit == "Flower "), "flowering_unit"] <- "Flower"
    data_traits[(!is.na(data_traits$flowering_unit) & data_traits$flowering_unit == "flower this year"), "flowering_unit"] <- "Flower"
    data_traits[(!is.na(data_traits$flowering_unit) & data_traits$flowering_unit == "flowers"), "flowering_unit"] <- "Flower"
    data_traits[(!is.na(data_traits$flowering_unit) & data_traits$flowering_unit == "Flowers"), "flowering_unit"] <- "Flower"
    data_traits[(!is.na(data_traits$flowering_unit) & data_traits$flowering_unit == "inflorescence"), "flowering_unit"] <- "Inflorescence"
    data_traits[(!is.na(data_traits$flowering_unit) & data_traits$flowering_unit == "inflorescence "), "flowering_unit"] <- "Inflorescence"
    data_traits[(!is.na(data_traits$flowering_unit) & data_traits$flowering_unit == "Inflorescence "), "flowering_unit"] <- "Inflorescence"
    data_traits[(!is.na(data_traits$flowering_unit) & data_traits$flowering_unit == "Inflorescence count"), "flowering_unit"] <- "Inflorescence"
    data_traits[(!is.na(data_traits$flowering_unit) & data_traits$flowering_unit == "inflorescence lenght"), "flowering_unit"] <- "Inflorescence length"
    data_traits[(!is.na(data_traits$flowering_unit) & data_traits$flowering_unit == "inflorescence length"), "flowering_unit"] <- "Inflorescence length"
    data_traits[(!is.na(data_traits$flowering_unit) & data_traits$flowering_unit == "inflorescence length "), "flowering_unit"] <- "Inflorescence length"
    data_traits[(!is.na(data_traits$flowering_unit) & data_traits$flowering_unit == "Inflorescence length "), "flowering_unit"] <- "Inflorescence length"
    data_traits[(!is.na(data_traits$flowering_unit) & data_traits$flowering_unit == "Inflorescene"), "flowering_unit"] <- "Inflorescence"
    data_traits[(!is.na(data_traits$flowering_unit) & data_traits$flowering_unit == "inflorescens"), "flowering_unit"] <- "Inflorescence"
    data_traits[(!is.na(data_traits$flowering_unit) & data_traits$flowering_unit == "Inflorescense"), "flowering_unit"] <- "Inflorescence"
      ## clean spelling mistakes
    data_traits$flowering_unit <- factor(data_traits$flowering_unit)
      ## drop not used levels
    match_species_flowering_unit <-  data_traits %>% select(plant_species, flowering_unit)
    match_species_flowering_unit <- unique(match_species_flowering_unit)
    match_species_flowering_unit <- match_species_flowering_unit[order(match_species_flowering_unit$plant_species),]
      ## create a list to match species ID and flowering unit
    match_species_flowering_unit <- match_species_flowering_unit[!is.na(match_species_flowering_unit$flowering_unit),]
      ## drop rows with NA (since they have not been filled out in the field)
    rownames(match_species_flowering_unit) <- NULL
      ## reset rownames
    match_species_flowering_unit <- match_species_flowering_unit[-c(19, 25, 37, 40, 42, 56, 58:59, 65, 101, 106, 111, 117, 120),]
      ## manually choosing correct flowering unit type - some conflicts remaining (e.g. more than one measurement type per species, need to be tested manually)
    match_species_flowering_unit_del <- match_species_flowering_unit[duplicated(match_species_flowering_unit$plant_species), "plant_species"]
      ## species that need to be checked manually
    rownames(match_species_flowering_unit) <- NULL
    match_species_flowering_unit <- match_species_flowering_unit[-c(59:60, 67:70, 75:76),]
      ## delete species with two different flowering unit types out of list (to avoid conflicts with merging) - still stored in match_species_flowering_unit_del
    match_species_flowering_unit[match_species_flowering_unit$plant_species == "Scabiosa lucida", "flowering_unit"] <- "Inflorescence"
    match_species_flowering_unit[match_species_flowering_unit$plant_species == "Carum carvi", "flowering_unit"] <- "Inflorescence"
    match_species_flowering_unit[match_species_flowering_unit$plant_species == "Buphthalmum salicifolium", "flowering_unit"] <- "Capitulum"  
    match_species_flowering_unit[match_species_flowering_unit$plant_species == "Carex montana", "flowering_unit"] <- "Inflorescence length"
      ## manually change wrong wrong flowering unit type
    data_traits <- data_traits %>%
      left_join(match_species_flowering_unit, by = "plant_species")
      ## merges data sets, now twice the same column (flowering_unit.x and flowering_unit.y)
    data_traits$flowering_unit.y[is.na(data_traits$flowering_unit.y)] <- 
      data_traits$flowering_unit.x[is.na(data_traits$flowering_unit.y)]
    data_traits <- data_traits[, c(1:14, 34, 16:33)]
    colnames(data_traits)[15] <- "flowering_unit"
      ## rearrange columns, delete flowering_unit.x and rename flowering_unit.y
    
    ## 2.3.9. Phenology state
    levels(data_traits$phenology_state)
    table(is.na(data_traits$phenology_state))
    data_traits[(!is.na(data_traits$phenology_state) & data_traits$phenology_state == ""), "phenology_state"] <- NA
      ## "" to NA
    data_traits[(!is.na(data_traits$phenology_state) & data_traits$phenology_state == "BF "), "phenology_state"] <- "BF"
    data_traits[(!is.na(data_traits$phenology_state) & data_traits$phenology_state == "BFD "), "phenology_state"] <- "BFD"
    data_traits[(!is.na(data_traits$phenology_state) & data_traits$phenology_state == "BFDSE "), "phenology_state"] <- "BFDSE"
    data_traits[(!is.na(data_traits$phenology_state) & data_traits$phenology_state == "DSE "), "phenology_state"] <- "DSE"
    data_traits[(!is.na(data_traits$phenology_state) & data_traits$phenology_state == "F "), "phenology_state"] <- "F"
    data_traits[(!is.na(data_traits$phenology_state) & data_traits$phenology_state == "FD "), "phenology_state"] <- "FD"
    data_traits[(!is.na(data_traits$phenology_state) & data_traits$phenology_state == "FDSE "), "phenology_state"] <- "FDSE"
    data_traits[(!is.na(data_traits$phenology_state) & data_traits$phenology_state == "Inflorescence length"), "phenology_state"] <- NA
    data_traits[(!is.na(data_traits$phenology_state) & data_traits$phenology_state == "SE "), "phenology_state"] <- "SE"
    data_traits[(!is.na(data_traits$phenology_state) & data_traits$phenology_state == "SK"), "phenology_state"] <- "S"
      ## clear mistakes
    data_traits$phenology_state <- factor(data_traits$phenology_state)
    ## get rid of levels that are not used
    
    ## 2.3.10. Further counts/outliers
    boxplot(data_traits$flowering_unit_count ~ data_traits$plant_species)
    outliers <- data_traits %>%
      group_by(plant_species) %>%
      mutate(across(where(is.numeric),
                    ~. < quantile(., 0.25, na.rm=TRUE) - 3 * IQR(., na.rm=TRUE) |
                      . > quantile(., 0.75, na.rm=TRUE) + 3 * IQR(., na.rm=TRUE),
                    .names = "outlier_{.col}"))
    ## checking for unrealistic outliers
    outliers_check <- outliers[outliers$outlier_vegetative_height,]
      ## all values seem plausible
    outliers_check <- outliers[outliers$outlier_reproductive_height,]
      ## all values seem plausible
    outliers_check <- outliers[outliers$outlier_diameter_long,]
      ## all values seem plausible
    outliers_check <- outliers[outliers$diameter_90,]
      ## all values seem plausible
    outliers_check <- outliers[outliers$outlier_vegetative_height_stretched,]
      ## all values seem plausible
    outliers_check <- outliers[outliers$outlier_reproductive_height_stretched,]
      ## all values seem plausible
    outliers_check <- outliers[outliers$outlier_flowering_unit_count,]
      data_traits[data_traits$plant_id == 1187, "flowering_unit_count"] <- NA
      data_traits[data_traits$plant_id == 1187, "comment_field"] <- "The original inflorescence count was too high for the inflorescence length, so this value was removed."
      data_traits[data_traits$plant_id == 1673, "flowering_unit_count"] <- NA
      data_traits[data_traits$plant_id == 1673, "grass_spikelet_count"] <- 901
      data_traits[data_traits$plant_id == 1673, "comment_field"] <- "The original inflorescence count was too high for the inflorescence length, so this value was moved to the ‘grass_spikelet_count’ column."
      dauc_corr_inflorescence_count <- data_traits %>% filter(plant_id >= 1066 & plant_id <= 1070) %>% select(flowering_unit_count)
      dauc_corr_inflorescence_count <- as.vector(dauc_corr_inflorescence_count)
      dauc_corr_flowering_unit_count <- data_traits %>% filter(plant_id >= 1066 & plant_id <= 1070) %>% select(inflorescence_count)
      dauc_corr_flowering_unit_count <- as.vector(dauc_corr_flowering_unit_count)
      data_traits[data_traits$plant_id >= 1066 & data_traits$plant_id <= 1070, "flowering_unit_count"] <- dauc_corr_flowering_unit_count
      data_traits[data_traits$plant_id >= 1066 & data_traits$plant_id <= 1070, "inflorescence_count"] <- dauc_corr_inflorescence_count
        ## swapping confoundet columns
    outliers_check <- outliers[outliers$outlier_inflorescence_count,]
    outliers_check <- outliers_check[rowSums(is.na(outliers_check)) != ncol(outliers_check), ]
      ## all values seem plausible
    outliers_check <- outliers[outliers$outlier_grass_spikelet_count,]
    outliers_check <- outliers_check[rowSums(is.na(outliers_check)) != ncol(outliers_check), ]
      ## all values seem plausible
    outliers_check <- outliers[outliers$outlier_leaves_sampled,]
    outliers_check <- outliers_check[rowSums(is.na(outliers_check)) != ncol(outliers_check), ]
      ## all values seem plausible
    outliers_check <- outliers[outliers$outlier_biomass_fresh_weight,]
    outliers_check <- outliers_check[rowSums(is.na(outliers_check)) != ncol(outliers_check), ]
      ## all values seem plausible
    outliers_check <- outliers[outliers$outlier_flowering_unit_fresh_weight,]
    outliers_check <- outliers_check[rowSums(is.na(outliers_check)) != ncol(outliers_check), ]
      ## all values seem plausible
    outliers_check <- outliers[outliers$outlier_leaves_fresh_weight,]
    outliers_check <- outliers_check[rowSums(is.na(outliers_check)) != ncol(outliers_check), ]
      ## all values seem plausible
    outliers_check <- outliers[outliers$outlier_biomass_dry_weight,]
    outliers_check <- outliers_check[rowSums(is.na(outliers_check)) != ncol(outliers_check), ]
      ## all values seem plausible
    outliers_check <- outliers[outliers$outlier_flowering_unit_dry_weight,]
    outliers_check <- outliers_check[rowSums(is.na(outliers_check)) != ncol(outliers_check), ]
      ## all values seem plausible
    outliers_check <- outliers[outliers$outlier_leaves_dry_weight,]
    outliers_check <- outliers_check[rowSums(is.na(outliers_check)) != ncol(outliers_check), ]
      ## all values seem plausible
    data_traits_fresh_dry <- data_traits
      ## create new data file to check wether dry weight is always smaller than fresh weight
    data_traits_fresh_dry$biomass_weight_diff <- data_traits_fresh_dry$biomass_fresh_weight - data_traits_fresh_dry$biomass_dry_weight
    data_traits_fresh_dry$flowering_unit_weight_diff <- data_traits_fresh_dry$flowering_unit_fresh_weight - data_traits_fresh_dry$flowering_unit_dry_weight
    data_traits_fresh_dry$leaves_weight_diff <- data_traits_fresh_dry$leaves_fresh_weight - data_traits_fresh_dry$leaves_dry_weight
    summary(data_traits_fresh_dry$biomass_weight_diff)
    summary(data_traits_fresh_dry$flowering_unit_weight_diff)
    summary(data_traits_fresh_dry$leaves_weight_diff)
      ## there are some where the difference between fresh and dry weigh is negative and an error occured - needs to have a warning in the comments
    biomass_diff_index <- na.omit(data_traits_fresh_dry[data_traits_fresh_dry$biomass_weight_diff < 0, "plant_id"])
    flowering_unit_diff_index <- na.omit(data_traits_fresh_dry[data_traits_fresh_dry$flowering_unit_weight_diff < 0, "plant_id"])
    leaves_diff_index <- na.omit(data_traits_fresh_dry[data_traits_fresh_dry$leaves_weight_diff < 0, "plant_id"])
    diff_index <- unique(c(biomass_diff_index, flowering_unit_diff_index, leaves_diff_index))
      ## getting indices of the ones with negative differences to make a comment
    data_traits[data_traits$plant_id %in% diff_index, "comment_lab"] <- "The difference between the fresh weight and the dry weight is negative for at least one sample category. The cause of the error cannot be determined; possible causes include incorrect measurement, mixing up of samples or a malfunction of the scales."
      ## add a comment to all the ones with negative differences in fresh/dry weigh
    
    ## 2.3.11. Leaves from another plant
    levels(data_traits$leaves_other)
    table(data_traits$leaves_other)
    data_traits[(!is.na(data_traits$leaves_other) & data_traits$leaves_other == ""), "leaves_other"] <- "no"
    data_traits[(!is.na(data_traits$leaves_other) & data_traits$leaves_other == "No"), "leaves_other"] <- "no"
    data_traits[(!is.na(data_traits$leaves_other) & data_traits$leaves_other == "Yes"), "leaves_other"] <- "yes"
    data_traits[(is.na(data_traits$leaves_other)), "leaves_other"] <- "no"
    data_traits$leaves_other <- factor(data_traits$leaves_other)
    
    ## 2.3.12. Biomass collected
    colnames(data_traits)[22] <- "biomass_collected"
    ## rename column
    data_traits <- data_traits %>% select(!biomass_collected)
    ## after cleaning its visible, that one "no" is wrong and the others are at least partially filled out - delete this column, not needed anymore
    
    ## 2.3.13. Rosette count
    colnames(data_traits)[22] <- "rosette_count"
    data_traits$rosette_count <- as.numeric(data_traits$rosette_count)
    summary(data_traits$rosette_count)
    
    ## 2.3.15. CN
    data_traits <- data_traits %>% select(!cn)
      ## delete this column, not needed
    
    ## 2.3.16. Leaf scans
    levels(data_traits$leaves_scans)
    table(data_traits$leaves_scans)
    data_traits[(!is.na(data_traits$leaves_scans) & data_traits$leaves_scans == " yes"), "leaves_scans"] <- "yes"
    data_traits[(!is.na(data_traits$leaves_scans) & data_traits$leaves_scans == "Yes"), "leaves_scans"] <- "yes"
    data_traits[(!is.na(data_traits$leaves_scans) & data_traits$leaves_scans == "ÿes"), "leaves_scans"] <- "yes"
    data_traits[(!is.na(data_traits$leaves_scans) & data_traits$leaves_scans == "ys"), "leaves_scans"] <- "yes"
    data_traits[(!is.na(data_traits$leaves_scans) & data_traits$leaves_scans == "ywa"), "leaves_scans"] <- "yes"
    data_traits[(!is.na(data_traits$leaves_scans) & data_traits$leaves_scans == ""), "leaves_scans"] <- NA
    data_traits$leaves_scans <- factor(data_traits$leaves_scans)
    
    ## 2.3.17. Adressing weighing mistake in lab
      ## On several days, there was a mistake in processing biomass/flowering fresh weigh (in the IDs 1606-1675 and 1886-2025)
      ## Flowering unit with those IDs were cut before weighing biomass by mistake
      ## So the flowering unit weight needs to be added to biomass weight for those IDs
    data_traits$sum_test <- rowSums(data_traits[,c("biomass_fresh_weight", "flowering_unit_fresh_weight")], na.rm=TRUE) 
      ## create column with sum of biomass fresh weight and flowering unit fresh weight for all columns
    data_traits[data_traits$plant_id %in% c(1606:1675,1886:2025), "biomass_fresh_weight"] <- data_traits[data_traits$plant_id %in% c(1606:1675,1886:2025), "sum_test"] 
      ## overwrite biomass value with sum of biomass fresh weight and flowering unit fresh weight just in rows with known IDs
    data_traits$sum_test <- NULL
      ## remove temporary column
    
    ## 2.3.18. Add leaf weight to total biomass (fresh and dry) 
      ## just if "no" in column "leaves_other" (since a "yes" means)
    data_traits$test_sum_fresh <- NA
    data_traits$test_sum_dry <- NA
      ## create temporary columns for sums of total biomass and leaf biomass (one fresh and one dry)
    data_traits$test_sum_fresh <- rowSums(data_traits[,c("biomass_fresh_weight", "leaves_fresh_weight")], na.rm=TRUE) 
    data_traits$test_sum_dry <- rowSums(data_traits[,c("biomass_dry_weight", "leaves_dry_weight")], na.rm=TRUE) 
    data_traits[data_traits$leaves_other == "no", "biomass_fresh_weight"] <- data_traits[data_traits$leaves_other == "no", "test_sum_fresh"]
    data_traits[data_traits$leaves_other == "no", "biomass_dry_weight"] <- data_traits[data_traits$leaves_other == "no", "test_sum_dry"]
      ## transfer the sums into the column with the total biomass (just for the rows with "no" in "leaves_other")
    data_traits$test_sum_fresh <- NULL
    data_traits$test_sum_dry <- NULL
      ## delete temporary columns
    
## 3. Merging ====
    
    ## 3.1. Merge files
    data_traits <- left_join(data_traits, data_CN, by = "plant_id")
    
    ## 3.2. Adjust column names
    colnames(data_traits)[5] <- "date_field"
    colnames(data_traits)[32] <- "date_CN"
    colnames(data_traits)[33] <- "weight_CN_analysis"
    colnames(data_traits)[40] <- "comment_CN"

## 4. Exporting ====
    
    write.csv(data_traits, file = "/2024_TRAITS_CleanData/2025_TRAITS_clean_20251023.csv", row.names = FALSE)
                    
    
    