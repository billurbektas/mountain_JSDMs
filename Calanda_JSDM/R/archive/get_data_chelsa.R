library("terra")
library("tidyverse")
library("readxl")
library("lubridate")
library("e1071")
library("DBI")
library("RSQLite")
library("visNetwork")
library("sf")
library("devtools")
#install_git("https://gitlabext.wsl.ch/karger/rchelsa.git")
library(Rchelsa)
library(reticulate)
py_install("chelsa-cmip6", pip=T)
chelsa_cmip6 <- import('chelsa_cmip6')

#https://gitlabext.wsl.ch/karger/chelsa_cmip6/-/tree/master/



get_df = function(df){
  df = df %>%
    pivot_longer(
      cols = -time)%>%
    mutate(
      lon = as.numeric(str_extract(name, "(?<=lon:)[0-9.-]+")),  # Extract longitude
      lat = as.numeric(str_extract(name, "(?<=lat:)[0-9.-]+"))  # Extract latitude
    )%>%
    dplyr::select(-name)%>%
    left_join(metadata %>% dplyr::select(gradient, destSiteID, lon, lat, elevation, YearEstablished))%>%
    rename(site_name = destSiteID)
  return(df)
}


# GET coordinates ----

coords = data.frame(lon = metadata$lon, lat = metadata$lat)

az_coords = metadata %>% filter(gradient == "US_Arizona")%>%
  select(gradient, destSiteID, lon, lat)

co_coords = metadata %>% filter(gradient == "US_Colorado")%>%
  select(gradient, destSiteID, lon, lat)

mon_coords = metadata %>% filter(gradient == "US_Montana")%>%
  select(gradient, destSiteID, lon, lat)

eu_coords = metadata %>% filter(Country %in% c("Germany", "France", "Italy", "Switzerland"))%>%
  select(gradient, destSiteID, lon, lat)

no_coords = metadata %>% filter(Country %in% c("Norway"))%>%
  select(gradient, destSiteID, lon, lat)

sw_coords = metadata %>% filter(Country %in% c("Sweden"))%>%
  select(gradient, destSiteID, lon, lat)

# Get Future - USA ----
get_future_clim(coords = co_coords, scenario = 'ssp370', 
                output = "co_chelsa", year = 2022:2100)

get_future_clim(coords = az_coords, scenario = 'ssp370', 
                output = "az_chelsa", year = 2021:2100)

get_future_clim(coords = mon_coords, scenario = 'ssp370', 
                output = "mon_chelsa", year = 2021:2100)

# Get Future - EU ----
get_future_clim(coords = eu_coords, scenario = 'ssp370',
                output = "eu_chelsa", year = 2021:2100)

# Get Future - Scandinavia ----
get_future_clim(coords = no_coords, scenario = 'ssp370',
                output = "no_chelsa", year = 2091:2100)

get_future_clim(coords = sw_coords, scenario = 'ssp370',
                output = "sw_chelsa", year = 2093:2100)

# MERGE Future files ----
chelsa.future = list("eu", "no", "az", "mon", "co", "sw")%>%
  map(get_future)

get_future("az")
get_future("mon")
get_future("co")
get_future("sw")

# GET Present ----
## Temperature ----
tas = getChelsa('tas',coords=coords, startdate=as.Date("2005-1-1"), enddate=as.Date("2021-1-1"))
tasc = getChelsa('tas',coords=coords, startdate=as.Date("2021-1-1"), enddate=as.Date("2021-1-31"))

tasx = 
  list(tas, tasc) %>% 
  map(get_df) %>% 
  bind_rows()
  
tasx = 
  left_join(get.experiments, tasx)%>%
  distinct()

check.tas =
  tasx %>%
  group_by(experiment, level, site_name, lon, lat )%>%
  summarize(count = n()) # seems OK

write.csv(tasx, file = "output/tas_temperature_daily_present.csv")
saveRDS(tasx, file = "output/tas_temperature_daily_present.rds")

## Precipitation ----
pr = getChelsa('pr',coords=coords, startdate=as.Date("2005-1-1"), enddate=as.Date("2020-1-1"))
prc = getChelsa('pr',coords=coords, startdate=as.Date("2020-1-1"), enddate=as.Date("2020-1-4"))
prc2 = getChelsa('pr',coords=coords, startdate=as.Date("2020-1-6"), enddate=as.Date("2020-3-1"))
prc3 = getChelsa('pr',coords=coords, startdate=as.Date("2020-3-3"), enddate=as.Date("2020-4-4"))
prc4 = getChelsa('pr',coords=coords, startdate=as.Date("2020-4-6"), enddate=as.Date("2020-5-1"))
prc5 = getChelsa('pr',coords=coords, startdate=as.Date("2020-5-3"), enddate=as.Date("2020-6-1"))
prc6 = getChelsa('pr',coords=coords, startdate=as.Date("2020-6-3"), enddate=as.Date("2020-11-4"))
prc7 = getChelsa('pr',coords=coords, startdate=as.Date("2020-11-6"), enddate=as.Date("2020-12-31"))

prx = 
  list(pr, prc, prc2, prc3, prc4, prc5, prc6, prc7) %>% 
  map(get.df) %>% 
  bind_rows()

prx = 
  left_join(get.experiments, prx)%>%
  distinct()

check.pr =
  prx %>%
  group_by(experiment, level, site_name, lon, lat )%>%
  summarize(count = n()) # seems OK

write.csv(prx, file = "output/pr_precipitation_daily_present.csv")
saveRDS(prx, file = "output/pr_precipitation_daily_present.rds")

