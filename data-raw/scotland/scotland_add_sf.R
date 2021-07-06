library(sf)
library(dplyr)
library(Hmisc)
library(usethis)
library(SpatialEpi)

scotland_sf <- st_as_sf(scotland$spatial.polygon)
scotland_df <- scotland$data
scotland_sf <- merge(scotland_sf, scotland_df, by="row.names", all.x=TRUE) %>%
  select(county.names,cases,expected,AFF,geometry)

usethis::use_data(scotland_sf,overwrite = TRUE)
