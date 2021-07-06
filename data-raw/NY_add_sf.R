library(sf)
library(dplyr)
library(Hmisc)
library(usethis)
library(SpatialEpi)


NYleukemia_sf <- st_as_sf(NYleukemia$spatial.polygon)
ny_df <- NYleukemia$data
NYleukemia_sf <- merge(NYleukemia_sf, ny_df, by="row.names", all.x=TRUE) %>%
  select(population,cases,censustract.FIPS,geometry)

usethis::use_data(NYleukemia_sf,overwrite = TRUE)
