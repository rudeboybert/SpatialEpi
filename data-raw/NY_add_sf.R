library(sf)
library(dplyr)
library(Hmisc)
library(usethis)
library(SpatialEpi)


NYleukemia_sf <- st_as_sf(NYleukemia$spatial.polygon)
ny_df <- NYleukemia$data
NYleukemia_sf <- merge(NYleukemia_sf, ny_df, by="row.names", all.x=TRUE) %>%
  select(population,cases,censustract.FIPS,geometry)

# Set projection: 
# https://stackoverflow.com/questions/61286108/error-in-cpl-transformx-crs-aoi-pipeline-reverse-ogrcreatecoordinatetrans
st_crs(NYleukemia_sf) <- 4326

usethis::use_data(NYleukemia_sf,overwrite = TRUE)
