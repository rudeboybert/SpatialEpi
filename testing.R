# Testing define_zones
load("/Users/rudeboybert/Dropbox/SpatialEpi/data/NYleukemia.rda")
geo <- NYleukemia$geo[,2:3] %>% as.data.frame()
population <- NYleukemia$data$population
pop_upper_bound <- 0.15

testing <- define_zones(geo, population, pop_upper_bound)
testing2 <- zones(geo, population, pop_upper_bound)



# Save NYleukemia data as SpatialPolygonsDataFrame
load("/Users/rudeboybert/Dropbox/SpatialEpi/data/NYleukemia.rda")

sp <- NYleukemia$spatial.polygon
data <- NYleukemia$data %>% 
  rename(FIPS = censustract.FIPS) %>% 
  mutate(ID = 1:nrow(NYleukemia$data)) %>% 
  select(ID, FIPS, population, cases)

NYleukemia <- SpatialPolygonsDataFrame(sp, data)
devtools::use_data(NYleukemia, NYleukemia, overwrite=TRUE)

library(sp)
testing <- define_zones(coordinates(NYleukemia), NYleukemia$population, 0.15)
