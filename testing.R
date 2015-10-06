# Save NYleukemia data as SpatialPolygonsDataFrame
load("/Users/rudeboybert/Dropbox/SpatialEpi/data/NYleukemia.rda")
sp <- NYleukemia$spatial.polygon
data <- NYleukemia$data %>% 
  rename(FIPS = censustract.FIPS) %>% 
  mutate(ID = 1:nrow(NYleukemia$data)) %>% 
  select(ID, FIPS, population, cases)

NYleukemia <- SpatialPolygonsDataFrame(sp, data)
devtools::use_data(NYleukemia, NYleukemia, overwrite=TRUE)



# Testing define_single_zones
data("NYleukemia")
centroids <- coordinates(NYleukemia)
population <- NYleukemia$population
pop_upper_bound <- 0.15
testing <- define_single_zones(centroids, population, pop_upper_bound)



# Testing create_geo_objects
testing <- create_geo_objects(centroids, population, pop_upper_bound)



# Testing kulldorff
data("NYleukemia")
centroids <- coordinates(NYleukemia)
cases <- NYleukemia$cases
population <- NYleukemia$population
expected_cases <- NULL
pop_upper_bound <- 0.15
alpha_level <- 0.05
n_sim <- 9999
testing <- kulldorff(centroids, cases, population, expected_cases, pop_upper_bound, 
                     alpha_level=0.05)