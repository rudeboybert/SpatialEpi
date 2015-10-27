# Save NYleukemia data as SpatialPolygonsDataFrame
library(rgeos)
load("/Users/rudeboybert/Dropbox/SpatialEpi/data/NYleukemia.rda")
sp <- NYleukemia$spatial.polygon
sp <- gSimplify(sp, tol = 0.00001)
sp <- gBuffer(sp, byid=TRUE, width=0)
# this is a well known R / GEOS hack (usually combined with the above) to 
# deal with "bad" polygons
# http://gis.stackexchange.com/questions/163445/r-solution-for-topologyexception-input-geom-1-is-invalid-self-intersection-er
data <- NYleukemia$data %>% 
  rename(FIPS = censustract.FIPS) %>% 
  mutate(ID = 1:nrow(NYleukemia$data)) %>% 
  select(ID, FIPS, population, cases)
NYleukemia <- SpatialPolygonsDataFrame(sp, data)
# devtools::use_data(NYleukemia, NYleukemia, overwrite=TRUE)


# Two different ways of loading in shapefiles
options(stringsAsFactors=FALSE)
library(maptools)
shapefile_name <- system.file("extdata/sids/sids.shp", package="SpatialEpi")
nc_sids <- maptools::readShapePoly(fn = shapefile_name, ID = "FIPSNO",
  proj4string = CRS("+proj=longlat +ellps=clrk66")
)

library(rgdal)
shapefile_dir <- "/Library/Frameworks/R.framework/Versions/3.2/Resources/library/spdep/etc/shapes/"
nc_sids <- rgdal::readOGR(dsn=shapefile_dir, layer="sids", p4s = "+proj=longlat +ellps=clrk66",
  verbose = FALSE
)



# Testing plot kulldorff
data("NYleukemia")
centroids <- sp::coordinates(NYleukemia)
kulldorff_output <- kulldorff(centroids = centroids,
                    cases = NYleukemia$cases, population = NYleukemia$population,
                    expected_cases = NULL, pop_upper_bound = 0.15, alpha_level=0.05,
                    n_sim=999)
sp_obj <- NYleukemia


