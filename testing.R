# Save NYleukemia data as SpatialPolygonsDataFrame
load("/Users/rudeboybert/Dropbox/SpatialEpi/data/NYleukemia.rda")
sp <- NYleukemia$spatial.polygon
data <- NYleukemia$data %>% 
  rename(FIPS = censustract.FIPS) %>% 
  mutate(ID = 1:nrow(NYleukemia$data)) %>% 
  select(ID, FIPS, population, cases)

NYleukemia <- SpatialPolygonsDataFrame(sp, data)
devtools::use_data(NYleukemia, NYleukemia, overwrite=TRUE)


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