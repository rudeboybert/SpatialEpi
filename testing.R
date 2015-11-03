# Save scotland data as SpatialPolygonsDataFrame
library(rgdal)
shapefile_dir <- "./inst/extdata/scotland/"
scotland <- rgdal::readOGR(dsn=shapefile_dir, layer="scot", p4s = "+proj=longlat +ellps=clrk66",
                           verbose = FALSE
)

data <- read.table("./inst/extdata/scotland/scotland.dat", skip=1)
colnames(data) <- c("District", "Observed", "Expected", "AFF", "Latitude", "Longitude")

scotland@data <- inner_join(scotland@data, data, by=c("ID"="District"))
scotland$SMR <- scotland$Observed/scotland$Expected
map_variable(scotland, "SMR", type="jenks")




library(stringr)
load("/Users/aykim/Dropbox/SpatialEpi/data/scotland.rda")
sp <- scotland$spatial.polygon

getSpPPolygonsIDSlots(sp)
names <- sapply(slot(sp, "polygons"), function(x) slot(x, "ID")) %>% 
  str_replace_all("\\.", "_")

for(i in 1:length(sp)){
  slot(slot(sp, "polygons")[[i]], "ID") <- as.character(i)
}

data <- scotland$data %>% 
  mutate(county_names=str_replace_all(county.names, "\\.", "_")) %>% 
  select(county_names, cases, expected, AFF)

scotland <- SpatialPolygonsDataFrame(sp, data)
scotland$SMR <- scotland$cases/scotland$expected
map_variable(scotland, "SMR", type="equidistant")
# devtools::use_data(scotland, scotland, overwrite=TRUE)



# Save NYleukemia data as SpatialPolygonsDataFrame
library(rgeos)
load("/Users/rudeboybert/Dropbox/SpatialEpi/data/NYleukemia.rda")
sp <- NYleukemia$spatial.polygon
sp <- gSimplify(sp, tol = 0.00001)
sp <- gBuffer(sp, byid=TRUE, width=0)
# this is a well known R / GEOS hack (usually combined with the above) to 
# deal with "bad" polygons
# http://gis.stackexchange.com/questions/163445/
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

library(rgdal)
shapefile_dir <- "/Users/aykim/Downloads/"
scotland <- rgdal::readOGR(dsn=shapefile_dir, layer="scot", p4s = "+proj=longlat +ellps=clrk66",
                          verbose = FALSE
)



