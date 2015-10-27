## ---- eval=FALSE---------------------------------------------------------
#  library(rgdal)
#  library(rgeos)

## ---- message=FALSE, collapse=TRUE, tidy=TRUE, tidy.opts=list(width.cutoff=40), cache=TRUE----
shapefile_name <- system.file("extdata/sids/sids.shp", package="SpatialEpi")
proj4string <- sp::CRS("+proj=longlat +ellps=clrk66")
nc_sids <- maptools::readShapePoly(fn = shapefile_name, ID = "FIPSNO",
    proj4string = proj4string)
map_variable(nc_sids)
names(nc_sids)
nc_sids$RATE74 <- nc_sids$SID74/nc_sids$BIR74

