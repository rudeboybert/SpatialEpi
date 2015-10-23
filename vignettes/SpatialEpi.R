## ---- eval=FALSE---------------------------------------------------------
#  library(rgdal)
#  library(rgeos)

## ---- message=FALSE, collapse=TRUE, tidy=TRUE, tidy.opts=list(width.cutoff=40), cache=TRUE----
library(spdep)
data(nc.sids)

library(rgdal)
shapefile_dir <- "/Library/Frameworks/R.framework/Versions/3.2/Resources/library/spdep/etc/shapes/"
nc_sids <- readOGR(
  dsn=shapefile_dir, 
  layer="sids", 
  p4s = "+proj=longlat +ellps=clrk66",
  verbose = FALSE
  )

# Read in shapefile
options(stringsAsFactors=FALSE)
library(maptools)
shapefile_name <- system.file("etc/shapes/sids.shp", package="spdep")
shapefile_name <- "/Library/Frameworks/R.framework/Versions/3.2/Resources/library/spdep/etc/shapes/sids.shp"
nc_sids <- readShapePoly(
  fn = shapefile_name,
  ID = "FIPSNO",
  proj4string = CRS("+proj=longlat +ellps=clrk66")
  )

