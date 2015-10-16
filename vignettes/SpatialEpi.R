## ---- eval=FALSE---------------------------------------------------------
#  library(rgdal)
#  library(rgeos)

## ---- message=FALSE, collapse=TRUE, tidy=TRUE,tidy.opts=list(width.cutoff=40)----
library(maptools)
library(dplyr)
library(spdep)
# Awkward
data(nc.sids)
sids <- data.frame(
  Observed=nc.sids$SID74,
  Expected=nc.sids$BIR74*sum(nc.sids$SID74)/sum(nc.sids$BIR74),
  x=nc.sids$x, 
  y=nc.sids$y)

nc.sids <- readShapePoly(
  fn = system.file("etc/shapes/sids.shp", package="spdep")[1],
  ID="FIPSNO",
  proj4string=CRS("+proj=longlat +ellps=clrk66"))

# Create a copy to modify
nc.sids@data <- nc.sids@data %>% 
  mutate(
    Y = SID74,
    E = sum(Y) * BIR74/sum(BIR74),
    SMR74 = Y/E,
    EXP74 = E
  )
brks <- seq(0, 5)

## ---- echo=TRUE, collapse=TRUE,fig.height=4.5,fig.width=4, fig.cap="Map of SMRs for SIDS in 1974 in North Carolina", tidy.opts=list(width.cutoff=35)----
spplot(nc.sids, "SMR74", at=brks, col.regions=grey.colors(5,start=.9,end=.1))

## ---- fig.show='hold'----------------------------------------------------
plot(1:10)
plot(10:1)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(mtcars, 10))

