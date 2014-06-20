SpatialEpi: Methods and Data for Spatial Epidemiology
===========================================
We load the data and convert the coordinate system from latitude/longitude to a 
grid-based system.  

```r
require(SpatialEpi)
data(NYleukemia)
map <- NYleukemia$spatial.polygon
centroids <- latlong2grid(NYleukemia$geo[, 2:3])
population <- NYleukemia$data$population
cases <- NYleukemia$data$cases
```



There were 4 census tracts that were completely surrounded by another census tract (an enclave like Lesotho).  We merge them into their surrounding census tracts and modify the population/case counts and SpatialPolygons object accordingly. 

```r
remove <- NYleukemia$surrounded
add <- NYleukemia$surrounding
population[add] <- population[add] + population[remove]
population <- population[-remove]
cases[add] <- cases[add] + cases[remove]
cases <- cases[-remove]
centroids <- centroids[-remove, ]
map <- SpatialPolygons(map@polygons[-remove], proj4string=CRS("+proj=longlat"))
```

We plot the incidence of leukemia for each census tract.  

```r
plotmap(cases/population, map, log=TRUE, nclr=5)
points(grid2latlong(centroids), pch=4)
```

<img src="figuresunnamed-chunk-3.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" width="768" />

