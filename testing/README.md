
SpatialEpi
===========================================
Methods and Data for Spatial Epidemiology
-------------------------------------------
We load the data and convert the coordinate system from latitude/longitude to a 
grid-based system.  

```r
library(SpatialEpi)
```

```
## Loading required package: sp
```

```r
data(NYleukemia)
sp.obj <- NYleukemia$spatial.polygon
centroids <- latlong2grid(NYleukemia$geo[, 2:3])
population <- NYleukemia$data$population
cases <- NYleukemia$data$cases
```

We plot the incidence of leukemia for each census tract.  

```r
plotmap(cases/population, sp.obj, log=TRUE, nclr=5)
points(grid2latlong(centroids), pch=4)
```

![plot of chunk unnamed-chunk-2](./figure/unnamed-chunk-2.png) 

We run the Bayesian Cluster Detection method from Wakefield and Kim (2013)

```r
y <- cases
E <- expected(population, cases, 1)
max.prop <- 0.15
shape <- c(2976.3, 2.31)
rate <- c(2977.3, 1.31)
J <- 7
pi0 <- 0.95
n.sim.lambda <- 10^4
n.sim.prior <- 10^5
n.sim.post <- 10^5

# Compute output
output <- bayes_cluster(y, E, population, sp.obj, centroids, max.prop, 
 shape, rate, J, pi0, n.sim.lambda, n.sim.prior, n.sim.post)
```

```
## [1] "Algorithm started on: Wed Aug 13 07:49:25 2014"
## [1] "Geographic objects creation complete on: Wed Aug 13 07:49:36 2014"
```

```
## Warning: Walker's alias method used: results are different from R < 2.2.0
```

```
## [1] "Importance sampling of lambda complete on: Wed Aug 13 07:50:07 2014"
## [1] "Prior map MCMC complete on: Wed Aug 13 07:54:30 2014"
## [1] "Posterior estimation complete on: Wed Aug 13 08:08:16 2014"
```

```r
plotmap(output$post.map$high.area, sp.obj)
```

![plot of chunk unnamed-chunk-3](./figure/unnamed-chunk-3.png) 


