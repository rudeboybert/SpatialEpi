# SpatialEpi 1.2.8

* Removed dependencies on `rgdal`, `rgeos`, and `maptools` https://r-spatial.org/r/2022/04/12/evolution.html, including
    + `install_geo_libraries.Rmd` vignette
    + `maptools::leglabs()`: we copied this function into our package, with original author Roger Bivand's permission and giving him full attribution


# SpatialEpi 1.2.7

* Removed outdated vignette with instructions on installing geos and gdal geospatial libraries
* Removed non-vignettes from vignettes folder


# SpatialEpi 1.2.6

* Fixed GDAL/PROJ requiring UTM zone error https://github.com/rudeboybert/SpatialEpi/issues/34


# SpatialEpi 1.2.5

* Fixed errors with 1.2.4 CRAN submission


# SpatialEpi 1.2.4

* Added `sf` package `scotland_sf`, `pennLC_sf`, `NYleukemia_sf` versions of 
`scotland`, `pennLC`, `NYleukemia` maps + datasets
* Converted all `man` documentation to be built using `roxygen2`
* Started updating main package vignette


# SpatialEpi 1.2.3

* Fixed SpatialEpi-package.Rd file formatting issues