# SpatialEpi 1.2.2.9000

## R Style/Conventions

* Switched to building package with devtools as per http://r-pkgs.had.co.nz/
* All `.`'s in variable and function names are replaced with `_`'s.
* `dplyr` package and `magrittr` package `%>%` pipe opeators are used wherever
possible for code legibility. 

## Cluster Detection

* All cluster detection methods (`kulldorff`, `besag_newell`, and
`bayes_cluster`) no longer require `sp` spatial polygons objects; only a two-column
matrix of each area's centroid.
*



# SpatialEpi 1.2.2

* CRAN v 1.2.2 uploaded to GitHub