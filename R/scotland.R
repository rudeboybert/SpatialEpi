#' Lip Cancer in Scotland
#' 
#' County-level (n=56) data for lip cancer among males in Scotland between 1975-1980
#' 
#' 
#' @format 
#' List containing:
#' \describe{
#'   \item{geo}{a table of county IDs, x-coordinates (eastings) and y-coordinates (northings) of the geographic centroid of each county.}
#'   \item{data}{a table of county IDs, number of cases, population and strata information}
#'   \item{spatial.polygon}{a Spatial Polygons class (See \link[sp]{SpatialPolygons-class}) map of Scotland}
#'   \item{polygon}{a polygon map of Scotland (See \code{\link{polygon2spatial_polygon}}}
#' }
#' @source Kemp I., Boyle P., Smans M. and Muir C. (1985) Atlas of cancer in Scotland, 1975-1980, incidence and epidemiologic perspective \emph{International Agency for Research on Cancer} \bold{72}.
#'
#' @references Clayton D. and Kaldor J. (1987) Empirical Bayes estimates of age-standardized relative risks for use in disease mapping.  \emph{Biometrics}, \bold{43}, 671--681.
#' 
#' @examples
#' data(scotland)
#' data <- scotland$data
#' scotland.map <- scotland$spatial.polygon
#' SMR <- data$cases/data$expected
#' mapvariable(SMR,scotland.map)
#' 
#' 
#' 
#' 
#' 
"scotland"