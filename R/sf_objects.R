#' Pennsylvania Lung Cancer SF Object
#'
#' pennLC_sf is a spatial object that is used to create static or interactive maps for lunch cancer in Pennsylvania in 2002.
#' 
#' @format 
#' \describe{
#'  \item{geometry}{Counties in Pennsylvania}
#'  \item{cases}{Number of cases}
#'  \item{county}{Pennsylvania county}
#'  \item{population}{Population of the county}
#'  \item{race}{Race of the person (w = white and o = non-white) }
#'  \item{gender}{Gender of the person (f = female and m = male)}
#'  \item{age}{Age of the person}
#' 
#' }
#' 
#' 
#' 
#' 
#' 
#' @examples
#' Static Map
#' 
#' library(ggplot2)
#' ggplot() + 
#' geom_sf(data = pennLC_sf)
#' 
#' Interactive Map using Leaflet
#' 
#' library(leaflet)
#' leaflet(pennLC_sf) %>%
#' addPolygons()
#' 
"pennLC_sf"