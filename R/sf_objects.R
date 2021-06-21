#' 
#' @title Pennsylvania Lung Cancer SF Object
#' @name pennLC_sf
#' @description 
#' pennLC_sf is a spatial object that can be used to create static or interactive maps.
#' 
#' @examples
#' Static Map
#' 
#' ggplot() + 
#' geom_sf(data = pennLC_sf)
#' 
#' Interactive Map using Leaflet
#' 
#' leaflet(pennLC_sf) %>%
#' addPolygons()
#' 
#' @return \code{pennLC_sf}
NULL