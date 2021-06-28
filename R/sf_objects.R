#' Pennsylvania Lung Cancer SF Object
#'
#' pennLC_sf is a spatial object that is used to create static or interactive maps.
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
"pennLC_sf"