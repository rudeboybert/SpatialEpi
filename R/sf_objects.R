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
#' @source Population data was obtained from the 2000 decennial census, lung cancer and smoking data were obtained from the Pennsylvania Department of Health website:\url{http://www.dsf.health.state.pa.us/}.
"pennLC_sf"