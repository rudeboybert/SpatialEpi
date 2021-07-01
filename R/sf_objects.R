#' Pennsylvania Lung Cancer
#'
#' County-level (n=67) population/case data for lung cancer in 
#' Pennsylvania in 2002, stratified on race (white vs non-white), 
#' gender and age (Under 40, 40-59, 60-69 and 70+). Additionally, 
#' county-specific smoking rates.
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