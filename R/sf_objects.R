#' Pennsylvania Lung Cancer
#'
#' County-level (n=67) population/case data for lung cancer in 
#' Pennsylvania in 2002, stratified on race (white vs non-white), 
#' gender and age (Under 40, 40-59, 60-69 and 70+). Additionally, 
#' county-specific smoking rates.
#' 
#' @format An sf `POLYGON` data frame with 1072 rows = 67 counties x 2 race 
#' x 2 gender x 4 age bands
#' \describe{
#'  \item{county}{Pennsylvania county}
#'  \item{cases}{Number of cases per county split by strata}
#'  \item{population}{Population per county split by strata}
#'  \item{race}{Race (w = white and o = non-white)}
#'  \item{gender}{Gender (f = female and m = male)}
#'  \item{age}{Age (4 bands)}
#'  \item{geometry}{Geometric representation of counties in Pennsylvania}
#' }
#' @source Population data was obtained from the 2000 decennial census, lung cancer and smoking data were obtained from the Pennsylvania Department of Health website:\url{http://www.dsf.health.state.pa.us/}.
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
"pennLC_sf"