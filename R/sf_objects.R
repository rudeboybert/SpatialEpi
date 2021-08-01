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
#'  \item{smoking}{Overall county smoking rate (not broken down by strata)}
#'  \item{geometry}{Geometric representation of counties in Pennsylvania}
#' }
#' @source Population data was obtained from the 2000 decennial census, lung cancer and smoking data were obtained from the Pennsylvania Department of Health website:<https://www.health.pa.gov/Pages/default.aspx>.
#' @examples 
#' library(ggplot2)
#' library(dplyr)
#' # Sum cases & population for each county
#' lung_cancer_rate <- pennLC_sf %>% 
#'   group_by(county) %>% 
#'   summarize(cases = sum(cases), population = sum(population)) %>% 
#'   mutate(rate = cases/population)
#' 
#' # Static map of Pennsylvania lung cancer rates for each county
#' \dontrun{
#' ggplot() +
#'   geom_sf(data = lung_cancer_rate, aes(fill = rate))
#'   }
"pennLC_sf"
