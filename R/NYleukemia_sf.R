#' Upstate New York Leukemia 
#' 
#' 
#' Census tract level (`n=281`) leukemia data for the 8 counties in upstate New York from 1978-1982, paired with population data from the 1980 census.  
#' Note that 4 census tracts were completely surrounded by another unique census tract; 
#' when applying the Bayesian cluster detection model in [bayes_cluster()],
#' we merge them with the surrounding census tracts yielding `n=277` areas.
#' 
#' 
#' @format An sf 'POLYGON' data frame with 281 rows  and 4 variables:
#' \describe{
#'   \item{geometry}{Geometric representation of 8 counties in upstate New York }
#'   \item{cases}{Number of cases per county}
#'   \item{population}{Population of each census tract}
#'   \item{censustract.FIPS}{11-digit Federal Information Processing System identification number for each county}
#'   
#' }
#' @source Turnbull, B. W. et al (1990) Monitoring for clusters of disease: application to leukemia incidence in upstate New York *American Journal of Epidemiology*, **132**, 136--143
#' 
#' @examples 
#' 
#' # Static map of NY Leukemia rate per county
#' library(ggplot2)
#' \dontrun{
#' ggplot(NYleukemia_sf) + 
#'   geom_sf(aes(fill= cases/population)) + 
#'   scale_fill_gradient(low = "white", high = "red")
#'   }
"NYleukemia_sf"
