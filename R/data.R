#' Leukemia in Upstate New York.
#'
#' Leukemia data for 281 census tracts in 7 counties in Upstate New York.
#'
#' @format object of class \code{\link[sp]{SpatialPolygonsDataFrame-class}} with a data frame of
#' 281 observations and 4 variables:
#' \itemize{
#'   \item ID: id variable
#'   \item FIPS: Federal Information Processing Standards code of census tract
#'   \item population: population size (1980 U.S. Census)
#'   \item cases: number of cases 1978-1982
#' }
#' @references Waller, L. and C. Gotway (2004) Applied Spatial Statistics for
#' Public Health Data. New York: John Wiley and Sons.
#' @source \url{http://web1.sph.emory.edu/users/lwaller/ch9index.htm}
"NYleukemia"


#' Pennsylvania Lung Cancer.
#'
#' County-level (n=67) population/case data for lung cancer in Pennsylvania in 
#' 2002, stratified on race (white vs non-white), gender and age (Under 40, 
#' 40-59, 60-69 and 70+).  Additionally, county-specific smoking rates. 
#'
#' @format object of class \code{\link[sp]{SpatialPolygonsDataFrame-class}} with a data frame of
#' 281 observations and 4 variables:
#' \itemize{
#'   \item ID: id variable
#'   \item FIPS: Federal Information Processing Standards code of census tract
#'   \item population: population size (1980 U.S. Census)
#'   \item cases: number of cases 1978-1982
#' }
#' @source Population data was obtained from the 2000 decennial census, lung 
#' cancer and smoking data were obtained from the Pennsylvania Department of 
#' Health website: \url{http://www.dsf.health.state.pa.us/}
"pennLC"


#' Lip Cancer in Scotland.
#'
#' County-level (n=56) data for lip cancer among males in Scotland between 
#' 1975-1980.
#'
#' @format object of class \code{\link[sp]{SpatialPolygonsDataFrame-class}} with a data frame of
#' 281 observations and 4 variables:
#' \itemize{
#'   \item ID: id variable
#'   \item FIPS: Federal Information Processing Standards code of census tract
#'   \item population: population size (1980 U.S. Census)
#'   \item cases: number of cases 1978-1982
#' }
#' @source Kemp I., Boyle P., Smans M. and Muir C. (1985) Atlas of cancer in 
#' Scotland, 1975-1980, incidence and epidemiologic perspective 
#' \emph{International Agency for Research on Cancer} \bold{72}.
#' @references Clayton D. and Kaldor J. (1987) Empirical Bayes estimates of 
#' age-standardized relative risks for use in disease mapping.
#' \emph{Biometrics}, \bold{43}, 671--681
"scotland"


