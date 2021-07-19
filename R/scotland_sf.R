#' Lip Cancer in Scotland
#' 
#' County-level (n=56) data for lip cancer among males in Scotland between 1975-1980
#' 
#' 
#' @format A data frame with 56 rows representing counties and 5 variables:
#' \describe{
#'   \item{geometry}{Geometric representation of counties in Scotland}
#'   \item{cases}{Number of Lip Cancer cases per county}
#'   \item{county.names}{Scotland County name}
#'   \item{AFF}{Proportion of the population who work in agricultural fishing and farming}
#'   \item{expected}{Expected number of lip cancer cases}
#' }
#' @source Kemp I., Boyle P., Smans M. and Muir C. (1985) Atlas of cancer in Scotland, 1975-1980, incidence and epidemiologic perspective \emph{International Agency for Research on Cancer} \bold{72}.
#'
#' @references Clayton D. and Kaldor J. (1987) Empirical Bayes estimates of age-standardized relative risks for use in disease mapping.  \emph{Biometrics}, \bold{43}, 671--681.
#' 
#' @examples
#' library(ggplot2)
#' ggplot() +
#' geom_sf(data = scotland_sf, aes(fill= cases))
#' 
#' 
#' 
"scotland_sf"