#' SpatialEpi: Data and methods for spatial epidemiology.
#'
#' This package provides methods and Data for Spatial Epidemiology.
#'
#' @docType package
#' @name SpatialEpi
#' @useDynLib SpatialEpi
#' @importFrom Rcpp sourceCpp
#' @examples
#' # Example usage
#' library(SpatialEpi)
#' data(scotland)
#' map <- scotland$spatial.polygon
#' y <- scotland$data$cases
#' E <- scotland$data$expected
#' SMR <- y/E
#' # Plot SMR
#' plotmap(SMR, map, nclr=9, location="topleft") 
NULL
#> NULL