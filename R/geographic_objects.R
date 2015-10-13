# When using dplyr, R CMD CHECK returns NOTE that there is "no visible binding 
# for global variable", hence this hack fix.
# no visible binding for global variable
utils::globalVariables(c("distance", "prop_pop", "cum_prop", "."))

#' Define a study region's single zones.
#' 
#' Using the centroids of the n areas in a study region, define the set of 
#' single zones based on the upper bound \code{pop_upper_bound} of the 
#' proportion of the study region's population each single zone can contain.
#' 
#' @param centroids A data frame with 2 columns of the centroids of the n areas 
#'   in the study region.
#' @param population A vector of length n of the corresponding population 
#'   counts.
#' @param pop_upper_bound The upper bound of the proportion of the study 
#'   region's population the single zones can contain.
#'   
#' @return A list with two objects \describe{ 
#'   \item{\code{nearest_neighbors}}{A list of length n of the single zones 
#'   for each area} 
#'   \item{\code{cluster_coords}}{A data frame with 2 columns of the centering 
#'   and radial area for each single zone.} 
#'   }
#' @references Kulldorff, M. (1997) A spatial scan statistic. 
#'   \emph{Communications in Statistics: Theory and Methods}, \bold{26}, 
#'   1481--1496.
#' @export
#'   
#' @examples 
#' data(NYleukemia)
#' centroids <- sp::coordinates(NYleukemia)
#' single_zones <- define_single_zones(centroids, NYleukemia$population, 0.15)
define_single_zones <- function(centroids, population, pop_upper_bound) {
  
  # Number of areas
  n <- nrow(centroids)
  # Interpoint distance matrix
  dist_matrix <- 
    centroids %>% 
    dist(upper = TRUE, diag = TRUE) %>% 
    as.matrix()
  
  # Area indices and proportion of the population each area contains
  areas <- data.frame(
    area = 1:n,
    prop_pop= population/sum(population)
  )
  
  # For each area, we compute the list of neighbors in order of distance such
  # that the proportion of the study population that is included is less than or
  # equal to pop_upper_bound. Save the results in a list
  nearest_neighbors <- vector(mode = "list", length = n)
  
  for (i in 1:n) {
    nearest_neighbors[[i]] <- 
      areas %>% 
      mutate(distance = dist_matrix[, i]) %>% 
      arrange(distance) %>% 
      mutate(cum_prop = cumsum(prop_pop)) %>% 
      filter(cum_prop <= pop_upper_bound) %>% 
      .[["area"]]
  }
  
  # Save results in n_zones by 2 data frame denote the center and radial area
  # for each single zone
  cluster_coords <- data.frame(
    center = rep(1:n, times = sapply(nearest_neighbors, length)),
    radius = unlist(nearest_neighbors)
  )
  
  # Output results
  results <- list(
    nearest_neighbors = nearest_neighbors, 
    cluster_coords = cluster_coords
    )
  return(results)
}


#' Obtain all areas included in a single zone
#'
#' @param cluster_index integer index of single zone
#' @param zone_info output of \code{\link{define_single_zones}} function
#'
#' @return vector of areas
#' @export
#'
#' @examples
#' 1+1
return_single_zone_areas <- function(cluster_index, zone_info){
  # Obtain single zone center and radial areas
  center <- zone_info$cluster_coords[cluster_index, 1]
  radial <- zone_info$cluster_coords[cluster_index, 2]
  
  # Obtain all areas (in order of distance) from center to radial area
  cluster <- zone_info$nearest_neighbors[[center]]
  cluster <- cluster[1:which(cluster == radial)]
  
  return(cluster)
}



#' Create geo objects needed for Bayesian method
#'
#' @inheritParams define_single_zones
#'
#' @return A list with three objects \describe{ 
#'   \item{\code{presence}}{List of length n, for each area, the single zones it
#'   is present in.}
#'   \item{\code{cluster_list}}{List of length n_zones indicating the areas
#'   included in each single zone}
#'   \item{\code{cluster_coords}}{A data frame with 2 columns of the centering 
#'   and radial area for each single zone.}
#'   }
#' @export
#'
#' @examples
#' data(NYleukemia)
#' centroids <- sp::coordinates(NYleukemia)
#' geo_objects <- create_geo_objects(centroids, NYleukemia$population, 0.15)
create_geo_objects <- function(centroids, population, pop_upper_bound){
  
  # Number of areas
  n <- nrow(centroids)
  
  # Catch error
  #     if(max.prop < max(normalize(population))){
  #       print(paste("max.prop needs to be at least", max(normalize(population))))
  #     }
  
  # Define single zones
  zone_info <- define_single_zones(centroids, population, pop_upper_bound)
  nearest_neighbors <- zone_info$nearest_neighbors
  cluster_coords <- zone_info$cluster_coords
  n_zones <- nrow(cluster_coords)
  
  # 1. Create list of length n_zones indicating the component areas for each
  # zone
  cluster_list <- vector(mode="list", length=n_zones)
  counter <- 1
  for(i in 1:n) {
    nn <- nearest_neighbors[[i]]
    for(j in 1:length(nn)) {
      cluster_list[[counter]] <- nn[1:j] 
      counter <- counter + 1  
    } 
  }
  
  # 2. Generate overlap object which tracks the overlap between single zones
  # For each area, list all single zones that it is included in
  presence <- vector(mode="list", length=n)
  for(i in 1:n){
    presence[[i]] <- cluster_list %>% 
      sapply(function(x){is.element(i, x)}) %>% 
      which()
  }
  
  # Output results
  results <- list(
    presence = presence, 
    cluster_list = cluster_list,
    cluster_coords = cluster_coords
  )
  return(results)
}


