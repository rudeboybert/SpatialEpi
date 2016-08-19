# R CMD CHECK returns NOTE that there is "no visible binding 
# for global variable" for variables used by dplyr/ggplot2, hence this hack fix.
utils::globalVariables(c("distance", "prop_pop", "cum_prop", "."))



#' Define a study region's single zones.
#' 
#' Using the centroids of the n areas in a study region, define the set of 
#' single zones based on the upper bound \code{pop_upper_bound} of the 
#' proportion of the study region's population each single zone can contain.
#' 
#' @param centroids data frame with 2 columns of the centroids of the n areas 
#'   in the study region.
#' @param population vector of length n of the corresponding population 
#'   counts.
#' @param pop_upper_bound upper bound of the proportion of the study 
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
#' @seealso \code{\link{create_geo_objects}}
#' @export
#' @importFrom magrittr %>%
#'   
#' @examples 
#' data(NYleukemia)
#' centroids <- sp::coordinates(NYleukemia$spatial.polygon)
#' single_zones <- define_single_zones(centroids, NYleukemia$data$population, 0.15)
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
      dplyr::mutate(distance = dist_matrix[, i]) %>% 
      dplyr::arrange(distance) %>% 
      dplyr::mutate(cum_prop = cumsum(prop_pop)) %>% 
      dplyr::filter(cum_prop <= pop_upper_bound) %>% 
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
#' @param single_zone_index integer index of single zone
#' @param zone_info output of \code{\link{define_single_zones}} defining single
#'   zones
#'
#' @return vector of areas
#' @export
#'
#' @examples
#' 1+1
return_single_zone_areas <- function(single_zone_index, zone_info){
  # Obtain single zone center and radial areas
  center <- zone_info$cluster_coords[single_zone_index, 1]
  radial <- zone_info$cluster_coords[single_zone_index, 2]
  
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
#' @seealso \code{\link{define_single_zones}}
#' @importFrom magrittr %>%
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
  
  # Create list of length n_zones indicating the component areas for each zone
  cluster_list <- vector(mode="list", length=n_zones)
  counter <- 1
  for(i in 1:n) {
    nn <- nearest_neighbors[[i]]
    for(j in 1:length(nn)) {
      cluster_list[[counter]] <- nn[1:j]
      counter <- counter + 1
    }
  }
  
  # Generate overlap object which tracks the overlap between single zones. For
  # each area, list all single zones that it is included in
  presence <- vector(mode="list", length=n)
  for(i in 1:n){
    presence[[i]] <- cluster_list %>%
      sapply(function(x){is.element(i, x)}) %>%
      which()
  }
  
  # Output results
  results <- list(
    overlap = list(presence=presence, cluster_list=cluster_list),
    presence = presence,
    cluster_list = cluster_list,
    cluster_coords = cluster_coords
  )
  return(results)
}



# Deprecated Functions ----------------------------------------------------

#' Title
#'
#' @param geo 
#' @param cluster.center 
#' @param cluster.end 
#'
#' @return
#' @export
#'
#' @examples
circle <- function(geo, cluster.center, cluster.end){
  # Compute interpoint distance
  distance <- dist(as.matrix(geo), upper=TRUE, diag=TRUE)
  distance <- as.matrix(distance)[cluster.center, cluster.end]
  
  # For drawing radius of cluster on map
  polar <- seq(0, 2*pi, length=1000)
  
  cluster.radius <- data.frame(cbind(
    x=distance*cos(polar)+geo$x[cluster.center],
    y=distance*sin(polar)+geo$y[cluster.center]
  ))
  return(cluster.radius)
}



#' Title
#'
#' @param input 
#'
#' @return
#' @export
#'
#' @examples
grid2latlong <- function(input){
  toradians <- atan(1)/45
  radiusearth <- 0.5*(6378.2+6356.7)
  sine51 <- sin( 51.5*toradians )
  
  
  #-------------------------------------------------------------------------------
  # If a Spatial Polygons
  #-------------------------------------------------------------------------------
  if(is(input)[1] == "SpatialPolygons"){
    for(i in 1:length(input@polygons)){
      # for all Polygons's in polygon
      for(j in 1:length(input@polygons[[i]]@Polygons)){
        # Convert coordinates
        new.coords <- as.matrix(cbind(
          input@polygons[[i]]@Polygons[[j]]@coords[,1]/(toradians*radiusearth*sine51),
          input@polygons[[i]]@Polygons[[j]]@coords[,2]/(toradians*radiusearth)
        ))
        colnames(new.coords) <- NULL
        rownames(new.coords) <- NULL
        
        # Update Polygons
        input@polygons[[i]]@Polygons[[j]]@coords <- new.coords
        input@polygons[[i]]@Polygons[[j]] <- Polygon(input@polygons[[i]]@Polygons[[j]])
      }
      # Update polygons
      input@polygons[[i]] <- Polygons(
        input@polygons[[i]]@Polygons,
        ID=input@polygons[[i]]@ID
      )	
    }
    output <- SpatialPolygons(input@polygons,proj4string=CRS("+proj=utm"))
    
    
    #-------------------------------------------------------------------------------
    # else return numeric
    #-------------------------------------------------------------------------------
  }else{
    output <- data.frame(cbind(
      x=input[, 1]/(toradians*radiusearth*sine51),
      y=input[, 2]/(toradians*radiusearth)
    ))		
  }
  
  return(output)	
}



#' Title
#'
#' @param poly 
#' @param coordinate.system 
#' @param area.names 
#' @param nrepeats 
#'
#' @return
#' @export
#'
#' @examples
polygon2spatial_polygon <- function(poly, coordinate.system, area.names = NULL, nrepeats = NULL){
  #-------------------------------------------------------------------------------
  # Deal with non-specified values
  #-------------------------------------------------------------------------------
  if(missing(coordinate.system)){
    stop("Coordinate system must be specified: '+proj=utm' or '+proj=longlat'.")
  }
  
  if(is.null(nrepeats)){
    nrepeats <- rep(1, sum(is.na(poly[,1]))+1 )	
  }
  
  if(is.null(area.names)){
    area.names <- as.character( 1:length(nrepeats) )
  }
  
  
  #-------------------------------------------------------------------------------
  # Create list of all polygon objects
  #-------------------------------------------------------------------------------
  na.index <- which(is.na(poly[,1]))
  n <- length(nrepeats)
  list.polygon <- NULL
  
  # First Case
  list.polygon <-	list(Polygon(poly[1:(na.index[1]-1),], hole=FALSE))					
  
  # Middle cases
  for(i in 1:(length(na.index)-1)){
    list.polygon <- c(list.polygon,list(Polygon(
      poly[(na.index[i]+1):(na.index[i+1]-1),], hole=FALSE)))	
  }
  
  # Last case
  list.polygon <-	c(list.polygon,list(Polygon(
    poly[(na.index[i+1]+1):length(poly[,1]),], hole=FALSE)
  ))
  
  
  #-------------------------------------------------------------------------------
  # From list of polygon objects, create "polygon" objects, that has one element
  # for each county.  A county can consist of several polygon as indicated by
  # nrepeats
  #-------------------------------------------------------------------------------
  list.polygons <- NULL
  
  start <- 1
  for( i in 1:length(nrepeats) ){
    end <- start + nrepeats[i] - 1
    
    temp.polygon <- NULL
    for(j in start:end){
      temp.polygon <- c(temp.polygon, list(list.polygon[[j]]))
    }
    
    list.polygons <- c(list.polygons, list(
      Polygons(temp.polygon, ID=area.names[i])
    ))
    start <- end + 1	
  }
  
  
  #-------------------------------------------------------------------------------
  # Output spatial polygons object
  #-------------------------------------------------------------------------------
  Spatial.Polygon <- 
    SpatialPolygons(list.polygons, proj4string=CRS(coordinate.system))
  
  return(Spatial.Polygon)
}



#' Title
#'
#' @param input 
#'
#' @return
#' @export
#'
#' @examples
latlong2grid <- function(input){
  toradians <- atan(1)/45
  radiusearth <- 0.5*(6378.2+6356.7)
  sine51 <- sin( 51.5*toradians )
  
  
  #-------------------------------------------------------------------------------
  # If a Spatial Polygon
  #-------------------------------------------------------------------------------
  if(is(input)[1] == "SpatialPolygons"){  
    for( i in 1:length(input@polygons) ){
      # for all Polygons's in polygon
      for( j in 1:length(input@polygons[[i]]@Polygons) ){
        # Convert coordinates
        new.coords <- cbind(
          (input@polygons[[i]]@Polygons[[j]]@coords[,1]*toradians)*radiusearth*sine51,
          (input@polygons[[i]]@Polygons[[j]]@coords[,2]*toradians)*radiusearth
        )
        new.coords <- as.matrix(new.coords)
        colnames(new.coords) <- NULL
        rownames(new.coords) <- NULL
        
        # Update Polygons
        input@polygons[[i]]@Polygons[[j]]@coords <- new.coords
        input@polygons[[i]]@Polygons[[j]] <- Polygon(input@polygons[[i]]@Polygons[[j]])
      }
      # Update polygons
      input@polygons[[i]] <- Polygons(
        input@polygons[[i]]@Polygons,
        ID=input@polygons[[i]]@ID
      )	
    }
    
    output <- SpatialPolygons(input@polygons,proj4string=CRS("+proj=utm"))
    
    
    #-------------------------------------------------------------------------------  
    # else return numeric
    #-------------------------------------------------------------------------------
  }else{
    output <- data.frame(cbind(
      x=(input[,1]*toradians)*radiusearth*sine51,
      y=(input[,2]*toradians)*radiusearth
    ))		
  }
  
  return(output)	
}



#' Title
#'
#' @param geo 
#' @param population 
#' @param pop.upper.bound 
#'
#' @return
#' @export
#'
#' @examples
zones <- function(geo, population, pop.upper.bound){
  # number of areas
  n <- nrow(geo)
  # total population
  total.pop <- sum(population)
  # Interpoint distance matrix
  dist <- as.matrix(dist(as.matrix(geo), upper=TRUE, diag=TRUE))
  
  
  #-------------------------------------------------------------------------------
  # For each area, list of closest neighbors up until pop.upper.bound of
  # population is met.  We count the number of candidate zones in n.zones 
  # Note:  for each county, the list of nearest neighbors form the total number of
  # candidate zones
  #-------------------------------------------------------------------------------
  nearest.neighbors <- vector(mode="list", length=n)
  
  n.zones <- 0
  vector.cutoffs <- rep(0, n+1)
  
  for(i in 1:n) {	
    # Sort the areas by distance, then include them one-by-one until cluster bound
    # is met
    neighbors <- order(dist[,i])
    
    # include only up until pop.upper.bound is hit
    neighbors <- neighbors[ which( 
      cumsum(population[neighbors])/total.pop <= pop.upper.bound
    )] 
    
    nearest.neighbors[[i]] <- neighbors
    
    # Update total number of zones
    n.zones <- n.zones + length(neighbors)
    
    # internal
    vector.cutoffs[i+1] <- n.zones	
  }
  
  cluster.coords <- cbind(rep(1:n, times=diff(vector.cutoffs)), 
                          unlist(nearest.neighbors))
  
  
  #-------------------------------------------------------------------------------
  # Output results
  #-------------------------------------------------------------------------------
  results <- list(
    nearest.neighbors = nearest.neighbors, 
    cluster.coords = cluster.coords, 
    dist=dist)
  
  return(results)
}

