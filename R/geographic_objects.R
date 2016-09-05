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
#' in the study region.
#' @param population vector of length n of the corresponding population 
#' counts.
#' @param pop_upper_bound upper bound of the proportion of the study 
#' region's population the single zones can contain.
#'   
#' @return A list with two objects
#' \item{nearest_neighbors}{A list of length n of the single zones for each area} 
#' \item{cluster_coords}{A data frame with 2 columns of the centering 
#' and radial area for each single zone.}
#' @references Kulldorff, M. (1997) A spatial scan statistic. 
#' \emph{Communications in Statistics: Theory and Methods}, \bold{26}, 1481--1496.
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



#' Create geographic objects needed for Bayesian method
#'
#' @inheritParams define_single_zones
#'
#' @return A list with three objects
#' \item{presence}{List of length n, for each area, the single zones it is present in.}
#' \item{cluster_list}{List of length n_zones indicating the areas included in each single zone}
#' \item{cluster_coords}{A data frame with 2 columns of the centering and radial area for each single zone.}
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


#' Compute cartesian coordinates of a cluster center and radius
#' 
#' This function is used for plotting purposes.
#'
#' @param geo A \code{n x 2} table of the x-coordinate and y-coordinates of the 
#' centroids of each area
#' @param cluster.center The area index (an integer between \code{1} and 
#' \code{n}) indicating the center of the circle
#' @param cluster.end The area index (an integer between \code{1} and \code{n}) 
#' indicating the area at the end of the circle
#'
#' @return blah
#' \item{cluster.radius}{A data frame that you can plot}
#' @export
#'
#' @examples
#' data(pennLC)
#' geo <- pennLC$geo[,2:3]
#' plot(geo, type='n')
#' text(geo, labels=1:nrow(geo))
#' lines(circle(geo, 23, 46), col = "red") 
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








#' Convert a Polygon to a Spatial Polygons Object
#' 
#' Converts a polygon (a matrix of coordinates with NA values to separate 
#' subpolygons) into a Spatial Polygons object. Just as when plotting with the 
#' \code{\link[graphics]{polygon}} function, it is assumed that each subpolygon 
#' is to be closed by joining the last point to the first point.  In the matrix 
#' \code{poly}, NA values separate complete subpolygons. 
#'
#' @param poly A 2-column matrix of coordinates, where each complete subpolygon 
#' is separated by NA's
#' @param coordinate.system the coordinate system to use. Must be either 
#' \code{'+proj=utm'} or \code{'+proj=longlat'}.  In the case with an area 
#' consists of more than one separate closed polygon, \code{nrepeats} specifies 
#' the number of closed polygons associated with each area.
#' @param area.names names of all areas
#' @param nrepeats number of subpolygons for each area
#' @references Bivand, R. S., Pebesma E. J., and Gomez-Rubio V. (2008) \emph{Applied Spatial Data Analysis with R}.  Springer Series in Statistics.
#' @references E. J. Pebesma and R. S. Bivand. (2005) Classes and methods for spatial data in R. \emph{R News}, \bold{5}, 9--13.  
#' @return An object of class SpatialPolygons (See \link[sp]{SpatialPolygons-class} from the \pkg{sp} package).
#' @export
#' @examples
#' data(scotland)
#' polygon <- scotland$polygon$polygon
#' coord.system <- '+proj=utm'
#' names <- scotland$data$county.names
#' nrepeats <- scotland$polygon$nrepeats
#' 
#' spatial.polygon <- polygon2spatial_polygon(polygon,coord.system,names,nrepeats)
#' 
#' par(mfrow=c(1,2))
#' # plot using polygon function
#' plot(polygon,type='n',xlab="Eastings (km)",ylab="Northings (km)",main="Polygon File")
#' polygon(polygon)
#' 
#' # plot as spatial polygon object
#' plot(spatial.polygon,axes=TRUE)
#' title(xlab="Eastings (km)",ylab="Northings (km)",main="Spatial Polygon")
#' 
#' # Note that area 23 (argyll-bute) consists of 8 separate polygons
#' nrepeats[23]
#' plot(spatial.polygon[23],add=TRUE,col="red")
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





#' Convert Coordinates from Grid to Latitude/Longitude
#' 
#' Convert geographic coordinates from Universal Transverse Mercator system to 
#' Latitude/Longitude. Longitude/latitudes are not a grid-based coordinate 
#' system:  latitudes are equidistant but the distance between longitudes varies.
#'
#' @param input A data frame with columns named \code{x} and \code{y} of the 
#' UTM coordinates to convert
#'
#' @return Either a data frame with the corresponding longitude and latitude, or 
#' a SpatialPolygons object with the coordinates changed.
#' @export
#' @note Rough conversion of US lat/long to km (used by GeoBUGS):  (see also 
#' forum.swarthmore.edu/dr.math/problems/longandlat.html).  Radius of earth: 
#' r = 3963.34 (equatorial) or 3949.99 (polar) mi = 6378.2 or 6356.7 km, which 
#' implies: km per mile  = 1.609299 or 1.609295 a change of 1 degree of latitude 
#' corresponds to the same number of km, regardless of longitude.  
#' arclength=r*theta, so the multiplier for coord\$y should probably be just the 
#' radius of earth. On the other hand, a change of 1 degree in longitude 
#' corresponds to a different distance, depending on latitude.  (at N pole, the 
#' change is essentially 0.  at the equator, use equatorial radius.
#' 
#' @author Lance A. Waller
#' @seealso \code{\link{latlong2grid}}
#' @examples
#' coord <- data.frame(rbind(
#' # Montreal, QC
#' c(-6414.30, 5052.849),
#' # Vancouver, BC
#' c(-122.6042, 45.6605)
#' ))
#' grid2latlong(coord)
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




#' Convert Coordinates from Latitude/Longitude to Grid
#' 
#' Convert geographic latitude/longitude coordinates to kilometer-based grid 
#' coordinates.  Longitude/latitudes are not a grid-based coordinate system:  
#' latitudes are equidistant but the distance between longitudes varies.
#'
#' @param input either an \code{n x 2} matrix of longitude and latitude 
#' coordinates in decimal format or an object of class SpatialPolygons (See 
#' \link[sp]{SpatialPolygons-class})
#'
#' @return Either a data frame with the corresponding (x,y) kilometer-based grid 
#' coordinates, or a SpatialPolygons object with the coordinates changed.
#' @author Lance A. Waller
#' @seealso \code{\link{grid2latlong}}
#' @export
#' @note Rough conversion of US lat/long to km (used by GeoBUGS):  (see also 
#' forum.swarthmore.edu/dr.math/problems/longandlat.html).  Radius of earth: 
#' r = 3963.34 (equatorial) or 3949.99 (polar) mi = 6378.2 or 6356.7 km, which 
#' implies: km per mile  = 1.609299 or 1.609295 a change of 1 degree of latitude 
#' corresponds to the same number of km, regardless of longitude.  
#' arclength=r*theta, so the multiplier for coord\$y should probably be just the 
#' radius of earth. On the other hand, a change of 1 degree in longitude 
#' corresponds to a different distance, depending on latitude.  (at N pole, the 
#' change is essentially 0.  at the equator, use equatorial radius.

#' @examples
#' # Convert coordinates
#' coord <- data.frame(rbind(
#' # Montreal, QC:  Latitude: 45deg 28' 0" N (deg min sec), Longitude: 73deg 45' 0" W
#' c(-73.7500, 45.4667),
#' # Vancouver, BC:  Latitude: 45deg 39' 38" N (deg min sec), Longitude: 122deg 36' 15" W
#' c(-122.6042, 45.6605)
#' ))
#' latlong2grid(coord)
#'
#' # Convert SpatialPolygon
#' data(pennLC)
#' new <- latlong2grid(pennLC$spatial.polygon)
#' par(mfrow=c(1,2))
#' plot(pennLC$spatial.polygon,axes=TRUE)
#' title("Lat/Long")
#' plot(new,axes=TRUE)
#' title("Grid (in km)")
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






#' Create set of all single zones and output geographical information
#' 
#' Based on the population counts and centroid coordinates of each of \code{n} 
#' areas, output the set of \code{n.zones} single zones as defined by Kulldorff 
#' and other geographical information.
#'
#' @param geo \code{n x 2} table of the (x,y)-coordinates of the area centroids
#' @param population a vector of population counts of each area
#' @param pop.upper.bound maximum proportion of study region each zone can contain
#'
#' @return  A list containing
#' \item{nearest.neighbors}{list of \code{n} elements, where each element is a 
#' vector of the nearest neighbors in order of distance up until 
#' \code{pop.upper.bound} of the total population is attained}
#' \item{cluster.coords}{\code{n.zones x 2} table of the center and the radial 
#' area for each zone}
#' \item{dist}{\code{n x n} inter-point distance matrix of the centroids}
#' @export
#' @references Kulldorff, M. (1997) A spatial scan statistic. 
#' \emph{Communications in Statistics: Theory and Methods}, \bold{26}, 1481--1496.
#' @references Kulldorff M. and Nagarwalla N. (1995) Spatial disease clusters: 
#' Detection and Inference. \emph{Statistics in Medicine}, \bold{14}, 799--810.
#' @examples
#' data(pennLC)
#' geo <- pennLC$geo[,2:3]
#' geo <- latlong2grid(geo)
#' population <- tapply(pennLC$data$population, pennLC$data$county, sum)
#' pop.upper.bound <- 0.5
#' geo.info <- zones(geo, population, pop.upper.bound)
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

