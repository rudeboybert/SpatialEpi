#' Create set of all single zones and output geographical information
#'
#' @description Based on the population counts and centroid coordinates of each of `n` areas, output the set of `n.zones` single zones as defined by Kulldorff and other geographical information.
#' @param geo `n x 2` table of the (x,y)-coordinates of the area centroids
#' @param population a vector of population counts of each area
#' @param pop.upper.bound maximum proportion of study region each zone can contain
#'
#' @references Kulldorff, M. (1997) A spatial scan statistic. *Communications in Statistics: Theory and Methods*, **26**, 1481--1496.
#' Kulldorff M. and Nagarwalla N. (1995) Spatial disease clusters: Detection and Inference.
#' *Statistics in Medicine*, **14**, 799--810.
#' 
#' @author Albert Y. Kim
#' 
#' @return
#' A list containing
#' \item{nearest.neighbors}{list of `n` elements, where each element is a vector of the nearest neighbors in order of distance up until `pop.upper.bound` of the total population is attained}
#' \item{cluster.coords}{`n.zones x 2` table of the center and the radial area for each zone}
#' \item{dist}{`n x n` inter-point distance matrix of the centroids}
#'
#' 
#' @export
#' 
#' @examples
#' data(pennLC)
#' geo <- pennLC$geo[,2:3]
#' geo <- latlong2grid(geo)
#' population <- tapply(pennLC$data$population, pennLC$data$county, sum)
#' pop.upper.bound <- 0.5
#' geo.info <- zones(geo, population, pop.upper.bound)
zones <-
function(geo, population, pop.upper.bound){

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
