#-------------------------------------------------------------------------
# zones():	
# Arguments:
# -geo:  n x 2 table of (x,y) geographic coordinates of the centroid of each 
#	location.  Coordinate system must be a grid (i.e. Latitude and Longitude are not)
# -area.population:	n-length vector of population for each location
# -pop.upper.bound: scalar in [0,1] indicating the maximum proportion of the total
#	population we should consider
#
# Outputs:
# -nearest.neighbors:  list of nearest neighbors for each location
# -n.zones:  total number of zones
# -cluster.coords:  n.zones x 2 table of the center and boundary location for each zone
# -dist:  distance matrix

`zones` <- function(geo, area.population, pop.upper.bound){
	# number of areas
	n <- nrow(geo)
	# total population
	pop.N <- sum(area.population)
	
	# Interpoint distance matrix
	dist <- as.matrix(dist(as.matrix(geo),upper=TRUE,diag=TRUE))

	# Empty list
	nearest.neighbors <- vector(mode="list",length=n)

	#--------------------------------------------------------
	# List of closest neighbors for each of n areas, up until pop.upper.bound of population
	# We count the number of candidate zones in n.zones
	# Note:  for each county, the list of nearest neighbors form the total number of candidate zones
	n.zones <- 0
	vector.cutoffs <- rep(0,n+1)	# internal
		
	for(i in 1:n){	
		# Sort the counties by distance, then include them one-by-one until cluster bound is hit
		neighbors <- order(dist[,i])
		# include only up until pop.upper.bound is hit
		neighbors <- neighbors[ which( cumsum(area.population[neighbors])/pop.N <= pop.upper.bound)] 

		nearest.neighbors[[i]] <- neighbors
		
		# Update total number of zones
		n.zones <- n.zones + length( neighbors )
		
		# internal
		vector.cutoffs[i+1] <- n.zones	
	}
	rm(neighbors)
	
	nearest.neighbors.vector <- unlist(nearest.neighbors)
	cluster.coords <- cbind( rep(1:n,times=diff(vector.cutoffs)), nearest.neighbors.vector)


	#--------------------------------------------------------
	# Output Results
	results <- list(nearest.neighbors, n.zones, cluster.coords, dist)
	names(results) <- c("nearest.neighbors", "n.zones", "cluster.coords", "dist") 
	return(results)
}
