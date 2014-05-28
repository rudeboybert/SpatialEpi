#-------------------------------------------------------------------------
# circle():
# for plotting purposes.  Based on geo data object, output:
# -cluster.radius:  geographical coordinates of circle of zone 
# -cluster:  an array of all areas in cluster.  it necessarily has 
#  cluster.center as the first element, and cluster.end as the last
#

`circle` <- function(geo, cluster.center, cluster.end){

	# Compute interpoint distance
	distance <- as.matrix(dist(as.matrix(geo),upper=TRUE,diag=TRUE))[cluster.center, cluster.end]

	# For drawing radius of cluster on map
	polar <- seq(0,2*pi, length=1000)
	
	cluster.radius <- data.frame(cbind(
		x=distance*cos(polar)+geo$x[cluster.center],
		y=distance*sin(polar)+geo$y[cluster.center]
	))
	rm(polar)

	

return(cluster.radius)
}
