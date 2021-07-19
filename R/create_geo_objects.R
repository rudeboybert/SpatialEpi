

globalVariables(c(
  "poly2nb", "nb2mat"
))



#' Create geographical objects to be used in Bayesian Cluster Detection Method
#'
#' @description This internal function creates the geographical objects needed to run the Bayesian cluster detection method in \code{\link{bayes_cluster}}.  Specifically it creates all single zones based data objects, where single zones are the \emph{zones} defined by Kulldorff (1997).
#' 
#' @references Wakefield J. and Kim A.Y. (2013) A Bayesian model for cluster detection.\emph{Biostatistics}, \bold{14}, 752--765.
#' @author Albert Y. Kim
#' 
#' @param max.prop maximum proportion of study region's population each single zone can contain
#' @param population vector of length \code{n} of the population of each area
#' @param centroids \code{n x 2} table of the (x,y)-coordinates of the area centroids.  The coordinate system must be grid-based
#' @param sp.obj object of class SpatialPolygons (See \link[sp]{SpatialPolygons-class}) representing the study region
#'
#' @return
#' \item{overlap}{list with two elements: \code{1. presence} which lists for each area all the single zones it is present in and \code{2. cluster.list} for each single zone its component areas}
#' \item{cluster.coords}{\code{n.zones x 2} matrix of the center and radial area of each single zone}
#' 
#' 
#' 
#' 
#' @export
#' @importFrom spdep poly2nb
#' @importFrom spdep nb2mat
#' 
#'
#' @examples
#' data(pennLC)
#' max.prop <- 0.15
#' population <- tapply(pennLC$data$population, pennLC$data$county, sum)
#' centroids <- latlong2grid(pennLC$geo[, 2:3])
#' sp.obj <- pennLC$spatial.polygon
#' output <- create_geo_objects(max.prop, population, centroids, sp.obj)
#' ## number of single zones
#' nrow(output$cluster.coords)
create_geo_objects <-
function(max.prop, population, centroids, sp.obj){

# Number of areas
n <- nrow(centroids)


#-----------------------------------------------------------------------------
# Set up single zones
#-----------------------------------------------------------------------------
if(max.prop < max(normalize(population))){
  print(paste("max.prop needs to be at least", max(normalize(population))))
}

# Get basic single zone info
geoInfo <- zones(centroids, population, max.prop)
nearest.neighbors <- geoInfo$nearest.neighbors
cluster.coords <- geoInfo$cluster.coords
n.zones <- nrow(cluster.coords)

# Create list of length n.zones indicating the component areas for each zone
cluster.list <- vector(mode="list", length=n.zones)
counter <- 1
for(i in 1:n) {
  nn <- nearest.neighbors[[i]]
  for(j in 1:length(nn)) {
    cluster.list[[counter]] <- nn[1:j] 
    counter <- counter + 1  
  } 
}

#-----------------------------------------------------------------------------
# Generate overlap object which tracks the overlap between single zones
#-----------------------------------------------------------------------------
# For each area, get "presense":  list all unit zones that it is included in
presence <- vector(mode="list",length=n)
for(i in 1:n){
  presence[[i]] <- which(sapply(cluster.list, function(x){is.element(i, x)}))
}


#-------------------------------------------------------------------------------
# Compute "buffer" of areas between single zones.
#-------------------------------------------------------------------------------
# Compute adjacency matrix
adj <- poly2nb(sp.obj, queen=TRUE)
adj <- nb2mat(adj, zero.policy=TRUE)[1:n, 1:n]

# Convert Adjacency Matrix To A List
adj.new <- vector("list",length=n)
for(i in 1:n){
  adj.new[[i]] <- which(adj[i,]!=0)
}
adj <- adj.new
rm(adj.new)

# Update presence list to incorporate buffer.  We need to preserve original 
# presence list all the way thru, so temporarily store results here
presence_temp <- presence

# Loop thru all areas
for(i in 1:n){
  current.zones <- presence[[i]]
  adjacent.areas <- adj[[i]]
  adjacent.zones <- unique(unlist(presence[adjacent.areas]))
  presence_temp[[i]] <- sort(unique(c(adjacent.zones, current.zones)))
}

# restore presence list
presence <- presence_temp


#----------------------------------------------------------------
# Return list
#----------------------------------------------------------------
overlap <- list(presence=presence, cluster.list=cluster.list)
return(list(
  overlap=overlap, 
  cluster.coords=cluster.coords)
)
}
