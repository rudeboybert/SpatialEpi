create_geo_objects <-
function(max.prop, population, centroids, sp.obj, area=NULL){

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

# If not provided, compute geographic surface area for each single zone.  
if(is.null(area)){
  area <- rep(0,n)
  sp.obj.grid <- latlong2grid(sp.obj)
  for(i in 1:n)
    area[i] <- sp.obj.grid@polygons[[i]]@area
  rm(sp.obj.grid)
}
areaz <- sapply(cluster.list,function(x){sum(area[x])})


#-----------------------------------------------------------------------------
# Generate overlap object which tracks the overlap between single zones
#-----------------------------------------------------------------------------
# For each area, get "presense":  list all unit zones that it is included in
presence <- vector(mode="list",length=n)
for(i in 1:n){
  presence[[i]] <- which(sapply(cluster.list, function(x){is.element(i, x)}))
}


#----------------------------------------------------------------
# Return list
#----------------------------------------------------------------
overlap <- list(presence=presence, cluster.list=cluster.list)
return(list(
  overlap=overlap, 
  cluster.coords=cluster.coords, 
  areaz=areaz)
)
}
