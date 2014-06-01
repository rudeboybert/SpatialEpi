create_geo_objects <-
function(max.prop, pop, centroids, sp.obj, area=NULL){
  # Number of areas
  n <- nrow(centroids)
  
  #-----------------------------------------------------------------------------
  # Adjacency Stuff
  #-----------------------------------------------------------------------------
  # Compute adjacency matrix
  adj <- poly2nb(sp.obj, queen=TRUE)
  adj <- nb2mat(adj)[1:n,1:n]
  
  # Convert Adjacency Matrix To A List
  adj.new <- vector("list",length=n)
  for(i in 1:n){
    adj.new[[i]] <- which(adj[i,]!=0)
  }
  adj <- adj.new
  rm(adj.new)
  
  #-----------------------------------------------------------------------------
  # Set up Zones
  #-----------------------------------------------------------------------------
  if(max.prop < max(pop/sum(pop))){
    print(paste("max.prop needs to be at least ", 
                round(max(pop/sum(pop)),3),
                sep="")
          )
  }
  
  geoInfo <- zones(centroids, pop, max.prop)
  nn <- geoInfo$nearest.neighbors
  n.zones <- geoInfo$n.zones
  cluster.coords <- geoInfo$cluster.coords
  
  # Define areas each cluster
  cluster.list <- vector(mode="list",length=n.zones)
  counter <- 1
  for(i in 1:n)
    for(j in 1:length(nn[[i]])){
      cluster.list[[counter]] <- nn[[i]][1:j] 
      counter <- counter + 1  
    } 
  
  rm(counter)
  
  #-----------------------------------------------------------------------------
  # Zone Based Counts
  #-----------------------------------------------------------------------------
  # If not provided areas
  if(is.null(area)){
    area <- rep(0,n)
    sp.obj.grid <- latlong2grid(sp.obj)
    for(i in 1:n)
      area[i] <- sp.obj.grid@polygons[[i]]@area
    rm(sp.obj.grid)
  }
  areaz <- sapply(cluster.list,function(x){sum(area[x])})
  
  #-----------------------------------------------------------------------------
  # Generate Overlap
  #-----------------------------------------------------------------------------
  # For each area, get "presense":  list all unit zones that include it
  presence <- vector(mode="list",length=n)
  for(i in 1:n){
    presence[[i]] <- which( sapply(cluster.list,function(x){is.element(i,x)}) )
  }
  
  # If buffer wanted and regardless of overlap matrix or list, update presence
  # list need to preserve original presence list all the way thru, so
  # temporarily store results here
  presence_temp <- presence
  
  # Loop thru all areas
  for(i in 1:n){
    current.zones <- presence[[i]]
    adjacent.areas <- adj[[i]]
    adjacent.zones <- unique(unlist(presence[adjacent.areas]))
    
    presence_temp[[i]] <- sort(unique(c(adjacent.zones,current.zones)))
  }
  
  # restore presence list
  presence <- presence_temp
  rm(presence_temp, current.zones, adjacent.areas, adjacent.zones)
  
  #----------------------------------------------------------------
  # If data set is large, set overlap as a list
  #----------------------------------------------------------------
  overlap <- list(presence=presence, cluster.list=cluster.list)
  return(list(overlap=overlap, cluster.coords=cluster.coords, areaz=areaz))
}
