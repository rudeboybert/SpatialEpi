circle <-
  function(geo, cluster.center, cluster.end){
    
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





grid2latlong <-
  function(input){
    
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




polygon2spatial_polygon <-
  function(poly, coordinate.system, area.names = NULL, nrepeats = NULL){
    
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







latlong2grid <-
  function(input){
    
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

