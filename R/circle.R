circle <-
function(geo, cluster.center, cluster.end){
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
