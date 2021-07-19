

globalVariables(c(
  "dist"
))


#' Compute cartesian coordinates of a cluster center and radius
#' 
#' @description This function is used for plotting purposes
#' 
#' @author Albert Y. Kim
#' 
#' @param geo A \code{n x 2} table of the x-coordinate and y-coordinates of the centroids of each area
#' @param cluster.center The area index (an integer between \code{1} and \code{n}) indicating the center of the circle
#' @param cluster.end The area index (an integer between \code{1} and \code{n}) indicating the area at the end of the circle
#'
#' @return
#'  \item{cluster.radius}{A data frame that you can plot}
#' @export
#'
#' @examples data(pennLC)
#' geo <- pennLC$geo[,2:3]
#' plot(geo,type='n')
#' text(geo,labels=1:nrow(geo))
#' lines( circle(geo, 23, 46), col = "red" )
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
