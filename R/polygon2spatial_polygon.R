

globalVariables(c(
  "Polygon","Polygons","SpatialPolygons","CRS"
))




#' Convert a Polygon to a Spatial Polygons Object
#'
#' @description Converts a polygon (a matrix of coordinates with NA values to separate subpolygons) into a Spatial Polygons object.
#' @param poly a 2-column matrix of coordinates, where each complete subpolygon is separated by NA's
#' @param coordinate.system the coordinate system to use
#' @param area.names names of all areas
#' @param nrepeats number of sub polygons for each area
#'
#' @details Just as when plotting with the \code{\link[graphics]{polygon}} function, it is assumed that each subpolygon is to be closed by joining the last point to the first point.  In the matrix \code{poly}, NA values separate complete subpolygons. 
#' \code{coordinate.system} must be either \code{'+proj=utm'} or \code{'+proj=longlat'}.
#' In the case with an area consists of more than one separate closed polygon, \code{nrepeats} specifies the number of closed polygons associated with each area.
#' 
#' @references Bivand, R. S., Pebesma E. J., and Gomez-Rubio V. (2008) \emph{Applied Spatial Data Analysis with R}.  Springer Series in Statistics. E. J. Pebesma and R. S. Bivand. (2005) Classes and methods for spatial data in R. \emph{R News}, \bold{5}, 9--13.  
#' @author Albert Y. Kim
#' @return
#' An object of class SpatialPolygons (See \link[sp]{SpatialPolygons-class} from the \pkg{sp} package).
#' @export
#'
#' @examples
#' 
#' data(scotland)
#' 
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
