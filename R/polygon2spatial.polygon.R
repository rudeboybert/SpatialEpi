`polygon2spatial.polygon` <-
function(poly, coordinate.system, area.names = NULL, nrepeats = NULL){

#-----------------------------------------------------------------
# 1. Deal with non-specified values
#-----------------------------------------------------------------
if(missing(coordinate.system)){
	stop("Coordinate system must be specified: '+proj=utm' or '+proj=longlat'.")
}

if(is.null(nrepeats)){
	nrepeats <- rep(1, sum(is.na(poly[,1]))+1 )	
}

if(is.null(area.names)){
	area.names <- as.character( 1:length(nrepeats) )
}



#----------------------------------------------------------------
# 2. create list of all polygon objects
#-----------------------------------------------------------------
na.index <- which(is.na(poly[,1]))
n <- length(nrepeats)
		
list.polygon <- NULL

# First Case
list.polygon <-	list(Polygon(
	poly[1:(na.index[1]-1),]
,hole=FALSE))					

# Middle cases
for(i in 1:(length(na.index)-1)){
	list.polygon <- c(list.polygon,list(Polygon(
		poly[(na.index[i]+1):(na.index[i+1]-1),]
	,hole=FALSE)))	
}

# Last case
list.polygon <-	c(list.polygon,list(Polygon(
	poly[(na.index[i+1]+1):length(poly[,1]),]
,hole=FALSE)))



#-----------------------------------------------------------------
# 3.  from list of polygon objects, create "polygon" objects, that has
# one element for each county.  A county can consist of several polygon
# as indicated by nrepeats
#-----------------------------------------------------------------
list.polygons <- NULL

start <- 1
for( i in 1:length(nrepeats) ){
	end <- start + nrepeats[i] - 1
	
	temp.polygon <- NULL
	for(j in start:end){
		temp.polygon <- c(temp.polygon,list(list.polygon[[j]]))
	}
	
	list.polygons <- c(list.polygons, list(
		Polygons(temp.polygon,ID=area.names[i])
	))
	start <- end + 1	
}



#-----------------------------------------------------------------
# 4.  Output Spatial Polygons object
#-----------------------------------------------------------------
Spatial.Polygon <- SpatialPolygons(list.polygons,proj4string=CRS(coordinate.system))
Spatial.Polygon
}

