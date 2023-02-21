
globalVariables(c(
  "is","Polygon","Polygons","SpatialPolygons","CRS"
))



#' @title Convert Coordinates from Latitude/Longitude to Grid
#'
#' @description Convert geographic latitude/longitude coordinates to kilometer-based grid coordinates.
#' 
#' @param input either an `n x 2` matrix of longitude and latitude coordinates in decimal format or an object of class SpatialPolygons 
#'
#' @note Rough conversion of US lat/long to km (used by GeoBUGS):  (see
#'   also forum.swarthmore.edu/dr.math/problems/longandlat.html).
#'   Radius of earth: r = 3963.34 (equatorial) or 3949.99 (polar) mi =
#'   6378.2 or 6356.7 km, which implies: km per mile  = 1.609299 or
#'   1.609295 a change of 1 degree of latitude corresponds to the same
#'   number of km, regardless of longitude.  arclength=r*theta, so the
#'   multiplier for coord y should probably be just the radius of
#'   earth. On the other hand, a change of 1 degree in longitude
#'   corresponds to a different distance, depending on latitude.  (at N
#'   pole, the change is essentially 0.  at the equator, use equatorial
#'   radius.
#'
#' @details Longitude/latitudes are not a grid-based coordinate system:  latitudes are equidistant but the distance between longitudes varies.
#'
#' @author Lance A. Waller
#' 
#' @return Either a data frame with the corresponding (x,y) kilometer-based grid coordinates, or a SpatialPolygons object with the coordinates changed.
#' 
#' 
#' @export
#'
#' @examples
#' ## Convert coordinates
#' coord <- data.frame(rbind(
#'  # Montreal, QC:  Latitude: 45deg 28' 0" N (deg min sec), Longitude: 73deg 45' 0" W
#'  c(-73.7500, 45.4667),
#'  # Vancouver, BC:  Latitude: 45deg 39' 38" N (deg min sec), Longitude: 122deg 36' 15" W
#'  c(-122.6042, 45.6605)
#' ))
#' latlong2grid(coord)
#' ## Convert SpatialPolygon
#' data(pennLC)
#' new <- latlong2grid(pennLC$spatial.polygon)
#' par(mfrow=c(1,2))
#' plot(pennLC$spatial.polygon,axes=TRUE)
#' title("Lat/Long")
#' plot(new,axes=TRUE)
#' title("Grid (in km)")
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
  
  output <- SpatialPolygons(input@polygons,proj4string=CRS("+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))

  
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
