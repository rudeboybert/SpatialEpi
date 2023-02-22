
globalVariables(c(
  "is","Polygon","Polygons","SpatialPolygons","CRS"
))



#' Convert Coordinates from Grid to Latitude/Longitude
#'
#'
#' @description Convert geographic coordinates from Universal Transverse Mercator system to Latitude/Longitude.
#' 
#' @param input A data frame with columns named `x` and `y` of the UTM coordinates to convert or  an `n x 2` matrix of grid coordinates or an object of class SpatialPolygons (See [SpatialPolygons-class][sp::SpatialPolygons-class])
#'
#' @author Lance A. Waller
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
#'   radius.  Perhaps for U.S., might use an "average" latitude, 30 deg
#'   is roughly Houston, 49deg is most of N bdry of continental 48
#'   states.  0.5(30+49)=39.5 deg.  so use r approx 6378.2*sin(51.5)
#' 
#' @details Longitude/latitudes are not a grid-based coordinate system:  latitudes are equidistant but the distance between longitudes varies.
#' @return
#' Either a data frame with the corresponding longitude and latitude, or a SpatialPolygons object with the coordinates changed.
#' @export
#'
#' @examples
#' coord <- data.frame(rbind(
#' # Montreal, QC
#' c(-6414.30, 5052.849),
#' # Vancouver, BC
#' c(-122.6042, 45.6605)
#' ))
#'
#' grid2latlong(coord)
#'
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
  output <- SpatialPolygons(input@polygons,proj4string=CRS("+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))

  
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
