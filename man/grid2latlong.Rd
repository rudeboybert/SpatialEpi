\name{grid2latlong}
\alias{grid2latlong}
\title{Convert Coordinates from Grid to Latitude/Longitude}
\description{
Convert geographic coordinates from Universal Transverse Mercator system to Latitude/Longitude.
}
\usage{
grid2latlong(input)
}
\arguments{
  \item{input}{ A data frame with columns named \code{x} and \code{y} of the UTM coordinates to convert }
  \item{input}{either an \code{n x 2} matrix of grid coordinates or an object of class SpatialPolygons (See \link[sp]{SpatialPolygons-class})}
}
\note{ 
Rough conversion of US lat/long to km (used by GeoBUGS):  (see also forum.swarthmore.edu/dr.math/problems/longandlat.html).  Radius of earth: r = 3963.34 (equatorial) or 3949.99 (polar) mi = 6378.2 or 6356.7 km, which implies: km per mile  = 1.609299 or 1.609295 a change of 1 degree of latitude corresponds to the same number of km, regardless of longitude.  arclength=r*theta, so the multiplier for coord\$y should probably be just the radius of earth.
On the other hand, a change of 1 degree in longitude corresponds to a different distance, depending on latitude.  (at N pole, the change is essentially 0.  at the equator, use equatorial radius.  Perhaps for U.S., might use an "average" latitude, 30 deg is roughly Houston, 49deg is most of N bdry of continental 48 states.  0.5(30+49)=39.5 deg.  so use r approx 6378.2*sin(51.5)
}
\details{
Longitude/latitudes are not a grid-based coordinate system:  latitudes are equidistant but the distance between longitudes varies.
}
\value{
Either a data frame with the corresponding longitude and latitude, or a SpatialPolygons object with the coordinates changed.
}
%\references{}
\author{Lance A. Waller}
\seealso{\code{\link{latlong2grid}}}
\examples{
coord <- data.frame(rbind(
# Montreal, QC
c(-6414.30, 5052.849),
# Vancouver, BC
c(-122.6042, 45.6605)
))

grid2latlong(coord)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
