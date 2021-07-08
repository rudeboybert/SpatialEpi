\name{polygon2spatial_polygon}
\Rdversion{1.1}
\alias{polygon2spatial_polygon}
\title{
Convert a Polygon to a Spatial Polygons Object
}
\description{
Converts a polygon (a matrix of coordinates with NA values to separate subpolygons) into a Spatial Polygons object.
}
\usage{
polygon2spatial_polygon(poly, coordinate.system, area.names = NULL, nrepeats = NULL)
}
\arguments{
  \item{poly}{a 2-column matrix of coordinates, where each complete subpolygon is separated by NA's}
  \item{coordinate.system}{the coordinate system to use}
  \item{area.names}{names of all areas}
  \item{nrepeats}{number of subpolygons for each area}
}
\details{
Just as when plotting with the \code{\link[graphics]{polygon}} function, it is assumed that each subpolygon is to be closed by joining the last point to the first point.  In the matrix \code{poly}, NA values separate complete subpolygons. 

\code{coordinate.system} must be either \code{'+proj=utm'} or \code{'+proj=longlat'}.

In the case with an area consists of more than one separate closed polygon, \code{
nrepeats} specifies the number of closed polygons associated with each area.
}
\value{
An object of class SpatialPolygons (See \link[sp]{SpatialPolygons-class} from the \pkg{sp} package).
}
\references{
Bivand, R. S., Pebesma E. J., and Gomez-Rubio V. (2008) \emph{Applied Spatial Data Analysis with R}.  Springer Series in Statistics.

E. J. Pebesma and R. S. Bivand. (2005) Classes and methods for spatial data in R. \emph{R News}, \bold{5}, 9--13.  
}
\author{Albert Y. Kim}
%\seealso{}
\examples{
data(scotland)

polygon <- scotland$polygon$polygon
coord.system <- '+proj=utm'
names <- scotland$data$county.names
nrepeats <- scotland$polygon$nrepeats

spatial.polygon <- polygon2spatial_polygon(polygon,coord.system,names,nrepeats)

par(mfrow=c(1,2))
# plot using polygon function
plot(polygon,type='n',xlab="Eastings (km)",ylab="Northings (km)",main="Polygon File")
polygon(polygon)

# plot as spatial polygon object
plot(spatial.polygon,axes=TRUE)
title(xlab="Eastings (km)",ylab="Northings (km)",main="Spatial Polygon")

# Note that area 23 (argyll-bute) consists of 8 separate polygons
nrepeats[23]
plot(spatial.polygon[23],add=TRUE,col="red")
}
\keyword{file}
