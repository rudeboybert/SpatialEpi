\name{create_geo_objects}
\alias{create_geo_objects}
\title{Create geographical objects to be used in Bayesian Cluster Detection Method}
\description{
This internal function creates the geographical objects needed to run the Bayesian cluster detection method in \code{\link{bayes_cluster}}.  Specifically it creates all single zones based data objects, where single zones are the \emph{zones} defined by Kulldorff (1997).  
}
\usage{create_geo_objects(max.prop, population, centroids, sp.obj)}
\arguments{
  \item{max.prop}{maximum proportion of study region's population each single zone can contain}
  \item{population}{vector of length \code{n} of the population of each area}
  \item{centroids}{\code{n x 2} table of the (x,y)-coordinates of the area centroids.  The coordinate system must be grid-based} 
  \item{sp.obj}{object of class SpatialPolygons (See \link[sp]{SpatialPolygons-class}) representing the study region} 
}
\value{
  \item{overlap}{list with two elements: \code{1. presence} which lists for each area all the single zones it is present in and \code{2. cluster.list} for each single zone its component areas}
  \item{cluster.coords}{\code{n.zones x 2} matrix of the center and radial area of each single zone}
}
\references{Wakefield J. and Kim A.Y. (2013) A Bayesian model for cluster detection.
\emph{Biostatistics}, \bold{14}, 752--765.}
\author{Albert Y. Kim}
\seealso{\code{\link{latlong2grid}}, \code{\link{zones}}}
\examples{
data(pennLC)
max.prop <- 0.15
population <- tapply(pennLC$data$population, pennLC$data$county, sum)
centroids <- latlong2grid(pennLC$geo[, 2:3])
sp.obj <- pennLC$spatial.polygon
output <- create_geo_objects(max.prop, population, centroids, sp.obj)

## number of single zones
nrow(output$cluster.coords)
}
\keyword{file}
