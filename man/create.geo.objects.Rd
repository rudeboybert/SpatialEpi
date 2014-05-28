\name{create.geo.objects}
\alias{create.geo.objects}
\title{Create geographical objects to be used in \code{bayes.cluster} function.}
\description{
This internal function creates the geographical objects needed to run the Bayesian cluster detection method in \code{\link{bayes.cluster}}.  Specifically it creates all single zones based data objects, where single zones are the \emph{zones} defined by Kulldorff (1997).  
}
\usage{create.geo.objects(max.prop, pop, centroids, sp.obj, area = NULL)}
\arguments{
  \item{max.prop}{maximum proportion of study region's population each single zone can contain}
  \item{pop}{vector of length \code{n} of the population of each area}
  \item{centroids}{an \code{n x 2} table of the (x,y)-coordinates of the area centroids.  The coordinate system must be grid-based} 
  \item{sp.obj}{an object of class SpatialPolygons (See \link[sp]{SpatialPolygons-class}) representing the study region} 
  \item{area}{geographic area of each of the \code{n} areas. If not specified, this function will compute them based on \code{sp.obj}}
}
\value{
\item{overlap}{a list with two elements: \code{1. presence} which lists for each areas all the single zones it is present in and \code{2. cluster_list} for each single zone its component areas}
\item{cluster.coords}{\code{n.zones x 2} matrix of the centering and radial area of each single zone}
\item{areaz}{vector of the geographic surface area of each single zone}
}
\references{
  Kulldorff, M. (1997) A spatial scan statistic. \emph{Communications in Statistics: Theory and Methods}, \bold{26}, 1481--1496. %\url{http://www.satscan.org/papers/k-cstm1997.pdf}
}
\author{Albert Y. Kim}
\seealso{\code{\link{latlong2grid}}, \code{\link{zones}}}
\examples{
data(pennLC)
max.prop <- 0.15
pop <- tapply(pennLC$data$population, pennLC$data$county, sum)
centroids <- latlong2grid(pennLC$geo[, 2:3])
map <- pennLC$spatial.polygon
output <- create.geo.objects(max.prop, pop, centroids, map)

## number of single zones
nrow(output$cluster.coords)
length(output$areaz)
}
\keyword{file}
