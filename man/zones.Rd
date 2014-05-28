\name{zones}
\alias{zones}
\title{Create Set of All Zones and Output Geographical Information}
\description{
Based on the population counts and centroid coordinates of each area, output the set of zones as defined by Kulldorff and other geographical information.
}
\usage{
zones(geo, area.population, pop.upper.bound)
}
\arguments{
  \item{geo}{an \code{n x 2} table of the (x,y)-coordinates of the area centroids}
  \item{area.population}{a vector of population counts of each area}
  \item{pop.upper.bound}{maximum proportion of study region each zone can contain}
}
\value{
A list containing
  \item{nearest.neighbors}{a list of \code{n} elements, where each element is a vector of the nearest neighbors in order of distance up until \code{pop.upper.bound} of the total population is attained}
  \item{n.zones}{total number of zones}
  \item{cluster.coords}{an \code{n x 2} table of the center and the radial areas for each zone}
  \item{dist}{an \code{n x n} inter-point distance matrix}
}
\references{
  Kulldorff, M. (1997)
 A spatial scan statistic.
  \emph{Communications in Statistics: Theory and Methods}, \bold{26}, 1481--1496. %\url{http://www.satscan.org/papers/k-cstm1997.pdf}

  Kulldorff M. and Nagarwalla N. (1995) Spatial disease clusters: Detection and Inference.
\emph{Statistics in Medicine}, \bold{14}, 799--810.  %\url{http://www.satscan.org/papers/k-sm1995.pdf}
}
\author{Albert Y. Kim}
\examples{
data(pennLC)
geo <- pennLC$geo[,2:3]
geo <- latlong2grid(geo)
area.population <- tapply(pennLC$data$population, pennLC$data$county, sum)
pop.upper.bound <- 0.5

geo.info <- zones(geo, area.population, pop.upper.bound)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ file }
