\name{besag_newell_internal}
\alias{besag_newell_internal}
\title{Internal Besag Newell method}
\description{Internal}
\usage{
besag_newell_internal(observedCases, expectedCases, nearestNeighborsList, 
nZones, k)
}
\arguments{
  \item{observedCases}{observed cases for each area}
  \item{expectedCases}{Expected cases for each area}
  \item{nearestNeighbors}{See \code{\link{create_geo_objects}}}
  \item{nZones}{number of single zones}
  \item{k}{max number of cases in each}
}
\references{Besag J. and Newell J. (1991)
The Detection of Clusters in Rare Diseases \emph{Journal of the Royal Statistical Society. Series A (Statistics in Society)}, \bold{154}, 143--155}
\author{Albert Y. Kim}
\seealso{\code{\link{besag_newell}}}
\keyword{internal}
