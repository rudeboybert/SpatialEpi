\name{computeAllLogLkhd}
\alias{computeAllLogLkhd}
\title{Compute log likelihood for all single zones}
\description{For all n.zones single zones, compute Likelihood ratio test statistic}
\usage{
computeAllLogLkhd(observedCases, expectedCases, nearestNeighborsList, nZones, logLkhdType)
}
\arguments{
  \item{observedCases}{Vector of length n of cases}
  \item{expectedCases}{Vector of length n of expected cases}
  \item{nearestNeighborsList}{See \code{\link{create_geo_objects}}}
  \item{nZones}{number of single zones}
  \item{logLkhdType}{type of likelihood:  either "poisson" or "binomial"}
}
\value{Vector of length \code{nZones} with LRT for each single zone}
\author{Albert Y. Kim}
\seealso{
\code{\link{poissonLogLkhd}}, \code{\link{binomialLogLkhd}}
}
\keyword{internal}