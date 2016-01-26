\name{kulldorffMC}
\alias{kulldorffMC}
\title{Compute permutation distribution for kulldorff method}
\description{Given a matrix of permuted cases, return the likelihood ratio test statistic for each permuation}
\usage{kulldorffMC(permutedCaseMatrix, expectedCases, nearestNeighbors, nZones, logLkhdType)}
\arguments{
  \item{permutedCaseMatrix}{matrix of \code{n.sim} rows of permuted cases}
  \item{expectedCases}{expected cases for each area}
  \item{nearestNeighbors}{see \code{\link{create_geo_objects}}}
  \item{nZones}{number of single zones}
  \item{logLkhdType}{type of lkhd to use: either "poisson" or "binomial"}
}
\value{Vector of length \code{n.sim} of LRT statistics}
\references{ 
 SatScan:  Software for the spatial, temporal, and space-time scan statistics \url{http://www.satscan.org/} Kulldorff, M. (1997) A spatial scan statistic. \emph{Communications in Statistics: Theory and Methods}, \bold{26}, 1481--1496.
 
 Kulldorff M. and Nagarwalla N. (1995) Spatial disease clusters: Detection and Inference.
\emph{Statistics in Medicine}, \bold{14}, 799--810.
}
\author{Albert Y. Kim}
\seealso{\code{\link{kulldorff}}}
\keyword{internal}
