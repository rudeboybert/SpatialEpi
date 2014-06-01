\name{kulldorffMC}
\alias{kulldorffMC}
\title{Compute permutation distribution for kulldorff method}
\description{
Given a matrix of permuted cases, return for each permuation the LRT statistic
}
\usage{
kulldorffMC(permutedCaseMatrix, expectedCases, nearestNeighbors, nZones, logLkhdType)
}
\arguments{
  \item{permutedCaseMatrix}{Matrix of permuted cases}
  \item{expectedCases}{Expected cases for each area}
  \item{nearestNeighbors}{See \code{\link{create_geo_objects}}}
  \item{nZones}{number of single zones}
  \item{logLkhdType}{Type of lkhd to use:  either "poisson" or "binomial"}
}
\value{Vector of length \code{n.sim} of LRT statistics}
\author{Albert Y. Kim}\seealso{
\code{\link{kulldorff}}
}
\keyword{internal}
