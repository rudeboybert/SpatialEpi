\name{estimate_lambda}
\alias{estimate_lambda}
\title{Estimate q values for prior prob of j clusters/anti-clusters}
\description{Internal function to estimate values}
\usage{estimate_lambda(n.sim, J, prior.z, overlap, pi0)}
\arguments{
  \item{n.sim}{number of importance sampling iterations}
  \item{J}{maximum number of clusters/anti-clusters to consider}
  \item{prior.z}{prior probability of each single zone}
  \item{overlap}{output of \code{\link{create_geo_objects}}: a list with two elements: \code{1. presence} which lists for each areas all the single zones it is present in and \code{2. cluster_list} for each single zone its component areas}
  \item{pi0}{prior probability of no clusters}
}
\value{estimates of lambda and prior.j}
\references{Wakefield and Kim}
\seealso{\code{\link{estimate_lambda}}}
\keyword{internal}