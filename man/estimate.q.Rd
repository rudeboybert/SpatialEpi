\name{estimate.q}
\alias{estimate.q}
\title{Estimate q values for prior prob of j clusters/anti-clusters}
\description{Internal function to estimate values}
\usage{estimate.q(n.sim, J, prior.z, overlap)}
\arguments{
  \item{n.sim}{number of importance sampling iterations}
  \item{J}{maximum number of clusters/anti-clusters to consider}
  \item{prior.z}{prior probability of each single zone}
  \item{overlap}{output of \code{\link{create.geo.objects}}: a list with two elements: \code{1. presence} which lists for each areas all the single zones it is present in and \code{2. cluster_list} for each single zone its component areas}
}
\value{estimates of q}
\references{Wakefield and Kim}
\seealso{\code{\link{estimate.lambda}}}
\keyword{internal}