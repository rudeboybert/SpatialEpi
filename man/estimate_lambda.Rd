\name{estimate_lambda}
\alias{estimate_lambda}
\title{Estimate lambda values}
\description{Internal function to estimate values of lambda needed for \code{\link{MCMC_simulation}} and prior probability of k clusters/anti-clusters for k=0,...,J}
\usage{estimate_lambda(n.sim, J, prior.z, overlap, pi0, equal.non.null=FALSE)}
\arguments{
  \item{n.sim}{number of importance sampling iterations}
  \item{J}{maximum number of clusters/anti-clusters to consider}
  \item{prior.z}{prior probability of each single zone}
  \item{overlap}{output of \code{\link{create_geo_objects}}: list with two elements: \code{presence} which lists for each area all the single zones it is present in and \code{cluster_list} for each single zone its component areas}
  \item{pi0}{prior probability of no clusters}
  \item{equal.non.null}{Flag to set equal prior probability for each of 1 through J clusters/anti-clusters}
}
\value{estimates of lambda, prior.j, and q}
\references{Wakefield J. and Kim A.Y. (2013) A Bayesian model for cluster detection. \emph{Biostatistics}, \bold{14}, 752--765.}
\seealso{\code{\link{MCMC_simulation}}}
\keyword{internal}