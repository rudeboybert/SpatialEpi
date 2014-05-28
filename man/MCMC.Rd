\name{MCMC}
\alias{MCMC}
\title{MCMC Wrapper Function}
\description{Internal wrapper function for C call to MCMC method.}
\usage{
MCMC(n.sim, overlap, cluster.coords, J, values, lambda, theta.init)
}
\arguments{
  \item{n.sim}{number of MCMC simulations}
  \item{overlap}{output of \code{\link{create.geo.objects}}: a list with two elements: \code{1. presence} which lists for each areas all the single zones it is present in and \code{2. cluster_list} for each single zone its component areas}
  \item{cluster.coords}{output of \code{\link{create.geo.objects}}: \code{n.zones x 2} matrix of the centering and radial area of each single zone}
  \item{J}{Maximum number of clusters/anti-clusters to consider}
  \item{values}{values associated with each single zone to use in Metropolis-Hastings ratio}
  \item{lambda}{lambda from definition of prior on single zones}
  \item{theta.init}{initial configuration}
}
\value{
\item{sample}{MCMC chain}
\item{move_trace}{trace of proposed move at each iteration}
\item{accpt_trace}{trace of acceptance at each iteration}
}
\seealso{\code{\link{process.MCMC.chain}}}
\keyword{internal}
