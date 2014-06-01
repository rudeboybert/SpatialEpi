\name{MCMC_simulation}
\alias{MCMC_simulation}
\title{Perform main MCMC simulation to estimate posterior quantities}
\description{
Main simulation loop
}
\usage{
MCMC_simulation(n_sim, pattern, theta_init, overlap, cluster_coords, p_moves_orig,
K, lkhd_z, lambda)
}
\arguments{
  \item{n_sim}{Number of MCMC iterations}
  \item{pattern}{alternating pattern between unif and prop prior on single zones}
  \item{theta_init}{initial configuration}
  \item{overlap}{See \code{\link{create_geo_objects}}}
  \item{cluster_coords}{See \code{\link{create_geo_objects}}}
  \item{p_moves_orig}{Probability of each of the 5 possible moves}
  \item{K}{Maximum number of single zones in a configuration}
  \item{lkhd_z}{Unnormalized likelihood in question}
  \item{lambda}{Weights}
}
\value{
  \item{sample }{Sampled configurations}
  \item{move_trace }{trace of moves}
  \item{accpt_trace }{trace of acceptance}
  \item{ratio_trace }{trace of Metropolis-Hastings ratio}  
}
%\references{}
\author{Albert Y. Kim}
%\examples{}
\keyword{internal}
