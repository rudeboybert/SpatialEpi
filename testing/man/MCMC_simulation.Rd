\name{MCMC_simulation}
\alias{MCMC_simulation}
\title{MCMC simulation to sample configurations}
\description{MCMC simulation to estimate prior and posterior quantities by sampling 
configurations.}
\usage{MCMC_simulation(n_sim, pattern, theta_init, overlap, cluster_coords, 
p_moves_orig, J, lkhd_z, lambda)
}
\arguments{
  \item{n_sim}{number of MCMC iterations}
  \item{pattern}{alternating pattern between unif and prop prior on single zones}
  \item{theta_init}{initial configuration}
  \item{overlap}{output of \code{\link{create_geo_objects}}: a list with two elements: \code{1. presence} which lists for each areas all the single zones it is present in and \code{2. cluster_list} for each single zone its component areas}
  \item{cluster_coords}{output of \code{\link{create_geo_objects}}: \code{n.zones x 2} matrix of the centering and radial area of each single zone}
  \item{p_moves_orig}{probability of sampling each of the 5 possible moves to explore sample space}
  \item{J}{maximum number of clusters/anti-clusters to consider}
  \item{lkhd_z}{values associated with each single zone to use in Metropolis-Hastings ratio}
  \item{lambda}{lambda from definition of prior on single zones}
}
\value{
  \item{sample}{sampled configurations}
  \item{move_trace}{trace of moves (1 = growth, 2 = trim, 3 = recenter, 4 = death, 5 = birth)}
  \item{accpt_trace}{trace of acceptance (0 = not accepted)}
  \item{ratio_trace}{trace of Metropolis-Hastings ratio}  
}
\references{Wakefield J. and Kim A.Y. (2013) A Bayesian model for cluster detection. \emph{Biostatistics}, \bold{14}, 752--765.}
\author{Albert Y. Kim}
\seealso{\code{\link{create_geo_objects}}, \code{\link{process_MCMC_sample}}}
%\examples{}
\keyword{internal}
