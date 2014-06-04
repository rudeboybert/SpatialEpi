\name{coeff}
\alias{coeff}
\title{Compute log Bayes Factors}
\description{Compute log Bayes Factors for each single zone based on observed and expected counts from \code{n} areas}
\usage{coeff(y_vector, E_vector, a_values, b_values, cluster_list)}
\arguments{
  \item{y_vector}{Vector of length n of cases}
  \item{E_vector}{Vector of length n of expected cases}
  \item{a_values}{Vector of length 2 of shape parameters (wide, narrow)}
  \item{b_values}{Vector of length 2 of rate parameters}
  \item{cluster_list}{Output of \code{\link{create_geo_objects}}: list of length \code{n.zones} listing, for each single zone, its component areas}
}
\value{Vector of length \code{n.zones} with log Bayes Factor for each single zone.  }
\references{Wakefield J. and Kim A.Y. (2013) A Bayesian model for cluster detection. \emph{Biostatistics}, \bold{14}, 752--765.}\author{Albert Y. Kim}
\keyword{interal}
