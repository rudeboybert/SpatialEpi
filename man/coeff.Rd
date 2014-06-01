\name{coeff}
\alias{coeff}
\title{Compute Bayes Factors for each single zone}
\description{Compute BF}
\usage{coeff(y_vector, E_vector, a_values, b_values, cluster_list)}
\arguments{
  \item{y_vector}{Vector of length n (areas) of cases}
  \item{E_vector}{Vector of length n (areas) of expected cases}
  \item{a_values}{Vector of length 2 of shape parameters (wide, narrow)}
  \item{b_values}{Vector of length 2 of rate parameters}
  \item{cluster_list}{See \code{\link{create_geo_objects}}}
}
\value{
Vector of length \code{n.zones} with Bayes Factor for each single zone.  
}
%\references{}
\author{Albert Y. Kim}
%\seealso{}
%\examples{}
\keyword{interal}
