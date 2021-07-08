\name{return_birth_moves}
\alias{return_birth_moves}
\title{Return all possible birth moves}
\description{Given a current configuration theta, return all possible
non-overlapping birth move configuration}
\usage{return_birth_moves(theta, overlap)}
\arguments{
  \item{theta}{current configuration of length k}
  \item{overlap}{see \code{\link{create_geo_objects}}}
}
\value{\code{NumericMatrix} with \code{k+1} rows and each column corresponding to a new possible
configuration}
\author{Albert Y. Kim}
\seealso{\code{\link{return_death_moves}}, \code{\link{return_local_moves}}}
\keyword{internal}
