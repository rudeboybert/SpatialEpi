\name{check_overlap}
\alias{check_overlap}
\title{Check if proposed configurations overlap}
\description{Given a matrix \code{config} of proposed configurations, using
the single zone information in \code{overlap}, see if each proposed configuration
is valid in that the single zones don't overlap with each other.}  
\usage{check_overlap(config, overlap)}
\arguments{
  \item{config}{\code{n} by \code{k} matrix of configurations where \code{n} is
  the number of proposed configurations and \code{k} is the configuration size
  (in terms of the number of single zones)}
  \item{overlap}{See \code{\link{create_geo_objects}}}
}
\value{A vector of length \code{n} indicating if each of the proposed configs
overlap (0 if they do overlap, 1 if they don't)}
\author{Albert Y. Kim}
%\examples{}
\keyword{internal}
