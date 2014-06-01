\name{return_local_moves}
\alias{return_local_moves}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Return all possible local moves: growth, trim, recenter}
\description{
Given a current configuration theta, return all posssible growth, trim and
recenter moves
}
\usage{return_local_moves(theta, overlap, cluster_coords)}
\arguments{
  \item{theta}{Current configuration}
  \item{overlap}{See \code{\link{create_geo_objects}}}
  \item{cluster_coords}{See \code{\link{create_geo_objects}}}
}
\value{
List with 3 NumericMatrix's, corresponding to each of the possible moves
}
% \references{}
\author{Albert Y. Kim}
\seealso{\code{\link{return_death_moves}}, \code{\link{return_birth_moves}}}
%\examples{}
\keyword{internal}
