\name{return_local_moves}
\alias{return_local_moves}
\title{Return all possible local moves: growth, trim, recenter}
\description{Given a current configuration theta of size \code{k} single zones, return all configurations where we either grow, trim, or recenter one of the \code{k} single zones.}
\usage{return_local_moves(theta, overlap, cluster_coords)}
\arguments{
  \item{theta}{current configuration}
  \item{overlap}{see \code{\link{create_geo_objects}}}
  \item{cluster_coords}{see \code{\link{create_geo_objects}}}
}
\value{list of 3 \code{NumericMatrix}, corresponding to each of the possible moves (\code{k} rows, each possible move being a column)}
\author{Albert Y. Kim}
\seealso{\code{\link{return_death_moves}}, \code{\link{return_birth_moves}}}
\keyword{internal}
