\name{return_death_moves}
\alias{return_death_moves}
\title{Return all death moves}
\description{Given a proposed configuration return all possible death moves}
\usage{return_death_moves(theta)}
\arguments{
  \item{theta}{vector of length \code{k} single zones corresponding to configuration}
}
\value{\code{NumericMatrix} with \code{k-1} rows (i.e. one less single zone) and each column
being a new possible config}
\author{Albert Y. Kim}
\seealso{\code{\link{return_birth_moves}},\code{\link{return_local_moves}}}
%\examples{}
\keyword{internal}
