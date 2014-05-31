\name{clean_moves_matrix}
\alias{clean_moves_matrix}
\title{Clean up proposed moves matrix}
\description{Cleans up proposed move matrix}
\usage{
clean_moves_matrix(theta, moves, n_zones)
}
\arguments{
  \item{theta}{Original configuration}
  \item{moves}{Proposed moves}
  \item{n_zones}{Total number of single zones in region}
}
\value{
A matrix with kstar rows, and each column being a proposed configuration, with
the component single zones being sorted.  
}
\author{Albert Y. Kim}
%\examples{}
\keyword{internal}
