\name{clean_moves_matrix}
\alias{clean_moves_matrix}
\title{Clean up proposed moves matrix}
\description{Cleans up proposed move matrix}
\usage{clean_moves_matrix(theta, moves, n_zones)}
\arguments{
  \item{theta}{original configuration}
  \item{moves}{proposed moves}
  \item{n_zones}{total number of single zones in region}
}
\value{
A matrix with \code{kstar} rows, and each column being a proposed configuration, with
the component single zone indices being sorted (i.e. each column).  
}
\author{Albert Y. Kim}
\keyword{internal}
