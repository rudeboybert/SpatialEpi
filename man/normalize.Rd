\name{normalize}
\alias{normalize}
\title{Normalize a vector}
\description{Function to Normalize a vector, possibly over a restricted set of indices}
\usage{normalize(p, include.indices = NA)}
\arguments{
  \item{p}{vector to normalize}
  \item{include.indices}{indices to normalize over. Returns 0 for indices not included}
}
\value{
\item{p.prime}{normalized vector}
}
\examples{
x <- c(1:5)
normalize(x, c(1,1,1,0,0))
}
\keyword{file}
