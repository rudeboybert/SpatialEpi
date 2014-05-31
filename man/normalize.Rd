\name{normalize}
\alias{normalize}
\title{
Normalize vector to sum to 1.
}
\description{
Divide each element in the vector by the sum of the vector elements.  
}
\usage{
normalize(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector to be normalized}
}
\details{}
\value{Normalized vector}
\author{Albert Y. Kim}
\examples{
x <- rep(1, 10)
x <- normalize(x)
x
}
\keyword{file}