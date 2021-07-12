\name{ldmultinom}
\alias{ldmultinom}
\title{log multinomial Density}
\description{Compute log of density of multinomial}
\usage{
ldmultinom(x, prob)
}
\arguments{
  \item{x}{vector of multinomial counts}
  \item{prob}{vector of bin probabilities}
}
\value{log of density of multinomial}
\author{Albert Y. Kim}
\examples{
p <- rep(0.2, 5)
x <- rmultinom(1, 10, p)
dmultinom(x, prob=p, log=TRUE)
ldmultinom(x, p)
}
\keyword{file}
