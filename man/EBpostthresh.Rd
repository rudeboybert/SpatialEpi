\name{EBpostthresh}
\alias{EBpostthresh}
\title{Produce the probabilities of exceeding a threshold given a posterior gamma distribution.}
\description{
This function produces the posterior probabilities of exceeding a threshold given a gamma distributions with parameters (alpha+Y, (alpha+E*mu)/mu) where mu = exp(x beta). This model arises from Y being Poisson with mean theta times E where theta is the relative risk and E are the expected numbers. The prior on theta is gamma with parameters alpha and beta. The parameters alpha and beta may be estimated using empirical Bayes. 
}
\usage{
EBpostthresh(Y, E, alpha, beta, Xrow = NULL, rrthresh)
}
\arguments{
  \item{Y}{
  observed disease counts
}
  \item{E}{
  expected disease counts
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{beta}{
%%     ~~Describe \code{beta} here~~
}
  \item{Xrow}{
%%     ~~Describe \code{Xrow} here~~
}
  \item{rrthresh}{
%%     ~~Describe \code{rrthresh} here~~
}
}
\value{
	Posterior probabilities of exceedence are returned.
}
\author{Jon Wakefield}
\seealso{\code{\link{eBayes}}}
\examples{
data(scotland)
Y <- scotland$data$cases
E <- scotland$data$expected
ebresults <- eBayes(Y,E)
# Find probabilities of exceedence of 3
thresh3 <- EBpostthresh(Y, E, alpha=ebresults$alpha, beta=ebresults$beta,
rrthresh=3)
mapvariable(thresh3, scotland$spatial.polygon)
}
\keyword{file}
