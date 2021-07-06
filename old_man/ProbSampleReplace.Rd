\name{ProbSampleReplace}
\alias{ProbSampleReplace}
\title{
Draw a sample of size 1 given probabilities
}
\description{
This uses the RcppArmadillo implementation of the sample() command to draw a
sample of size 1 from 0:(n-1) given probability vector p of length n. This is 
used for sampling interal to an Rcpp method, and its start index is 0 to match
vector indexing in C++.
}
\usage{
ProbSampleReplace(prob)
}
\arguments{
  \item{prob}{Probability vector}
}
\value{Value between 1 and n}
\author{Albert Y. Kim}
\examples{
p <- rep(0.2, 5)
ProbSampleReplace(p)
}
\keyword{internal}
