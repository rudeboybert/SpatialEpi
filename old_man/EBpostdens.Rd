\name{EBpostdens}
\alias{EBpostdens}
\title{Produce plots of emprical Bayes posterior densities when the data Y are Poisson with expected number E and relative risk theta, with the latter having a gamma distribution with known values alpha and beta, which are estimated using empirical Bayes.}
\description{
This function produces plots of empirical Bayes posterior densities which are gamma distributions 
with parameters (alpha+Y, (alpha+E*mu)/mu) where mu = exp(x beta). The SMRs are drawn on for comparison.
}
\usage{
EBpostdens(Y, E, alpha, beta, Xrow = NULL, lower = NULL, upper = NULL, main = "")
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
  \item{lower}{
%%     ~~Describe \code{lower} here~~
}
  \item{upper}{
%%     ~~Describe \code{upper} here~~
}
  \item{main}{
%%     ~~Describe \code{main} here~~
}
}
\value{
A plot containing the gamma posterior distribution 
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{
Jon Wakefield
}
\seealso{\code{\link{EBpostthresh}}, \code{\link{eBayes}}}
\examples{
data(scotland)
Y <- scotland$data$cases
E <- scotland$data$expected
ebresults <- eBayes(Y,E)
EBpostdens(Y[1], E[1], ebresults$alpha, ebresults$beta, lower=0, upper=15,
main="Area 1")
}
\keyword{file}
