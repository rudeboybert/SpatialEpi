\name{GammaPriorCh}
\alias{GammaPriorCh}
\title{Compute Parameters to Calibrate a Gamma Distribution}
\description{
	Compute parameters to calibrate the prior distribution of a relative risk that has a gamma distribution.}
\usage{
GammaPriorCh(theta, prob, d)
}
\arguments{
  \item{theta}{upper quantile}
  \item{prob}{upper quantile}
  \item{d}{degrees of freedom}
}
%\details{}
\value{
A list containing
  \item{a}{shape parameter}
  \item{b}{rate parameter}
}
%\references{}
\author{Jon Wakefield}
%\note{}
\seealso{\code{\link{LogNormalPriorCh}}}
\examples{
param <- GammaPriorCh(5, 0.975,1)
curve(dgamma(x,shape=param$a,rate=param$b),from=0,to=6,n=1000,ylab="density")
}
\keyword{ file }
