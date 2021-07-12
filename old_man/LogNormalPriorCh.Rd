\name{LogNormalPriorCh}
\alias{LogNormalPriorCh}
\title{Compute Parameters to Calibrate a Log-normal Distribution}
\description{
	Compute parameters to calibrate the prior distribution of a relative risk that has a log-normal distribution.
}
\usage{
LogNormalPriorCh(theta1, theta2, prob1, prob2)
}
\arguments{
  \item{theta1}{lower quantile}
  \item{theta2}{upper quantile}
  \item{prob1}{lower probability}
  \item{prob2}{upper probability}
}
%\details{}
\value{
A list containing:
  \item{mu}{mean of log-normal distribution}
  \item{sigma}{variance of log-normal distribution}
}
%\references{}
\author{Jon Wakefield}
\seealso{\code{\link{GammaPriorCh}}}
\examples{
# Calibrate the log-normal distribution s.t. the 95% confidence interval is [0.2, 5]
param <- LogNormalPriorCh(0.2, 5, 0.025, 0.975)
curve(dlnorm(x,param$mu,param$sigma), from=0, to=6, ylab="density")
}
\keyword{file}
