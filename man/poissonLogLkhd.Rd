\name{poissonLogLkhd}
\alias{poissonLogLkhd}
\title{Compute Poisson Likelihoods}
\description{Compute poisson likelihood ratio test statistic for Kulldorff method}
\usage{poissonLogLkhd(cz, nz, N, C)}
\arguments{
  \item{cz}{Count inside zone}
  \item{nz}{Expected count inside zone}
  \item{N}{Total expected count in region}
  \item{C}{Total number of cases in region}
}
\value{LRT statistic}
\author{Albert Y. Kim}
\seealso{\code{\link{binomialLogLkhd}}
}
%\examples{}
\keyword{internal}
