\name{binomialLogLkhd}
\alias{binomialLogLkhd}
\title{Compute Binomial Likelihoods}
\description{Compute binomial likelihood ratio test statistic for Kulldorff method}
\usage{binomialLogLkhd(cz, nz, N, C)}
\arguments{
  \item{cz}{Count inside zone}
  \item{nz}{Expected count inside zone}
  \item{N}{Total expected count in region}
  \item{C}{Total number of cases in region}
}
\value{LRT statistic}
\author{Albert Y. Kim}
\seealso{\code{\link{poissonLogLkhd}}}
\keyword{internal}
