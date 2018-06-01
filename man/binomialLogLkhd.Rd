\name{binomialLogLkhd}
\alias{binomialLogLkhd}
\title{Compute Binomial Likelihoods}
\description{Compute binomial likelihood ratio test statistic for Kulldorff method}
\usage{binomialLogLkhd(cz, nz, N, C)}
\arguments{
  \item{cz}{count inside zone}
  \item{nz}{expected count inside zone}
  \item{N}{total expected count in region}
  \item{C}{total number of cases in region}
}
\value{Binomial likelihood ratio test statistic}
\references{ Kulldorff M. and Nagarwalla N. (1995) Spatial disease clusters: Detection and Inference.\emph{Statistics in Medicine}, \bold{14}, 799--810.}
\author{Albert Y. Kim}
\seealso{\code{\link{poissonLogLkhd}}, \code{\link{kulldorff}}}
\keyword{internal}
