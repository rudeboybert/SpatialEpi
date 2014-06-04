\name{poissonLogLkhd}
\alias{poissonLogLkhd}
\title{Compute Poisson Likelihoods}
\description{Compute poisson likelihood ratio test statistic for Kulldorff method}
\usage{poissonLogLkhd(cz, nz, N, C)}
\arguments{
  \item{cz}{count inside zone}
  \item{nz}{expected count inside zone}
  \item{N}{total expected count in region}
  \item{C}{total number of cases in region}
}
\value{Poisson likelihood ratio test statistic}
\references{ Kulldorff M. and Nagarwalla N. (1995) Spatial disease clusters: Detection and Inference.\emph{Statistics in Medicine}, \bold{14}, 799--810.}
\author{Albert Y. Kim}
\seealso{\code{\link{binomialLogLkhd}}, \code{\link{kulldorff}}}
\keyword{internal}
