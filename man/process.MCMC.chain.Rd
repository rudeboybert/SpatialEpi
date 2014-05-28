\name{process.MCMC.chain}
\alias{process.MCMC.chain}
\title{Process MCMC Chain}
\description{Takes the output of \code{\link{MCMC}} and produce area-by-area summaries}
\usage{process.MCMC.chain(chain, param, RR.area, cluster.list, cutoffs)}
\arguments{
  \item{chain}{output of \code{\link{MCMC}} function}
  \item{param}{list of wide/narrow shape/rate gamma prior parameters}
  \item{RR.area}{Relative risk associated with each of the \code{n} area}
  \item{cluster.list}{geo objects}
  \item{cutoffs}{cutoffs used to declare highs (clusters) and lows (anti-clusters)}
}
\value{
\item{high.area}{Probability of cluster membership for each area}
\item{low.area}{Probability of anti-cluster membership for each area}
\item{RR.est.area}{Smoothed relative risk estimates for each area}
}
\seealso{\code{\link{MCMC}}}
\keyword{internal}
