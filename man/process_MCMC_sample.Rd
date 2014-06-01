\name{process_MCMC_sample}
\alias{process_MCMC_sample}
\title{Process MCMC Sample}
\description{Takes the output of  and produce area-by-area summaries}
\usage{process_MCMC_sample(sample, param, RR.area, cluster.list, cutoffs)}
\arguments{
  \item{sample}{output of function}
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
%\seealso{\code{\link{MCMC}}}
\keyword{internal}
