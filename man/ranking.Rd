\name{ranking}
\alias{ranking}
\title{Get rankings of non-overlapping elements}
\description{Gives zones in order of non-overlap}
\usage{
ranking(variable, cluster.list)
}
\arguments{
  \item{variable}{ test statistic }
  \item{cluster.list}{list with areas in each zone}
}
%\details{}
\value{
  \item{ranks}{rank}
  \item{indices}{index}
}
\references{Kulldorff}
\author{Albert Y. Kim}
%\note{}
%\seealso{}
\examples{
data(pennLC)
}
% R documentation directory.
\keyword{internal}
