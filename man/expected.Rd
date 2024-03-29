% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expected.R
\name{expected}
\alias{expected}
\title{Compute Expected Numbers of Disease}
\usage{
expected(population, cases, n.strata)
}
\arguments{
\item{population}{a vector of population counts for each strata in each area}

\item{cases}{a vector of the corresponding number of cases}

\item{n.strata}{number of strata considered}
}
\value{
\item{expected.cases}{a vector of the expected numbers of disease for each area}
}
\description{
Compute the internally indirect standardized expected numbers of disease.
}
\details{
The \code{population} and \code{cases} vectors must be \emph{balanced}: all counts are sorted by area first, and then within each area the counts for all strata are listed (even if 0 count) in the same order.
}
\examples{
data(pennLC)
population <- pennLC$data$population
cases <- pennLC$data$cases
## In each county in Pennsylvania, there are 2 races, gender and 4 age bands 
## considered = 16 strata levels
pennLC$data[1:16,]
expected(population, cases, 16)


}
\references{
Elliot, P. et al. (2000) \emph{Spatial Epidemiology:  Methods and Applications}.  Oxford Medical Publications.
}
\author{
Albert Y. Kim
}
