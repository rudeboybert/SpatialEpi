\name{NumericVectorEquality}
\alias{NumericVectorEquality}
\title{Test if two numeric vectors are equal}
\description{Test if two numeric vectors are equal in their length and their 
elements}
\usage{
NumericVectorEquality(x, y)
}
\arguments{
  \item{x}{NumericVector}
  \item{y}{NumericVector}
}
\value{1 if equal, 0 if not.  }
\author{Albert Y. Kim}
\examples{
x <- c(1:10)
y <- rev(c(10:1))
NumericVectorEquality(x,y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{file}
