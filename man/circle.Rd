\name{circle}
\alias{circle}
\title{Compute cartesian coordinates of a cluster center and radius}
\description{
This function is used for plotting purposes
}
\usage{
circle(geo, cluster.center, cluster.end)
}
\arguments{
  \item{geo}{A \code{n x 2} table of the x-coordinate and y-coordinates of the centroids of each area}
  \item{cluster.center}{The area index (an integer between \code{1} and \code{n}) indicating the center of the circle}
  \item{cluster.end}{The area index (an integer between \code{1} and \code{n}) indicating the area at the end of the circle}
}
%\details{}
\value{
  \item{cluster.radius}{A data frame that you can plot}

}
\author{Albert Y. Kim}
%\note{}
\seealso{\code{\link{zones}} }
\examples{
data(pennLC)
geo <- pennLC$geo[,2:3]
plot(geo,type='n')
text(geo,labels=1:nrow(geo))
lines( circle(geo, 23, 46), col = "red" ) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
