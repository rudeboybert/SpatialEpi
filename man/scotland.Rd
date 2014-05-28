\name{scotland}
\alias{scotland}
\docType{data}
\title{Lip Cancer in Scotland}
\description{
County-level (n=56) data for lip cancer among males in Scotland between 1975-1980
}
\usage{data(scotland)}
\format{
List containing:
\tabular{ll}{
\code{geo} \tab a table of county IDs, x-coordinates (eastings) and y-coordinates (northings) of the geographic centroid of each county.\cr
\code{data} \tab a table of county IDs, number of cases, population and strata information \cr
\code{spatial.polygon} \tab a Spatial Polygons class (See \link[sp]{SpatialPolygons-class}) map of Scotland \cr
\code{polygon} \tab a polygon map of Scotland (See \code{\link{polygon2spatial.polygon}})\cr
}
}
%\details{}

\source{
Kemp I., Boyle P., Smans M. and Muir C. (1985) Atlas of cancer in Scotland, 1975-1980, incidence and epidemiologic perspective \emph{International Agency for Research on Cancer} \bold{72}.
}
\references{
Clayton D. and Kaldor J. (1987) Empirical Bayes estimates of age-standardized relative risks for use in disease mapping.  \emph{Biometrics}, \bold{43}, 671--681
}
\seealso{
\code{\link{mapvariable}}, \code{\link{polygon2spatial.polygon}}, \code{\link{pennLC}}, 
\code{\link{NYleukemia}}
}
\examples{
data(scotland)
data <- scotland$data
scotland.map <- scotland$spatial.polygon

SMR <- data$cases/data$expected
mapvariable(SMR,scotland.map)
}
\keyword{datasets}