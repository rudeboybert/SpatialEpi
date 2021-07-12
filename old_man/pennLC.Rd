\name{pennLC}
\alias{pennLC}
\docType{data}
\title{Pennsylvania Lung Cancer}
\description{
County-level (n=67) population/case data for lung cancer in Pennsylvania in 2002, stratified on race (white vs non-white), gender and age (Under 40, 40-59, 60-69 and 70+).  Additionally, county-specific smoking rates.  
}
\usage{data(pennLC)}
\format{
List of 3 items:
\tabular{ll}{
\code{geo} \tab a table of county IDs, longitude/latitude of the geographic centroid of each county \cr
\code{data} \tab a table of county IDs, number of cases, population and strata information \cr
\code{smoking} \tab a table of county IDs and proportion of smokers\cr
\code{spatial.polygon} \tab an object of class SpatialPolygons (See \link[sp]{SpatialPolygons-class}) \cr
}
}
%\details{}
\seealso{ 
	\code{\link{scotland}}, \code{\link{NYleukemia}}
}
\source{
Population data was obtained from the 2000 decennial census, lung cancer 
and smoking data were obtained from the Pennsylvania Department of Health website:  \url{http://www.dsf.health.state.pa.us/}
}
%\references{}
\examples{
data(pennLC)
pennLC$geo
pennLC$data
pennLC$smoking

# Map smoking rates in Pennsylvania
mapvariable(pennLC$smoking[,2], pennLC$spatial.polygon)
}
\keyword{datasets}
