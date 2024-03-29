% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pennLC.R
\docType{data}
\name{pennLC}
\alias{pennLC}
\title{Pennsylvania Lung Cancer}
\format{
List of 3 items
\describe{
\item{geo}{a table of county IDs, longitude/latitude of the geographic centroid of each county}
\item{data}{a table of county IDs, number of cases, population and strata information}
\item{smoking}{a table of county IDs and proportion of smokers}
\item{spatial.polygon}{an object of class SpatialPolygons}

}
}
\source{
Population data was obtained from the 2000 decennial census, lung cancer and smoking data were obtained from the Pennsylvania Department of Health website:  \url{https://www.health.pa.gov/Pages/default.aspx}
}
\usage{
pennLC
}
\description{
County-level (n=67) population/case data for lung cancer in Pennsylvania in 2002, stratified on race (white vs non-white), gender and age (Under 40, 40-59, 60-69 and 70+).  Additionally, county-specific smoking rates.
}
\examples{
data(pennLC)
pennLC$geo
pennLC$data
pennLC$smoking
# Map smoking rates in Pennsylvania
mapvariable(pennLC$smoking[,2], pennLC$spatial.polygon)




}
\keyword{datasets}
