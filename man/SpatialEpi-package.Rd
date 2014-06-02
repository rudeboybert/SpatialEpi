\name{SpatialEpi-package}
\alias{SpatialEpi}
\title{Methods and Data for Spatial Epidemiology}
\description{This package contains methods for cluster detection and disease 
mapping, as well as plotting methods using the \link[sp]{sp} package.}
\author{
	\tabular{ll}{
		Cici Chen \tab \email{cicichen@u.washington.edu}\cr
		Albert Y. Kim \tab \email{albert@stat.washington.edu}\cr
		Michelle E. Ross \tab \email{micher3@u.washington.edu}\cr
		Jon Wakefield \tab \email{jonno@stat.washington.edu}
	}
}
\details{
\tabular{ll}{
Package: \tab SpatialEpi\cr
Type: \tab Package\cr
Version: \tab 1.2.1\cr
Date: \tab 2014-06-03\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}
}
\keyword{ package }
\seealso{\code{\link[sp]{sp}}}
\examples{
data(pennLC)
data(scotland)
data(NYleukemia)
?kulldorff
?besag_newell
?bayes_cluster
}
