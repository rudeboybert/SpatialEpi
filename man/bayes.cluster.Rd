\name{bayes.cluster}
\alias{bayes.cluster}
\title{Bayesian Cluster Detection}
\description{
Implementation of the Bayesian Cluster detection model of Wakefield and Kim for a study region with \code{n} areas. The prior and posterior probabilities of each of the \code{n.zones} single zones being a cluster/anti-cluster are estimated using Markov chain Monte Carlo.
}
\usage{
bayes.cluster(E, cases, population, centroids, map, max.prop, k, shape, rate, J, 
pi0, n.sim.imp, n.sim.prior, n.sim.post)
}
\arguments{
  \item{E}{vector of length \code{n} of the expected number of disease in each area}
  \item{cases}{vector of length \code{n} of the observed number of disease in each area}
  \item{population}{vector of length \code{n} of the population in each area}
  \item{centroids}{an \code{n x 2} table of the (x,y)-coordinates of the area centroids.  The coordinate system must be grid-based} 
  \item{map}{an object of class SpatialPolygons (See \link[sp]{SpatialPolygons-class}) representing the study region}
  \item{max.prop}{maximum proportion of study region's population each single zone can contain}
  \item{k}{parameter to tune the prior probability of each single zone being a cluster/anti-cluster as a function of geographical area size}
  \item{shape}{narrow/wide shape parameter for gamma prior on relative risk}
  \item{rate}{narrow/wide rate parameter for gamma prior on relative risk}
  \item{J}{maximum number of clusters/anti-clusters considered}
  \item{pi0}{prior probability of no clusters/anti-clusters i.e. probability that J=0}
  \item{n.sim.imp}{number of importance sampling iterations to estimate lambda}
  \item{n.sim.prior}{number of MCMC iterations to estimate prior probabilities associated with each single zone}
  \item{n.sim.post}{number of MCMC iterations to estimate posterior probabilities associated with each single zone}
}
\value{
List containing
\item{prior.map}{A sublist containing, for each area: 1) \code{high.area} the prior probability of cluster membership, 2) \code{low.area} anti-cluster membership, and 3) \code{RR.est.area} smoothed prior estimate of relative risk}
\item{post.map}{A sublist containing, for each area: 1) \code{high.area} the posterior probability of cluster membership, 2) \code{low.area} anti-cluster membership, and 3) \code{RR.est.area} smoothed posterior estimates of the relative risk}
\item{pj.y}{posterior probability of j clusters/anti-clusters given y}
}
\references{Wakefield, J. and Kim, A.y. (2012) 
A Bayesian Model for Cluster Detection}
\author{Albert Y. Kim}
\seealso{\code{\link{kulldorff}}}
\examples{
## Note for the NYleukemia example, 4 census tracts were completely surrounded by another 
## unique census tract; when applying the Bayesian cluster detection model in 
## \code{\link{bayes.cluster}}, we merge them with the surrounding census tracts yielding 
## \code{n=277} areas.

## Load data and convert coordinate system from latitude/longitude to grid
data(NYleukemia)
map <- NYleukemia$spatial.polygon
population <- NYleukemia$data$population
cases <- NYleukemia$data$cases
centroids <- latlong2grid(NYleukemia$geo[, 2:3])

## Identify the 4 census tract to be merged into their surrounding census tracts.  
remove <- NYleukemia$surrounded
add <- NYleukemia$surrounding

## Merge population and case counts and geographical objects accordingly
population[add] <- population[add] + population[remove]
population <- population[-remove]
cases[add] <- cases[add] + cases[remove]
cases <- cases[-remove]
map <- SpatialPolygons(map@polygons[-remove], proj4string=CRS("+proj=longlat"))
centroids <- centroids[-remove, ]

## Expected numbers
E <- expected(population, cases, 1)

## Set Parameters
max.prop <- 0.15
k <- 0.00005
shape <- c(2976.3, 2.31); rate <- c(2977.3, 1.31)
J <- 7
pi0 <- 0.95
n.sim.imp <- 0.5*10^4
n.sim.prior <- 0.5*10^4
n.sim.post <- 0.5*10^5

## (Uncomment first) Compute output
# output <- bayes.cluster(E, cases, population, centroids, map, max.prop, k, 
# 	shape, rate, J, pi0, n.sim.imp, n.sim.prior, n.sim.post)
# plotmap(output$prior.map$high.area, map)
# plotmap(output$post.map$high.area, map)
# plotmap(output$post.map$RR.est.area, map, log=TRUE)
# barplot(output$pj.y, names.arg=0:J, xlab="j", ylab="P(j|y)")                   
}
\keyword{file}
