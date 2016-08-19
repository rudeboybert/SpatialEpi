\name{bayes_cluster}
\alias{bayes_cluster}
\title{Bayesian Cluster Detection Method}
\description{
Implementation of the Bayesian Cluster detection model of Wakefield and Kim 
(2013) for a study region with \code{n} areas. The prior and posterior 
probabilities of each of the \code{n.zones} single zones being a 
cluster/anti-cluster are estimated using Markov chain Monte Carlo. Furthermore,
the posterior probability of k clusters/anti-clusters is computed.  
}
\usage{
bayes_cluster(y, E, population, sp.obj, centroids, max.prop, shape, rate, J, pi0,
         n.sim.lambda, n.sim.prior, n.sim.post, burnin.prop = 0.1,
         theta.init = vector(mode="numeric", length=0))
}
\arguments{
  \item{y}{vector of length \code{n} of the observed number of disease in each area}
  \item{E}{vector of length \code{n} of the expected number of disease in each area}
  \item{population}{vector of length \code{n} of the population in each area}
  \item{sp.obj}{an object of class SpatialPolygons (See \link[sp]{SpatialPolygons-class}) representing the study region}
  \item{centroids}{\code{n x 2} table of the (x,y)-coordinates of the area centroids.  The coordinate system must be grid-based} 
  \item{max.prop}{maximum proportion of the study region's population each single zone can contain}
  \item{shape}{vector of length 2 of narrow/wide shape parameter for gamma prior on relative risk}
  \item{rate}{vector of length 2 of narrow/wide rate parameter for gamma prior on relative risk}
  \item{J}{maximum number of clusters/anti-clusters}
  \item{pi0}{prior probability of no clusters/anti-clusters}
  \item{n.sim.lambda}{number of importance sampling iterations to estimate lambda}
  \item{n.sim.prior}{number of MCMC iterations to estimate prior probabilities associated with each single zone}
  \item{n.sim.post}{number of MCMC iterations to estimate posterior probabilities associated with each single zone}
  \item{burnin.prop}{proportion of MCMC samples to use as burn-in}
  \item{theta.init}{Initial configuration used for MCMC sampling}
}
\value{
List containing
\item{prior.map}{A list containing, for each area: 1) \code{high.area} the prior probability of cluster membership, 2) \code{low.area} anti-cluster membership, and 3) \code{RR.est.area} smoothed prior estimates of relative risk}
\item{post.map}{A list containing, for each area: 1) \code{high.area} the posterior probability of cluster membership, 2) \code{low.area} anti-cluster membership, and 3) \code{RR.est.area} smoothed posterior estimates of the relative risk}
\item{pk.y}{posterior probability of k clusters/anti-clusters given y for k=0,...,J}
}
\references{Wakefield J. and Kim A.Y. (2013) A Bayesian model for cluster detection.
\emph{Biostatistics}, \bold{14}, 752--765.}
\author{Albert Y. Kim}
\seealso{\code{\link{kulldorff}}}
\examples{
## Note for the NYleukemia example, 4 census tracts were completely surrounded 
## by another unique census tract; when applying the Bayesian cluster detection 
## model in \code{\link{bayes_cluster}}, we merge them with the surrounding 
## census tracts yielding \code{n=277} areas.

## Load data and convert coordinate system from latitude/longitude to grid
data(NYleukemia)
sp.obj <- NYleukemia$spatial.polygon
population <- NYleukemia$data$population
cases <- NYleukemia$data$cases
centroids <- latlong2grid(NYleukemia$geo[, 2:3])

## Identify the 4 census tract to be merged into their surrounding census tracts 
remove <- NYleukemia$surrounded
add <- NYleukemia$surrounding

## Merge population and case counts and geographical objects accordingly
population[add] <- population[add] + population[remove]
population <- population[-remove]
cases[add] <- cases[add] + cases[remove]
cases <- cases[-remove]
sp.obj <-
  SpatialPolygons(sp.obj@polygons[-remove], proj4string=CRS("+proj=longlat +ellps=WGS84"))
centroids <- centroids[-remove, ]


## Set parameters
y <- cases
E <- expected(population, cases, 1)
max.prop <- 0.15
shape <- c(2976.3, 2.31)
rate <- c(2977.3, 1.31)
J <- 7
pi0 <- 0.95
n.sim.lambda <- 10^4
n.sim.prior <- 10^5
n.sim.post <- 10^5


## (Uncomment first) Compute output
#output <- bayes_cluster(y, E, population, sp.obj, centroids, max.prop, 
#  shape, rate, J, pi0, n.sim.lambda, n.sim.prior, n.sim.post)
#plotmap(output$prior.map$high.area, sp.obj)
#plotmap(output$post.map$high.area, sp.obj)
#plotmap(output$post.map$RR.est.area, sp.obj, log=TRUE)
#barplot(output$pk.y, names.arg=0:J, xlab="k", ylab="P(k|y)")                   
}
\keyword{file}
