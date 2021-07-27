
globalVariables(c(
  "rmultinom","hist","abline","legend"
))



#' Kulldorff Cluster Detection Method
#'
#' @description Kulldorff spatial cluster detection method for a study region with \code{n} areas.  The method constructs \emph{zones} by consecutively aggregating nearest-neighboring areas until a proportion of the total study population is included.  Given the observed number of cases, the likelihood of each zone is computed using either binomial or poisson likelihoods. The procedure reports the zone that is the \emph{most likely cluster} and generates significance measures via Monte Carlo sampling.  Further, \emph{secondary clusters}, whose Monte Carlo p-values are below the \eqn{\alpha}-threshold, are reported as well.
#' 
#' @param geo an \code{n x 2} table of the (x,y)-coordinates of the area centroids
#' @param cases aggregated case counts for all \code{n} areas
#' @param population aggregated population counts for all \code{n} areas
#' @param expected.cases expected numbers of disease for all \code{n} areas
#' @param pop.upper.bound the upper bound on the proportion of the total population each zone can include
#' @param n.simulations number of Monte Carlo samples used for significance measures
#' @param alpha.level alpha-level threshold used to declare significance
#' @param plot flag for whether to plot histogram of Monte Carlo samples of the log-likelihood of the most likely cluster
#' 
#' @details If \code{expected.cases} is specified to be \code{NULL}, then the binomial likelihood is used.  Otherwise, a Poisson model is assumed.  Typical values of \code{n.simulations} are \code{99}, \code{999}, \code{9999}
#' @return
#' List containing:
#' \item{most.likely.cluster}{information on the most likely cluster}
#' \item{secondary.clusters}{information on secondary clusters, if none \code{NULL} is returned}
#' \item{type}{type of likelihood}
#' \item{log.lkhd}{log-likelihood of each zone considered}
#' \item{simulated.log.lkhd}{\code{n.simulations} Monte Carlo samples of the log-likelihood of the most likely cluster}
#' 
#' 
#' @note The \code{most.likely.cluster} and \code{secondary.clusters} list elements are themselves lists reporting:\cr\cr
#' \tabular{ll}{
#'  \code{location.IDs.included} \tab ID's of areas in cluster, in order of distance\cr
#'  \code{population} \tab population of cluster\cr
#'  \code{number.of.cases} \tab number of cases in cluster\cr
#'  \code{expected.cases} \tab expected number of cases in cluster\cr
#'  \code{SMR} \tab estimated SMR of cluster\cr
#'  \code{log.likelihood.ratio} \tab log-likelihood of cluster\cr
#'  \code{monte.carlo.rank} \tab rank of lkhd of cluster within Monte Carlo simulated values\cr
#'  \code{p.value} \tab Monte Carlo \eqn{p}-value\cr
#' }
#' 
#' @references SatScan:  Software for the spatial, temporal, and space-time scan statistics \url{https://www.satscan.org/} Kulldorff, M. (1997) A spatial scan statistic. \emph{Communications in Statistics: Theory and Methods}, \bold{26}, 1481--1496.
#' Kulldorff M. and Nagarwalla N. (1995) Spatial disease clusters: Detection and Inference.
#' \emph{Statistics in Medicine}, \bold{14}, 799--810.
#' @author Albert Y. Kim
#' @export
#'
#' @examples
#' ## Load Pennsylvania Lung Cancer Data
#' data(pennLC)
#' data <- pennLC$data
#'
#' ## Process geographical information and convert to grid
#' geo <- pennLC$geo[,2:3]
#' geo <- latlong2grid(geo)
#'
#' ## Get aggregated counts of population and cases for each county
#' population <- tapply(data$population,data$county,sum)
#' cases <- tapply(data$cases,data$county,sum)
#'
#' ## Based on the 16 strata levels, computed expected numbers of disease
#' n.strata <- 16
#' expected.cases <- expected(data$population, data$cases, n.strata)
#'
#' ## Set Parameters
#' pop.upper.bound <- 0.5
#' n.simulations <- 999
#' alpha.level <- 0.05
#' plot <- TRUE
#'
#' ## Kulldorff using Binomial likelihoods
#' binomial <- kulldorff(geo, cases, population, NULL, pop.upper.bound, n.simulations, 
#'                      alpha.level, plot)
#' cluster <- binomial$most.likely.cluster$location.IDs.included
#'
#' ## plot
#' plot(pennLC$spatial.polygon,axes=TRUE)
#' plot(pennLC$spatial.polygon[cluster],add=TRUE,col="red")
#' title("Most Likely Cluster")
#'
#' ## Kulldorff using Poisson likelihoods
#' poisson <- kulldorff(geo, cases, population, expected.cases, pop.upper.bound, 
#'                     n.simulations, alpha.level, plot)
#' cluster <- poisson$most.likely.cluster$location.IDs.included
#'
#' ## plot
#' plot(pennLC$spatial.polygon,axes=TRUE)
#' plot(pennLC$spatial.polygon[cluster],add=TRUE,col="red")
#' title("Most Likely Cluster Controlling for Strata")
kulldorff <-
function(geo, cases, population, expected.cases=NULL, pop.upper.bound, 
         n.simulations, alpha.level, plot=TRUE){

#-------------------------------------------------------------------------------
# Initialization 
#-------------------------------------------------------------------------------
# Determine likelihood type: binomial or poisson
if(is.null(expected.cases)){
	type <- "binomial"
	denominator <- population
	expected.cases <- sum(cases) * (denominator/sum(denominator))
# poisson case	
}else{
	type <- "poisson"
	denominator <- expected.cases
}

# Get geographic information
geo.results <- zones(geo, population, pop.upper.bound)
nearest.neighbors <- geo.results$nearest.neighbors
cluster.coords <- geo.results$cluster.coords
n.zones <- nrow(cluster.coords)


#-------------------------------------------------------------------------------
# Observed statistic computation
#-------------------------------------------------------------------------------
lkhd <- computeAllLogLkhd(cases, denominator, nearest.neighbors, n.zones, type)

# Get areas included in most likely cluster
cluster.index <- which.max(lkhd)

# cluster center and radial area
center <- cluster.coords[cluster.index,1]
end <- cluster.coords[cluster.index,2]

# list of all areas included in cluster	
cluster <- nearest.neighbors[[center]]
cluster <- cluster[1:which(cluster == end)]


#-------------------------------------------------------------------------------
# Compute Monte Carlo randomized p-value
#-------------------------------------------------------------------------------
# Simulate cases under null hypothesis of no area effects i.e. conditioned on E
perm <- rmultinom(n.simulations, round(sum(cases)), prob=denominator)

# Compute simulated lambda's:  max log-lkhd in region
sim.lambda <- kulldorffMC(perm, denominator, nearest.neighbors, n.zones, type)

# Compute Monte Carlo p-value
combined.lambda <- c(sim.lambda, max(lkhd))
p.value <- 1-mean(combined.lambda < max(lkhd))

# Plot histogram
if(plot){
	hist(combined.lambda, 
       main="Monte Carlo Distribution of Lambda",
       xlab=expression(log(lambda)))
	abline(v=max(lkhd), col="red")
	legend("top",
		c(paste("Obs. log(Lambda) = ",round(max(lkhd),3),sep=""), 
			paste("p-value = ", round(p.value,log10(n.simulations + 1)),sep="")),
		lty=c(1, 1), 
		col=c("red","white"), 
		bty="n"
	)
}


#-------------------------------------------------------------------------------
# Create Most Likely Cluster Object
#-------------------------------------------------------------------------------
most.likely.cluster = list(
	location.IDs.included = cluster,
	population = sum(population[cluster]),	
	number.of.cases = sum(cases[cluster]),
	expected.cases = sum(expected.cases[cluster]),
	SMR = sum(cases[cluster])/sum(expected.cases[cluster]),
	log.likelihood.ratio = lkhd[cluster.index],
	monte.carlo.rank = sum(combined.lambda >= lkhd[cluster.index]),
	p.value = p.value
)


#-------------------------------------------------------------------------------
# Investigate Secondary Clusters
#-------------------------------------------------------------------------------
# running list of areas already covered to test overlap
current.cluster <- cluster
secondary.clusters <- NULL

# go through log-likelihoods in decreasing order, skipping largest which 
# corresponds to most likely cluster
indices <- order(lkhd, decreasing=TRUE)

for(i in 2:length(indices)){
	#---------------------------------------------------
	# Get areas included in new potential cluster
	new.cluster.index <- indices[i]
	new.center <- cluster.coords[new.cluster.index, 1]
	new.end <- cluster.coords[new.cluster.index, 2]
	
	new.cluster <- nearest.neighbors[[new.center]]
	new.cluster <- new.cluster[1:which(new.cluster == new.end)]	
	#---------------------------------------------------
	# if there is no overlap between existing clusters and the new potential cluster
	if(length(intersect(new.cluster,current.cluster)) == 0){		
		new.secondary.cluster <- list(
			location.IDs.included = new.cluster,
			population = sum(population[new.cluster]),
			number.of.cases = sum(cases[new.cluster]),
			expected.cases = sum(expected.cases[new.cluster]),
			SMR = sum(cases[new.cluster])/sum(expected.cases[new.cluster]),
			log.likelihood.ratio = lkhd[new.cluster.index],
			monte.carlo.rank = sum(combined.lambda >= lkhd[new.cluster.index]),
			p.value =  1 - mean(combined.lambda < lkhd[new.cluster.index])
		)

		# If no longer significant, end loop.
		if(new.secondary.cluster$p.value > alpha.level){
			break
		}
		# Otherwise add to existing clusters
		secondary.clusters <- c(secondary.clusters,list(new.secondary.cluster))
		current.cluster <- unique(c(current.cluster, new.cluster))
	}		
}


#-------------------------------------------------------------------------------
# Output results
#-------------------------------------------------------------------------------
results <- list(
	most.likely.cluster=most.likely.cluster,
	secondary.clusters=secondary.clusters,
	type = type,
	log.lkhd=lkhd,
	simulated.log.lkhd=sim.lambda
)
return(results)
}
