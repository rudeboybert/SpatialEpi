#' Bayesian Cluster Detection Method
#' @description Implementation of the Bayesian Cluster detection model of Wakefield and Kim (2013) for a study region with `n` areas. The prior and posterior probabilities of each of the `n.zones` single zones being a cluster/anti-cluster are estimated using Markov chain Monte Carlo. Furthermore, the posterior probability of k clusters/anti-clusters is computed.  
#'
#' @param y vector of length `n` of the observed number of disease in each area
#' @param E vector of length `n` of the expected number of disease in each area
#' @param population vector of length `n` of the population in each area
#' @param sp.obj an object of class SpatialPolygons
#' @param centroids `n x 2` table of the (x,y)-coordinates of the area centroids.  The coordinate system must be grid-based
#' @param max.prop maximum proportion of the study region's population each single zone can contain
#' @param shape vector of length 2 of narrow/wide shape parameter for gamma prior on relative risk
#' @param rate vector of length 2 of narrow/wide rate parameter for gamma prior on relative risk
#' @param J maximum number of clusters/anti-clusters
#' @param pi0 prior probability of no clusters/anti-clusters
#' @param n.sim.lambda number of importance sampling iterations to estimate lambda
#' @param n.sim.prior number of MCMC iterations to estimate prior probabilities associated with each single zone
#' @param n.sim.post number of MCMC iterations to estimate posterior probabilities associated with each single zone
#' @param burnin.prop proportion of MCMC samples to use as burn-in
#' @param theta.init Initial configuration used for MCMC sampling
#' 
#'
#' @return 
#' List containing
#' return(list(
#' prior.map=prior.map, 
#' post.map=post.map, 
#' pk.y=pk.y))
#' \item{prior.map}{A list containing, for each area: 1) `high.area` the prior probability of cluster membership, 2) `low.area` anti-cluster membership, and 3) `RR.est.area` smoothed prior estimates of relative risk}
#' \item{post.map}{A list containing, for each area: 1) `high.area` the posterior probability of cluster membership, 2) `low.area` anti-cluster membership, and 3) `RR.est.area` smoothed posterior estimates of the relative risk}
#' \item{pk.y}{posterior probability of k clusters/anti-clusters given y for k=0,...,J}
#'
#'
#'
#'
#' @references  Wakefield J. and Kim A.Y. (2013) A Bayesian model for cluster detection.
#' @author Albert Y. Kim  
#' 
#' @export
#' 
#' @examples 
#' ## Note for the NYleukemia example, 4 census tracts were completely surrounded 
#' ## by another unique census tract; when applying the Bayesian cluster detection 
#' ## model in [bayes_cluster()], we merge them with the surrounding 
#' ## census tracts yielding `n=277` areas.
#' 
#' ## Load data and convert coordinate system from latitude/longitude to grid
#' data(NYleukemia)
#' sp.obj <- NYleukemia$spatial.polygon
#' population <- NYleukemia$data$population
#' cases <- NYleukemia$data$cases
#' centroids <- latlong2grid(NYleukemia$geo[, 2:3])
#' 
#' ## Identify the 4 census tract to be merged into their surrounding census tracts 
#' remove <- NYleukemia$surrounded
#' add <- NYleukemia$surrounding
#' 
#' ## Merge population and case counts and geographical objects accordingly
#' population[add] <- population[add] + population[remove]
#' population <- population[-remove]
#' cases[add] <- cases[add] + cases[remove]
#' cases <- cases[-remove]
#' sp.obj <-
#'   SpatialPolygons(sp.obj@polygons[-remove], proj4string=CRS("+proj=longlat +ellps=WGS84"))
#' centroids <- centroids[-remove, ]
#' 
#' ## Set parameters
#' y <- cases
#' E <- expected(population, cases, 1)
#' max.prop <- 0.15
#' shape <- c(2976.3, 2.31)
#' rate <- c(2977.3, 1.31)
#' J <- 7
#' pi0 <- 0.95
#' n.sim.lambda <- 10^4
#' n.sim.prior <- 10^5
#' n.sim.post <- 10^5
#' 
#' ## (Uncomment first) Compute output
#' #output <- bayes_cluster(y, E, population, sp.obj, centroids, max.prop, 
#' #  shape, rate, J, pi0, n.sim.lambda, n.sim.prior, n.sim.post)
#' #plotmap(output$prior.map$high.area, sp.obj)
#' #plotmap(output$post.map$high.area, sp.obj)
#' #plotmap(output$post.map$RR.est.area, sp.obj, log=TRUE)
#' #barplot(output$pk.y, names.arg=0:J, xlab="k", ylab="P(k|y)")    
bayes_cluster <-
  function(y, E, population, sp.obj, centroids, max.prop, shape, rate, J, pi0,
           n.sim.lambda, n.sim.prior, n.sim.post, burnin.prop = 0.1,
           theta.init = vector(mode="numeric", length=0)
  ){
    print(paste("Algorithm started on:", date()))
    # MCMC parameters
    pattern <- c(0, 1)
    p.moves <- normalize(c(1, 1, 1, 1, 1))
    
    # Geographic info and create geographical objects to use
    geo.objects <- create_geo_objects(max.prop, population, centroids, sp.obj)
    overlap <- geo.objects$overlap
    cluster.coords <- geo.objects$cluster.coords    
    n <- length(E)
    n.zones <- nrow(cluster.coords)
    
    # Set cutoffs of the relative risk to determine low and high risk
    RR.range <- seq(0.7, 1/0.7, length=100000)
    pN <- dgamma(RR.range, shape[1], rate[1])
    pW <- dgamma(RR.range, shape[2], rate[2])
    cutoffs <- list(
      low=RR.range[min(which(pN/pW>1))],
      high=RR.range[max(which(pN/pW>1))]
    )
    
    
    #-------------------------------------------------------------------------------
    # Obtain values of lambda using importance sampling
    #-------------------------------------------------------------------------------
    # Set uniform prior on single zones.
    prior.z <- rep(1/n.zones, n.zones)
    log_prior.z <- log(prior.z) - log(sum(prior.z))
    
    # Estimate q and lambda
    results <- estimate_lambda(n.sim.lambda, J, prior.z, overlap, pi0)
    lambda <- results$lambda
    prior.j <- results$prior.j
    print(paste("Importance sampling of lambda complete on:", date()))
    
    
    #-------------------------------------------------------------------------------
    # Obtain prior map
    #-------------------------------------------------------------------------------
    # Generate MCMC samples from prior
    prior.chain <- MCMC_simulation(n.sim.prior, pattern, theta.init, overlap, 
                                   cluster.coords, p.moves, J, prior.z, lambda)
    
    # Trim burn-in
    burnin <- n.sim.prior * burnin.prop
    prior.sample <- prior.chain$sample[-c(1:burnin)]
    
    # Prior Probs of Cluster Membership for each Area
    param.prior.zone <- list(shape=rep(shape[2],n.zones), rate=rep(rate[2],n.zones))
    RR.prior.area <- rep(shape[1]/rate[1], n)
    prior.map <- process_MCMC_sample(prior.sample, param.prior.zone, RR.prior.area, 
                                     overlap$cluster.list, cutoffs)
    print(paste("Prior map MCMC complete on:", date()))
    
    
    #-------------------------------------------------------------------------------
    # Obtain posterior map
    #-------------------------------------------------------------------------------
    # Single Zone Values
    yz <- sapply(overlap$cluster.list, function(x){sum(y[x])})
    Ez <- sapply(overlap$cluster.list, function(x){sum(E[x])})
    log_BF.z <- coeff(y, E, shape, rate, overlap$cluster.list)
    BF.z <- exp(log_BF.z)
    log_post.z <- log_prior.z + log_BF.z
    # Do NOT normalize this quantity
    post.z <- exp(log_post.z)
    
    # Generate MCMC samples from prior
    post.chain <- MCMC_simulation(n.sim.post, pattern, theta.init, overlap, 
                                  cluster.coords, p.moves, J, post.z, lambda)
    
    # Trim burn-in
    burnin <- n.sim.post * burnin.prop
    post.sample <- post.chain$sample[-c(1:burnin)]
    
    # Prior Probs of Cluster Membership for each area
    param.post.zone <- list(shape=c(yz+shape[2]), rate=c(Ez+rate[2]))
    RR.post.area <- (y+shape[1])/(E+rate[1])
    post.map <- process_MCMC_sample(post.sample, param.post.zone, RR.post.area, 
                                    overlap$cluster.list, cutoffs)
    
    # Posterior probs of k cluster/anti-clusters for k=0,...,J
    k.vector <- table(sapply(post.chain$sample,length))
    k.names <- as.numeric(names(k.vector))
    k.vector <- normalize(k.vector)
    pk.y <- rep(0, J+1)
    pk.y[k.names+1] <- k.vector
    print(paste("Posterior estimation complete on:", date()))
    
    
    
    #-------------------------------------------------------------------------------
    # Output
    #-------------------------------------------------------------------------------
    return(list(
      prior.map=prior.map, 
      post.map=post.map, 
      pk.y=pk.y))
    
    }

