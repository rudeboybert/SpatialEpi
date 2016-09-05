#' Bayesian Cluster Detection Method
#' 
#' Implementation of the Bayesian Cluster detection model of Wakefield and Kim 
#' (2013) for a study region with \code{n} areas. The prior and posterior 
#' probabilities of each of the \code{n.zones} single zones being a 
#' cluster/anti-cluster are estimated using Markov chain Monte Carlo. Furthermore,
#' the posterior probability of k clusters/anti-clusters is computed.  
#'
#' @param y vector of length \code{n} of the observed number of disease in each area
#' @param E vector of length \code{n} of the expected number of disease in each area
#' @param population vector of length \code{n} of the population in each area
#' @param sp.obj an object of class SpatialPolygons (See 
#' \link[sp]{SpatialPolygons-class}) representing the study region
#' @param centroids \code{n x 2} table of the (x,y)-coordinates of the area 
#' centroids.  The coordinate system must be grid-based
#' @param max.prop maximum proportion of the study region's population each 
#' single zone can contain
#' @param shape vector of length 2 of narrow/wide shape parameter for gamma 
#' prior on relative risk
#' @param rate vector of length 2 of narrow/wide rate parameter for gamma prior 
#' on relative risk
#' @param J maximum number of clusters/anti-clusters
#' @param pi0 prior probability of no clusters/anti-clusters
#' @param n_sim.lambda number of importance sampling iterations to estimate 
#' lambda
#' @param n_sim.prior number of MCMC iterations to estimate prior probabilities 
#' associated with each single zone
#' @param n_sim_post number of MCMC iterations to estimate posterior 
#' probabilities associated with each single zon
#' @param burnin.prop proportion of MCMC samples to use as burn-in. Defaults to 
#' 10 percent
#' @param theta.init Initial configuration used for MCMC sampling. Defaults to 
#' \code{NULL}
#' 
#' @importFrom magrittr %>%
#'
#' @return List containing
#' \item{prior.map}{A list containing, for each area: 1) \code{high.area} the 
#' prior probability of cluster membership, 2) \code{low.area} anti-cluster 
#' membership, and 3) \code{RR.est.area} smoothed prior estimates of relative risk}
#' \item{post.map}{A list containing, for each area: 1) \code{high.area} the 
#' posterior probability of cluster membership, 2) \code{low.area} anti-cluster 
#' membership, and 3) \code{RR.est.area} smoothed posterior estimates of the 
#' relative risk}
#' \item{pk.y}{posterior probability of k clusters/anti-clusters given y for 
#' k=0,...,J}
#' @export
#' @references Wakefield J. and Kim A.Y. (2013) A Bayesian model for cluster 
#' detection. \emph{Biostatistics}, \bold{14}, 752--765.
#' @seealso \code{\link{kuldorff}}
#'
#' @examples
#' # Load data
#' data(NYleukemia)
#' sp.obj <- NYleukemia$spatial.polygon
#' centroids <- latlong2grid(NYleukemia$geo[, 2:3])
#' population <- NYleukemia$data$population
#' cases <- NYleukemia$data$cases
#' 
#' # Set parameters
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
#' # Run method (uncomment first)
#' # output <- bayes_cluster(y, E, population, sp.obj, centroids, max.prop,
#' #   shape, rate, J, pi0, n.sim.lambda, n.sim.prior,
#' #   n.sim.post)
#' # plotmap(output$post_map$high_area, sp.obj)
bayes_cluster <- function(y, E, population, sp.obj, centroids, max.prop, shape, 
                          rate, J, pi0,
                          n.sim.lambda, n.sim.prior, n.sim_post, 
                          burnin.prop = 0.1,
                          theta.init = vector(mode="numeric", length=0)){
  print(paste("Algorithm started on:", date()))
  
  # MCMC parameters
  pattern <- c(0, 1)
  p_moves <- c(1, 1, 1, 1, 1) %>% normalize()
  
  # Geographic info and create geographical objects to use
  geo_objects <- create_geo_objects(centroids, NYleukemia$data$population, 0.15)
  cluster_coords <- geo_objects$cluster_coords
  overlap <- geo_objects$overlap
  n <- nrow(centroids)
  n_zones <- nrow(cluster_coords)
  
  # Set cutoffs of the relative risk to determine low and high risk
  RR_range <- seq(0.2, 1/0.2, length=100000)
  pN <- dgamma(RR_range, shape[1], rate[1])
  pW <- dgamma(RR_range, shape[2], rate[2])
  cutoffs <- list(
    low=RR_range[min(which(pN/pW>1))], high=RR_range[max(which(pN/pW>1))]
  )
  
  
  #-------------------------------------------------------------------------------
  # Obtain values of lambda using importance sampling
  #-------------------------------------------------------------------------------
  # Set uniform prior on single zones.
  prior_z <- rep(1/n_zones, n_zones)
  log_prior_z <- log(prior_z) - log(sum(prior_z))
  
  # Estimate q and lambda
  results <- estimate_lambda(n.sim.lambda, J, prior_z, overlap, pi0)
  lambda <- results$lambda
  prior_j <- results$prior_j
  print(paste("Importance sampling of lambda complete on:", date()))
  
  
  #-------------------------------------------------------------------------------
  # Obtain prior map
  #-------------------------------------------------------------------------------
  # Generate MCMC samples from prior
  prior_chain <- MCMC_simulation(n.sim.prior, pattern, theta.init, overlap, 
                                 as.matrix(cluster_coords), p_moves, J, prior_z, 
                                 lambda)
  
  # Trim burn-in
  burnin <- n.sim.prior * burnin.prop
  prior_sample <- prior_chain$sample[-c(1:burnin)]
  
  # Prior Probs of Cluster Membership for each Area
  param_prior_zone <- list(shape=rep(shape[2],n_zones), rate=rep(rate[2],n_zones))
  RR_prior_area <- rep(shape[1]/rate[1], n)
  prior_map <- process_MCMC_sample(prior_sample, param_prior_zone, RR_prior_area, 
                                   overlap$cluster_list, cutoffs)
  print(paste("Prior map MCMC complete on:", date()))
  
  
  #-------------------------------------------------------------------------------
  # Obtain posterior map
  #-------------------------------------------------------------------------------
  # Single Zone Values
  yz <- overlap$cluster_list %>% sapply(function(x){sum(y[x])})
  Ez <- overlap$cluster_list %>% sapply(function(x){sum(E[x])})
  log_BF_z <- coeff(y, E, shape, rate, overlap$cluster_list)
  BF_z <- exp(log_BF_z)
  log_post_z <- log_prior_z + log_BF_z
  # Do NOT normalize this quantity
  post_z <- exp(log_post_z)
  
  # Generate MCMC samples from prior
  post_chain <- MCMC_simulation(n.sim.post, pattern, theta.init, overlap, 
                                as.matrix(cluster_coords), p_moves, J, post_z, 
                                lambda)
  
  # Trim burn-in
  burnin <- n.sim.post * burnin.prop
  post_sample <- post_chain$sample[-c(1:burnin)]
  
  # Prior Probs of Cluster Membership for each area
  param_post_zone <- list(shape=c(yz+shape[2]), rate=c(Ez+rate[2]))
  RR_post_area <- (y+shape[1])/(E+rate[1])
  post_map <- process_MCMC_sample(post_sample, param_post_zone, RR_post_area, 
                                  overlap$cluster_list, cutoffs)
  
  # Posterior probs of k cluster/anti-clusters for k=0,...,J
  k_vector <- post_chain$sample %>% sapply(length) %>% table()
  k_names <- names(k_vector) %>% as.numeric()
  k_vector <- normalize(k_vector)
  pk_y <- rep(0, J+1)
  pk_y[k_names+1] <- k_vector
  print(paste("Posterior estimation complete on:", date()))
  
  # Output
  output <- list(prior_map=prior_map, post_map=post_map, pk_y=pk_y)
  return(output)
}





#' Estimate lambda values
#' 
#' Internal function to estimate values of lambda needed for 
#' \code{\link{MCMC_simulation}} and prior probability of k 
#' clusters/anti-clusters for k=0,...,J
#'
#' @param n_sim number of importance sampling iterations
#' @param maximum number of clusters/anti-clusters to consider
#' @param prior_z prior probability of each single zone
#' @param overlap output of \code{\link{create_geo_objects}}: list with two 
#' elements: \code{presence} which lists for each area all the single zones it 
#' is present in and \code{cluster_list} for each single zone its component areas
#' @param pi0 prior probability of no clusters
#' @return Estimates of lambda and prior.j
#' @export
#' @seealso \code{\link{MCMC_simulation}}
#' @references Wakefield J. and Kim A.Y. (2013) A Bayesian model for cluster 
#' detection. \emph{Biostatistics}, \bold{14}, 752--765.
#' @examples
#' 1+1
estimate_lambda <- function(n_sim, J, prior_z, overlap, pi0){
  n_zones <- length(prior_z)
  
  # q_0 = q_1 = 1, since there is no concept of overlap for j=0, 1
  q <- c(1, 1, rep(0, J-1))
  var_q <- c(NA, NA, rep(0, J-1))
  
  
  #-------------------------------------------------------------------------------
  # All k's, sample configurations of k single zones, and see if the k single
  # zones overlap or not
  #-------------------------------------------------------------------------------
  for(k in 2:J){
    unit_zones <- sample(1:n_zones, k*n_sim, replace=TRUE, prob=prior_z)
    unit_zones <- matrix(unit_zones, ncol=k)
    no_overlap <- check_overlap(unit_zones, overlap)
    
    denom <- rep(factorial(k), n_sim)
    
    # Compute multinomial coefficient watching out for when certain single zones
    # are sampled more than once
    duplicates <- apply(unit_zones, 1, function(x){length(unique(x))})
    duplicates <- which(duplicates != k)
    for(i in duplicates) {
      denom[i] <-  factorial(k)/prod(factorial(table(unit_zones[i,])))
    }
    
    q[k+1] <- mean(no_overlap/denom)
    var_q[k+1] <- var(no_overlap/denom)/n_sim
  }
  
  
  #-------------------------------------------------------------------------------
  # Using q compute estimate of lambda and prior probability of j clusters
  #-------------------------------------------------------------------------------
  lambda <- (1-pi0)/( (1-pi0)*J + pi0*sum(q[-1]) )
  lambda <- rep(lambda, J)
  lambda <- c(1-sum(lambda), lambda)
  prior_j <- normalize(lambda*q)
  
  # Output
  output <- list(lambda=lambda, prior_j=prior_j)
  return(output)
}



#' Process Markov Chain Monte Carlo Sample
#' 
#' Take the output of sampled configurations from \code{\link{MCMC_simulation}} 
#' and produce area-by-area summaries
#'
#' @param sample list objects of sampled configurations
#' @param param mean relative risk associted with each of the \code{n.zones} 
#' single zones considering the wide prior
#' @param RR_area mean relative risk associated with each of the \code{n} areas 
#' considering the narrow prior
#' @param cluster_list list of length \code{n.zones} listing, for each single 
#' zone, its component areas
#' @param cutoffs cutoffs used to declare highs (clusters) and lows (anti-clusters)
#'
#' @return A list containing
#' \item{high.area}{Probability of cluster membership for each area}
#' \item{low.area}{Probability of anti-cluster membership for each area}
#' \item{RR.est.area}{Smoothed relative risk estimates for each area}
#' @references Wakefield J. and Kim A.Y. (2013) A Bayesian model for cluster 
#' etection. 
#' \emph{Biostatistics}, \bold{14}, 752--765.
#' @export
#' @seealso \code{\link{MCMC_simulation}}
#' @examples
#' 1+1
process_MCMC_sample <- function(sample, param, RR_area, cluster_list, cutoffs){
  n <- length(RR_area)
  n_sim <- length(sample)
  
  # Store outputs here. These are averages weighted by the probability of 
  # each configuration
  high_area <- low_area <- RR_est_area <- rep(0, n)
  
  # Zone based measures
  high_z <- pgamma(cutoffs$high, param$shape, param$rate, lower.tail=FALSE)
  low_z <- pgamma(cutoffs$low, param$shape, param$rate, lower.tail=TRUE)
  RR_z <- param$shape / param$rate
  
  #-------------------------------------------------------------------------------
  # Find probabilities of cluster and anti-cluster membership & estimated 
  # relative risk
  #-------------------------------------------------------------------------------
  # Get probabilities of each configuration
  unique <- unique(sample)
  indices <- match(sample, unique)
  prob <- table(indices)/n_sim
  prob_names <- names(prob) %>% as.numeric()
  
  # Over all sampled configurations
  for(i in 1:length(prob_names)){
    config <- unique[[prob_names[i]]]
    k <- length(config)
    
    if(k==0){
      RR_est_area <- RR_est_area + RR_area * prob[i]
    }else{
      # Used to identify all areas included in configuration
      config_areas <- NULL    
      for(j in 1:k){        
        # All areas in this single zone in this configuration
        zone <- config[j]
        areas <- cluster_list[[zone]]
        config_areas <- c(config_areas, areas)
        
        # Update all areas inside configuration
        high_area[areas] <- high_area[areas] + high_z[zone] * prob[i]
        low_area[areas] <- low_area[areas] + low_z[zone] * prob[i]
        RR_est_area[areas] <- RR_est_area[areas] + RR_z[zone] * prob[i] 
      }
      
      # Update all areas outside configuration
      outside_areas <- setdiff(1:n, config_areas)
      RR_est_area[outside_areas] <- 
        RR_est_area[outside_areas] + RR_area[outside_areas] * prob[i]
    }
  }
  
  # Output
  output <- list(high_area=high_area, low_area=low_area, RR_est_area=RR_est_area)
  return(output)
}









#' Compute Parameters to Calibrate a Gamma Distribution
#' 
#' Compute parameters to calibrate the prior distribution of a relative risk 
#' that has a gamma distribution.
#'
#' @param theta upper quantile
#' @param prob upper quantile
#' @param d mode
#'
#' @return A list containing
#' \item{a}{shape parameter}
#' \item{b}{rate parameter}
#' @export
#' @author Jon Wakefield
#' @seealso \code{\link{LogNormalPriorCh}}
#' @examples
#' param <- GammaPriorCh(5, 0.975, 1)
#' curve(dgamma(x,shape=param$a,rate=param$b),from=0,to=6,n=1000,ylab="density")
GammaPriorCh <- function(theta, prob, d){
  a <- d/2
  b <- 0.5*2*a*theta^2/qt(p=prob,df=2*a)^2
  cat("Gamma Parameters: ", a, b, "\n")
  list(a=a, b=b)
}



#' Compute Parameters to Calibrate a Log-normal Distribution
#' 
#' Compute parameters to calibrate the prior distribution of a relative risk 
#' that has a log-normal distribution.
#'
#' @param theta1 lower quantile
#' @param theta2 upper quantile
#' @param prob1 lower probability
#' @param prob2 upper probability
#'
#' @return A list containing:
#' \item{mu}{mean of log-normal distribution}
#' \item{sigma}{variance of log-normal distribution}
#' @export
#' @author Jon Wakefield
#' @seealso \code{\link{GammaPriorCh}}
#' @examples
#' # Calibrate the log-normal distribution s.t. the 95 percent confidence 
#' interval is [0.2, 5]
#' param <- LogNormalPriorCh(0.2, 5, 0.025, 0.975)
#' curve(dlnorm(x,param$mu,param$sigma), from=0, to=6, ylab="density")
LogNormalPriorCh <- function(theta1, theta2, prob1, prob2){	
  zq1 <- qnorm(prob1)
  zq2 <- qnorm(prob2)
  mu <- log(theta1)*zq2/(zq2-zq1) - log(theta2)*zq1/(zq2-zq1)
  sigma <- (log(theta1)-log(theta2))/(zq1-zq2)
  list(mu=mu,sigma=sigma)
}