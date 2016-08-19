#' Title
#'
#' @param y 
#' @param E 
#' @param population 
#' @param sp.obj 
#' @param centroids 
#' @param max.prop 
#' @param shape 
#' @param rate 
#' @param J 
#' @param pi0 
#' @param n_sim.lambda 
#' @param n_sim.prior 
#' @param n_sim_post 
#' @param burnin.prop 
#' @param theta.init 
#' 
#' @importFrom magrittr %>%
#'
#' @return
#' @export
#'
#' @examples
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
  # geo.objects <- create_geo_objects(max.prop, population, centroids, sp.obj)
  # overlap <- geo.objects$overlap
  # cluster_coords <- geo.objects$cluster_coords 
  # single_zones <- define_single_zones(centroids, NYleukemia$data$population, 0.15)
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



#' Title
#'
#' @param n_sim 
#' @param J 
#' @param prior_z 
#' @param overlap 
#' @param pi0 
#'
#' @return
#' @export
#'
#' @examples
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



#' Title
#'
#' @param sample 
#' @param param 
#' @param RR_area 
#' @param cluster_list 
#' @param cutoffs 
#'
#' @return
#' @export
#'
#' @examples
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



#' Title
#'
#' @param theta 
#' @param prob 
#' @param d 
#'
#' @return
#' @export
#'
#' @examples
GammaPriorCh <- function(theta, prob, d){
  a <- d/2
  b <- 0.5*2*a*theta^2/qt(p=prob,df=2*a)^2
  cat("Gamma Parameters: ",a,b,"\n")
  list(a=a,b=b)
}



#' Title
#'
#' @param theta1 
#' @param theta2 
#' @param prob1 
#' @param prob2 
#'
#' @return
#' @export
#'
#' @examples
LogNormalPriorCh <- function(theta1, theta2, prob1, prob2){	
  zq1 <- qnorm(prob1)
  zq2 <- qnorm(prob2)
  mu <- log(theta1)*zq2/(zq2-zq1) - log(theta2)*zq1/(zq2-zq1)
  sigma <- (log(theta1)-log(theta2))/(zq1-zq2)
  list(mu=mu,sigma=sigma)
}