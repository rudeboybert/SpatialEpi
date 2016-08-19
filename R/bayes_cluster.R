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






estimate_lambda <-
  function(n.sim, J, prior.z, overlap, pi0){
    
    n.zones <- length(prior.z)
    
    # q_0 = q_1 = 1, since there is no concept of overlap for j=0, 1
    q <- c(1, 1, rep(0, J-1))
    var.q <- c(NA, NA, rep(0, J-1))
    
    
    #-------------------------------------------------------------------------------
    # All k's, sample configurations of k single zones, and see if the k single
    # zones overlap or not
    #-------------------------------------------------------------------------------
    for(k in 2:J){
      unit.zones <- sample(1:n.zones, k*n.sim, replace=TRUE, prob=prior.z)
      unit.zones <- matrix(unit.zones, ncol=k)
      no.overlap <- check_overlap(unit.zones, overlap)
      
      denom <- rep(factorial(k), n.sim)
      
      # Compute multinomial coefficient watching out for when certain single zones
      # are sampled more than once
      duplicates <- apply(unit.zones, 1, function(x){length(unique(x))})
      duplicates <- which(duplicates != k)
      for(i in duplicates) {
        denom[i] <-  factorial(k)/prod(factorial(table(unit.zones[i,])))
      }
      
      q[k+1] <- mean(no.overlap/denom)
      var.q[k+1] <- var(no.overlap/denom)/n.sim
    }
    
    
    #-------------------------------------------------------------------------------
    # Using q compute estimate of lambda and prior probability of j clusters
    #-------------------------------------------------------------------------------
    lambda <- (1-pi0)/( (1-pi0)*J + pi0*sum(q[-1]) )
    lambda <- rep(lambda, J)
    lambda <- c(1-sum(lambda), lambda)
    prior.j <- normalize(lambda*q)
    
    return(list(lambda=lambda, prior.j=prior.j))
  }




process_MCMC_sample <-
  function(sample, param, RR.area, cluster.list, cutoffs){
    
    n <- length(RR.area)
    n.sim <- length(sample)
    
    # Store outputs here. These are averages weighted by the probability of 
    # each configuration
    high.area <- low.area <- RR.est.area <- rep(0, n)
    
    # Zone based measures
    high.z <- pgamma(cutoffs$high, param$shape, param$rate, lower.tail=FALSE)
    low.z <- pgamma(cutoffs$low, param$shape, param$rate, lower.tail=TRUE)
    RR.z <- param$shape / param$rate
    
    
    #-------------------------------------------------------------------------------
    # Find probabilities of cluster and anti-cluster membership & estimated 
    # relative risk
    #-------------------------------------------------------------------------------
    # Get probabilities of each configuration
    unique <- unique(sample)
    indices <- match(sample, unique)
    prob <- table(indices)/n.sim
    prob.names <- as.numeric(names(prob))
    
    # Over all sampled configurations
    for(i in 1:length(prob.names)){
      config <- unique[[prob.names[i]]]
      k <- length(config)
      
      if(k==0){
        RR.est.area <- RR.est.area + RR.area * prob[i]
      }else{
        # Used to identify all areas included in configuration
        config.areas <- NULL    
        for(j in 1:k){        
          # All areas in this single zone in this configuration
          zone <- config[j]
          areas <- cluster.list[[zone]]
          config.areas <- c(config.areas, areas)
          
          # Update all areas inside configuration
          high.area[areas] <- high.area[areas] + high.z[zone] * prob[i]
          low.area[areas] <- low.area[areas] + low.z[zone] * prob[i]
          RR.est.area[areas] <- RR.est.area[areas] + RR.z[zone] * prob[i] 
        }
        
        # Update all areas outside configuration
        outside.areas <- setdiff(1:n, config.areas)
        RR.est.area[outside.areas] <- 
          RR.est.area[outside.areas] + RR.area[outside.areas]*prob[i]
      }
    }
    
    return(list(
      high.area=high.area,
      low.area=low.area,
      RR.est.area=RR.est.area)
    )
  }



















GammaPriorCh <-
  function(theta, prob, d){
    a <- d/2
    b <- 0.5*2*a*theta^2/qt(p=prob,df=2*a)^2
    cat("Gamma Parameters: ",a,b,"\n")
    list(a=a,b=b)
  }


LogNormalPriorCh <-
  function(theta1, theta2, prob1, prob2){	
    zq1 <- qnorm(prob1)
    zq2 <- qnorm(prob2)
    mu <- log(theta1)*zq2/(zq2-zq1) - log(theta2)*zq1/(zq2-zq1)
    sigma <- (log(theta1)-log(theta2))/(zq1-zq2)
    list(mu=mu,sigma=sigma)
  }



