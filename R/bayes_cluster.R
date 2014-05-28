bayes_cluster <-
function(E, cases, population, centroids, map, max.prop, k, shape, rate, J, pi0,
         n.sim.imp, n.sim.prior, n.sim.post
  ){
    burnin.prop <- 0.1
    theta.init.prior <- NULL
    theta.init.post <- NULL 
    
    #-----------------------------------------------
    # Create Geographical Objects to Use
    #-----------------------------------------------
    geo.objects <- create_geo_objects(max.prop, population, centroids, map)
    overlap <- geo.objects$overlap
    cluster.coords <- geo.objects$cluster.coords
    areaz <- geo.objects$areaz
    
    
    #-----------------------------------------------
    # Set Cutoffs for High and Low:  need to change the range
    #-----------------------------------------------
    theta <- seq(0.90, 1.1, length=50000)
    pN <- dgamma(theta, shape[1], rate[1])
    pW <- dgamma(theta, shape[2], rate[2])
    cutoffs <- list(
      high=theta[max(which(pN/pW>1))],
      low=theta[min(which(pN/pW>1))]
    )
    rm(pN, pW, theta)
    
    
    #-----------------------------------------------
    # prior on single zones
    #-----------------------------------------------
    n <- length(cases)
    n.zones <- length(areaz)
    cluster.list <- overlap$cluster_list
    
    prior.z <- normalize( exp(-k*areaz) )
    log_prior.z <- log(prior.z) - log(sum(prior.z))
    
    
    #-----------------------------------------------
    # Estimate q and lambda via importance sampling
    #-----------------------------------------------
    results <- estimate_q(n.sim.imp, J, prior.z, overlap)
    q <- results$q
    lambda <- estimate_lambda(pi0, q)
    prior.j <- normalize(lambda*q)
    
    
    #-----------------------------------------------
    # Estimate prior maps
    #-----------------------------------------------
    prior.chain.full <- MCMC(n.sim.prior, overlap, cluster.coords, J, prior.z, 
                             lambda, theta.init.prior)    
    # trim burn-in
    prior.chain <- prior.chain.full
    burnin <- n.sim.prior * burnin.prop
    
    prior.chain$sample <- prior.chain$sample[-c(1:burnin)]
    
    # Prior Probs of Cluster Membership for each Area
    param.prior.zone <- list(shape=rep(shape[2],n.zones), rate=rep(rate[2],n.zones))
    RR.prior.area <- rep(shape[1]/rate[1], n)
    prior.map <- process_MCMC_chain(prior.chain, param.prior.zone, 
                                    RR.prior.area, cluster.list, cutoffs)
    
    
    #-----------------------------------------------
    # Estimate posterior maps
    #-----------------------------------------------
    # Unit Zone Values
    yz <- sapply(cluster.list, function(x){sum(cases[x])})
    Ez <- sapply(cluster.list, function(x){sum(E[x])})
    log_BF.z <- coeff(cases, E, shape, rate, cluster.list)
    BF.z <- exp(log_BF.z)
    log_post.z <- log_prior.z + log_BF.z
    post.z <- exp(log_post.z)
    
    # Generate MCMC Chain
    post.chain.full <- MCMC(n.sim.post, overlap, cluster.coords, J, post.z, 
                            lambda, theta.init.prior) 
    
    # trim burn-in
    post.chain <- post.chain.full
    burnin <- n.sim.post * burnin.prop
    
    post.chain$sample <- post.chain$sample[-c(1:burnin)]
    
    # Posterior Probs of j cluster/anti-clusters
    pj.y.orig <- normalize(table(sapply(post.chain$sample,length)))
    pj.y <- rep(0, J+1)
    pj.y[as.numeric(names(pj.y.orig))+1] <- pj.y.orig
    
    # Posterior Probs of Cluster Membership for each Area
    param.post.zone <- list(shape=c(yz+shape[2]), rate=c(Ez+rate[2]))
    RR.post.area <- (cases+shape[1])/(E+rate[1])
    post.map <- process_MCMC_chain(post.chain, param.post.zone, RR.post.area,
                                   cluster.list, cutoffs)
    
    
    return(list(prior.map=prior.map, post.map=post.map, pj.y=pj.y))
  }
