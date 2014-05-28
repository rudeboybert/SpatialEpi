process.MCMC.chain <-
function(chain, param, RR.area, cluster.list, cutoffs){
  n <- length(RR.area)
  # Store Outputs Here
  high.area <- low.area <- RR.est.area <- rep(0, n)

  #--------------------------------------------------------------
  # Zone Based Measures
  #--------------------------------------------------------------
  high.z <- pgamma(cutoffs$high, param$shape, param$rate, lower.tail=FALSE)
  low.z <- pgamma(cutoffs$low, param$shape, param$rate, lower.tail=TRUE)
  RR.z <- param$shape / param$rate

  
  #--------------------------------------------------------------
  # Find Probability of Cluster Membership & Estimated Relative Risk
  #--------------------------------------------------------------
  #-----------------------------------------------
  # Get probabilities via frequencies
  #-----------------------------------------------
  n.sim <- length(chain$sample)

  unique <- unique(chain$sample)
  indices <- match(chain$sample, unique)
  
  prob <- table(indices)/n.sim
  prob.names <- as.numeric(names(prob))

  
  #-----------------------------------------------
  # Over all sampled configurations
  #-----------------------------------------------
  for(i in 1:length(prob.names)){
    config <- unique[[prob.names[i]]]
    k <- length(config)
    
    if(k==0){
      RR.est.area <- RR.est.area + RR.area * prob[i]
    }else{
      # Used to identify all areas included in configuration
      config.areas <- NULL
      
      for(j in 1:k){        
        # All areas in this single zone
        zone <- config[j]
        areas <- cluster.list[[zone]]
        config.areas <- c(config.areas,areas)

        # Update outputs
        high.area[areas] <- high.area[areas] + high.z[zone] * prob[i]
        low.area[areas] <- low.area[areas] + low.z[zone] * prob[i]
        RR.est.area[areas] <- RR.est.area[areas] + RR.z[zone] * prob[i] 
      }
      
      # All areas outside configuration
      outside.areas <- setdiff(1:n, config.areas)
      RR.est.area[outside.areas] <- RR.est.area[outside.areas] + RR.area[outside.areas]*prob[i]
    }
  }
  
  return(list(high.area=high.area, low.area=low.area, RR.est.area=RR.est.area))
}

