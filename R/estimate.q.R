estimate.q <-
function(n.sim, J, prior.z, overlap){
  n.zones <- length(prior.z)

  # q_0 = q_1 = 1, since no concept of overlap for j=0,1
  q <- c(1, 1, rep(0,J-1))
  var.q <- c(NA, NA, rep(0,J-1))

  for(k in 2:J){
    unit.zones <- sample(1:n.zones, k*n.sim, replace=TRUE, prob=prior.z)
    unit.zones <- matrix(unit.zones, ncol=k)

    if(is.matrix(overlap)){
      indicator <- as.numeric(apply(unit.zones,1,function(x){sum(overlap[x,x])-k==0}))
    }else{
      indicator <- .Call("check_overlap", t(unit.zones), overlap, PACKAGE="SpatialEpi")
    }
  
    denom <- rep(factorial(k),n.sim)
    # Special Multinomial Coefficient Case
    duplicates <- apply(unit.zones,1,function(x){length(unique(x))})
    duplicates <- which(duplicates != k)
    for(i in duplicates)
      denom[i] <-  factorial(k)/prod(factorial(table(unit.zones[i,])))
    
    q[k+1] <- mean(indicator/denom)
    var.q[k+1] <- var(indicator/denom)/n.sim

    print(k) 
  }
  
  return(list(q=q, var.q=var.q))
}

