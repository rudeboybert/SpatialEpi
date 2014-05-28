estimate.lambda <-
function(pi0, q){
  K <- length(q) - 1
  
  lambda <- (1-pi0)/( (1-pi0)*K + pi0*sum(q[-1]) )
  lambda <- rep(lambda, K)
  lambda <- c(1-sum(lambda), lambda)

  return(lambda)
}

