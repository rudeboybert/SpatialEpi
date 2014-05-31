MCMC <-
function(n.sim, overlap, cluster.coords, K, values, lambda, theta.init){
  pattern <- c(0,1)
  p.moves.orig <- normalize(c(1,1,1,1,1))
  
  chain <- MCMC_simulation(
    n.sim, pattern, theta.init, overlap, cluster.coords, p.moves.orig, K, 
    values, lambda)
  
  return(chain)
}
