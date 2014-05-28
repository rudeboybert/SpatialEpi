MCMC <-
function(n.sim, overlap, cluster.coords, J, values,
                 lambda, theta.init){
  pattern <- c(0,1)
  p.moves.orig <- normalize(c(1,1,1,1,1))

  chain <- .Call("MCMC", n.sim, pattern, theta.init, overlap,
                  cluster.coords, p.moves.orig, J, values,
                  lambda, PACKAGE="SpatialEpi")
  return(chain)
}

