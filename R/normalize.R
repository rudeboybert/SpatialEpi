normalize <-
function (p, include.indices=NA){
  p.prime <- p
  p.prime[!include.indices] <- 0
  p.prime <- p.prime/sum(p.prime)
  return(p.prime)
}

