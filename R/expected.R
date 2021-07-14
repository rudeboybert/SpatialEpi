#' Compute Expected Numbers of Disease
#'
#' @description Compute the internally indirect standardized expected numbers of disease.
#' @param population a vector of population counts for each strata in each area
#' @param cases a vector of the corresponding number of cases
#' @param n.strata number of strata considered
#' 
#' @details The \code{population} and \code{cases} vectors must be \emph{balanced}: all counts are sorted by area first, and then within each area the counts for all strata are listed (even if 0 count) in the same order.
#'
#' @references Elliot, P. et al. (2000) \emph{Spatial Epidemiology:  Methods and Applications}.  Oxford Medical Publications.
#' @author Albert Y. Kim
#' @return
#' \item{expected.cases}{a vector of the expected numbers of disease for each area}
#' @export
#'
#' @examples
#' data(pennLC)
#' population <- pennLC$data$population
#' cases <- pennLC$data$cases

#' ## In each county in Pennsylvania, there are 2 races, gender and 4 age bands 
#' ## considered = 16 strata levels
#' pennLC$data[1:16,]
#' expected(population, cases, 16)
#' 
#' 
expected <- 
  function(population, cases, n.strata){
    
n <- length(population)/n.strata
E <- rep(0, n)
qNum <- rep(0, n.strata)
qDenom <- rep(0, n.strata)
q <- rep(0, n.strata)


# Compute q: strata-specific rates. Numerator and denominator separately
for(i in 1:n.strata){
  indices <- rep(i, n) + seq(0, n-1)*n.strata
  qNum[i] <- qNum[i] + sum(cases[indices])
  qDenom[i] <- qDenom[i] + sum(population[indices])
}
q <- qNum/qDenom


# Compute E expected counts
for(i in 1:n) {
  indices <- 1:n.strata + (i-1)*n.strata
  E[i] <- sum(population[indices]*q)
}

return(E)
}





