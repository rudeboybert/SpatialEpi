
globalVariables(c(
  "pgamma"
))


#' Produce the probabilities of exceeding a threshold given a posterior gamma distribution.
#'
#' @description This function produces the posterior probabilities of exceeding a threshold given a gamma distributions with parameters (alpha+Y, (alpha+E*mu)/mu) where mu = exp(x beta). This model arises from Y being Poisson with mean theta times E where theta is the relative risk and E are the expected numbers. The prior on theta is gamma with parameters alpha and beta. The parameters alpha and beta may be estimated using empirical Bayes.
#' @author Jon Wakefield
#' @seealso  \code{\link{eBayes}}
#' 
#' @param Y observed disease counts
#' @param E expected disease counts
#' @param alpha x
#' @param beta x
#' @param Xrow x
#' @param rrthresh x
#'
#' @return
#' Posterior probabilities of exceedence are returned.
#' 
#' @export
#'
#' @examples
#' data(scotland)
#' Y <- scotland$data$cases
#' E <- scotland$data$expected
#' ebresults <- eBayes(Y,E)
#' #Find probabilities of exceedence of 3
#' thresh3 <- EBpostthresh(Y, E, alpha=ebresults$alpha, beta=ebresults$beta, rrthresh=3)
#' mapvariable(thresh3, scotland$spatial.polygon)

EBpostthresh <-
function(Y, E, alpha, beta, Xrow=NULL, rrthresh){
  if (is.null(Xrow)) Xrow <- matrix(rep(1,length(Y)),nrow=length(Y),ncol=1)
  mu <- as.numeric(exp(Xrow %*% beta))
  thresh <- 1-pgamma(rrthresh,alpha+Y,(alpha+E*mu)/mu)
  return(thresh=thresh)
}
