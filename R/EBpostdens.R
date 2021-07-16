
globalVariables(c(
  "dgamma","title","lines"
))



#' @title Produce plots of empirical Bayes posterior densities when the data Y are Poisson with expected number E and relative risk theta, with the latter having a gamma distribution with known values alpha and beta, which are estimated using empirical Bayes.
#'
#' @description This function produces plots of empirical Bayes posterior densities which are gamma distributions with parameters (alpha+Y, (alpha+E*mu)/mu) where mu = exp(x beta). The SMRs are drawn on for comparison.
#'
#' @author Jon Wakefield
#'
#' @param Y observed disease counts
#' @param E expected disease counts
#' @param alpha x
#' @param beta x
#' @param Xrow x
#' @param lower x
#' @param upper x
#' @param main x
#'
#' @return
#' A plot containing the gamma posterior distribution 
#' 
#' @export
#'
#' @examples
#' data(scotland)
#' Y <- scotland$data$cases
#' E <- scotland$data$expected
#' ebresults <- eBayes(Y,E)
#' EBpostdens(Y[1], E[1], ebresults$alpha, ebresults$beta, lower=0, upper=15,
#'           main="Area 1")
EBpostdens <-
function(Y, E, alpha, beta, Xrow=NULL, lower=NULL, upper=NULL, main=""){
  if (is.null(Xrow)) Xrow <- matrix(c(1),nrow=1,ncol=1)
  xvals <- seq(lower,upper,.01)
  mu <- as.numeric(exp(Xrow %*% beta))
  yvals <- dgamma(xvals,alpha+Y,(alpha+E*mu)/mu)
  plot(xvals,yvals,type="n",
       xlab=expression(theta),ylab="EB density",cex.lab=1.2)
  title(paste(main))
  lines(xvals,yvals)
  lines(c(Y/E,Y/E),c(0,max(yvals)),lty=2)
}
