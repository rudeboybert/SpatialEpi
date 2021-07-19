
globalVariables(c(
  "qt"
))



#' Compute Parameters to Calibrate a Gamma Distribution
#'
#' 
#' @description  Compute parameters to calibrate the prior distribution of a relative risk that has a gamma distribution.
#' @param theta upper quantile
#' @param prob upper quantile
#' @param d degrees of freedom
#' 
#' 
#' @author Jon Wakefield
#' @seealso LogNormalPriorCh
#' 
#' @return
#' List containing
#' \item{a}{shape parameter}
#' \item{b}{rate parameter}
#' 
#' @export

#' @examples 
#' param <- GammaPriorCh(5, 0.975,1)
#' curve(dgamma(x,shape=param$a,rate=param$b),from=0,to=6,n=1000,ylab="density")
#' 
GammaPriorCh <-
function(theta, prob, d){
	a <- d/2
	b <- 0.5*2*a*theta^2/qt(p=prob,df=2*a)^2
	cat("Gamma Parameters: ",a,b,"\n")
	list(a=a,b=b)
}

