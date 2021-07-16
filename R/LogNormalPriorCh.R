
globalVariables(c(
  "qnorm"
))



#' Compute Parameters to Calibrate a Log-normal Distribution
#'
#'
#' @description 	Compute parameters to calibrate the prior distribution of a relative risk that has a log-normal distribution.
#' @param theta1 lower quantile
#' @param theta2 upper quantile
#' @param prob1  lower probability
#' @param prob2  upper probability
#' @author Jon Wakefield
#' @return
#'  A list containing
#' \item{mu}{mean of log-normal distribution}
#' \item{sigma}{variance of log-normal distribution}
#' @export
#' @examples 
#' # Calibrate the log-normal distribution s.t. the 95% confidence interval is [0.2, 5]
#' param <- LogNormalPriorCh(0.2, 5, 0.025, 0.975)
#' curve(dlnorm(x,param$mu,param$sigma), from=0, to=6, ylab="density")
LogNormalPriorCh <-
function(theta1, theta2, prob1, prob2){	
zq1 <- qnorm(prob1)
zq2 <- qnorm(prob2)
mu <- log(theta1)*zq2/(zq2-zq1) - log(theta2)*zq1/(zq2-zq1)
sigma <- (log(theta1)-log(theta2))/(zq1-zq2)
list(mu=mu,sigma=sigma)
}