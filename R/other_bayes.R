#' Empirical Bayes Estimates of Relative Risk
#'
#' The computes empirical Bayes estimates of relative risk of study region with 
#' \code{n} areas, given observed and expected numbers of counts of disease and 
#' covariate information.  
#'
#' @param Y a length \code{n} vector of observed cases
#' @param E a length \code{n} vector of expected number of cases
#' @param Xmat \code{n x p} dimension matrix of covariates
#'
#' @return A list with 5 elements:
#' \item{RR}{the ecological relative risk posterior mean estimates}
#' \item{RRmed}{the ecological relative risk posterior mean estimates}
#' \item{beta}{the MLE's of the regression coefficients}
#' \item{alpha}{the MLE of negative binomial dispersion parameter}
#' \item{SMR}{the standardized mortality/morbidity ratio Y/E}
#' @export
#' @seealso \code{\link{scotland}}
#' @seealso \code{\link{mapvariable}}
#' @references Clayton D. and Kaldor J. (1987) Empirical Bayes estimates of 
#' age-standardized relative risks for use in disease mapping.\emph{Biometrics},
#' \bold{43}, 671--681
#' @examples
#' data(scotland)
#' data <- scotland$data
#' 
#' x <- data$AFF
#' Xmat <- cbind(x,x^2)
#' results <- eBayes(data$cases,data$expected,Xmat)
#' scotland.map <- scotland$spatial.polygon
#' mapvariable(results$RR, scotland.map)
eBayes <-
  function(Y, E, Xmat=NULL){
    # Check for covariates
    if (is.null(Xmat)) {
      mod <- glm.nb(Y ~ 1+offset(log(E)), link=log)
    }else{
      mod <- glm.nb(Y ~ Xmat+offset(log(E)), link=log)
    }
    
    # Get MLE's of parameters
    alpha <- mod$theta
    muhat <- mod$fitted/E
    
    # Generate weights between global and local SMR
    wgt <- E*muhat/(alpha+E*muhat); SMR <- Y/E
    RR <- as.numeric(wgt*SMR + (1-wgt)*muhat)
    RRmed <- qgamma(0.5,alpha+Y,(alpha+E*muhat)/muhat)
    
    # Output results
    list(RR=RR, RRmed=RRmed, beta=mod$coeff, alpha=alpha, SMR=SMR)
  }


#' Produce plots of emprical Bayes posterior densities when the data Y are 
#' Poisson with expected number E and relative risk theta, with the latter 
#' having a gamma distribution with known values alpha and beta, which are 
#' estimated using empirical Bayes.
#' 
#' This function produces plots of empirical Bayes posterior densities which are 
#' gamma distributions with parameters (alpha+Y, (alpha+E*mu)/mu) where 
#' mu = exp(x beta). The SMRs are drawn on for comparison.
#'
#' @param Y observed disease counts
#' @param E expected disease counts
#' @param alpha 
#' @param beta 
#' @param Xrow 
#' @param lower 
#' @param upper 
#' @param main 
#'
#' @return A plot containing the gamma posterior distribution 
#' @export
#' @author Jon Wakefield
#' @seealso \code{\link{EBpostthresh}}
#' @seealso \code{\link{eBayes}}
#' @examples
#' data(scotland)
#' Y <- scotland$data$cases
#' E <- scotland$data$expected
#' ebresults <- eBayes(Y,E)
#' EBpostdens(Y[1], E[1], ebresults$alpha, ebresults$beta, lower=0, upper=15, main="Area 1")
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


#' Produce the probabilities of exceeding a threshold given a posterior gamma
#' distribution.
#' 
#' This function produces the posterior probabilities of exceeding a threshold 
#' given a gamma distributions with parameters (alpha+Y, (alpha+E*mu)/mu) where 
#' mu = exp(x beta). This model arises from Y being Poisson with mean theta 
#' times E where theta is the relative risk and E are the expected numbers. The 
#' prior on theta is gamma with parameters alpha and beta. The parameters alpha 
#' and beta may be estimated using empirical Bayes. 
#'
#' @param Y observed disease counts
#' @param E expected disease counts
#' @param alpha 
#' @param beta 
#' @param Xrow 
#' @param rrthresh 
#'
#' @author Jon Wakefield
#' @return Posterior probabilities of exceedence are returned.
#' @export
#' @seealso \code{\link{eBayes}}
#'
#' @examples
#' data(scotland)
#' Y <- scotland$data$cases
#' E <- scotland$data$expected
#' ebresults <- eBayes(Y,E)
#' # Find probabilities of exceedence of 3
#' thresh3 <- EBpostthresh(Y, E, alpha=ebresults$alpha, beta=ebresults$beta, rrthresh=3)
#' mapvariable(thresh3, scotland$spatial.polygon)
EBpostthresh <-
  function(Y, E, alpha, beta, Xrow=NULL, rrthresh){
    if (is.null(Xrow)) Xrow <- matrix(rep(1,length(Y)),nrow=length(Y),ncol=1)
    mu <- as.numeric(exp(Xrow %*% beta))
    thresh <- 1-pgamma(rrthresh,alpha+Y,(alpha+E*mu)/mu)
    return(thresh=thresh)
  }
