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


EBpostthresh <-
  function(Y, E, alpha, beta, Xrow=NULL, rrthresh){
    if (is.null(Xrow)) Xrow <- matrix(rep(1,length(Y)),nrow=length(Y),ncol=1)
    mu <- as.numeric(exp(Xrow %*% beta))
    thresh <- 1-pgamma(rrthresh,alpha+Y,(alpha+E*mu)/mu)
    return(thresh=thresh)
  }
