
globalVariables(c(
  "qgamma"
))



#' Empirical Bayes Estimates of Relative Risk
#' @description The computes empirical Bayes estimates of relative risk of study region with `n` areas, given observed and expected numbers of counts of disease and covariate information.
#' 
#' @param Y a length `n` vector of observed cases
#' @param E a length `n` vector of expected number of cases
#' @param Xmat `n x p` dimension matrix of covariates
#'
#' @references Clayton D. and Kaldor J. (1987) Empirical Bayes estimates of age-standardized relative risks for use in disease mapping.  *Biometrics*, **43**, 671--681
#' 
#' 
#' @return
#' A list with 5 elements:
#'  \item{RR}{the ecological relative risk posterior mean estimates}
#'  \item{RRmed}{the ecological relative risk posterior median estimates}
#'  \item{beta}{the MLE's of the regression coefficients}
#'  \item{alpha}{the MLE of negative binomial dispersion parameter}
#'  \item{SMR}{the standardized mortality/morbidity ratio Y/E}
#'  
#' @export
#' @importFrom MASS glm.nb
#' @examples 
#' data(scotland)
#' data <- scotland$data
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
