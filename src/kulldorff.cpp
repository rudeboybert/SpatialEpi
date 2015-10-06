#include <Rcpp.h>
using namespace Rcpp;

//' Compute binomial likelihood ratio test statistic for Kulldorff method.
//' 
//' @param cz Disease count inside single zone.
//' @param nz Expected count inside single zone.
//' @param C Total disease count in region.
//' @param N Total expected count in region.
//' @return log of likelihood ratio test statistic.
//' @export
// [[Rcpp::export]]
double binomialLogLkhd(double cz, double nz, double C, double N) {
  double logLkhd = 0;
  
  if(cz/nz <= (C-cz)/(N-nz)) {
    logLkhd = 0;
  } else {
    logLkhd  =	
      N * (log(N-nz-C+cz) - log(N-C) + log(N) - log(N-nz)) +
      nz * (log(nz-cz) - log(N-nz-C+cz) + log(N-nz) - log(nz))+
      cz * (log(cz) - log(nz-cz) + log(N-nz-C+cz) - log(C-cz))+
      C * (log(C-cz) - log(C) + log(N-C) - log(N-nz-C+cz));
  }
  
  return logLkhd;
}

//' Compute poisson likelihood ratio test statistic for Kulldorff method.
//' 
//' @inheritParams binomialLogLkhd
//' @return log of likelihood ratio test statistic.
//' @export
// [[Rcpp::export]]
double poissonLogLkhd(double cz, double nz, double C, double N) {
  double logLkhd = 0;
  
  if(cz/nz <= (C-cz)/(N-nz)) {
    logLkhd = 0;
  } else {
    logLkhd =
      cz * log(cz / nz) +
      cz * log((N-nz)/(C-cz)) +
      C * log((C-cz)/(N-nz)) +
      C * log(N/C);
  }  
  return logLkhd;
}

//' Compute log likelihood ratio test statistic for all single zones.
//' 
//' @param observedCases
//' @param expectedCases
//' @param nearestNeighborsList
//' @param nZones
//' @param logLkhdType
//' @return Vector of log likelihood
//' @export
// [[Rcpp::export]]
NumericVector computeAllLogLkhd(NumericVector observedCases, 
                                NumericVector expectedCases, 
                                List nearestNeighborsList, 
                                String logLkhdType) {
  int nAreas = expectedCases.size(), 
    C = sum(observedCases), 
    N = sum(expectedCases), 
    incrementor = 0,
    nNeighbors = 0,
    nZones = 0;
  double cz = 0.0, 
    nz = 0.0;

  // Compute number of single zones
  // i.e. sapply(nearestNeighborsList, length) %>% sum()
  for (int i = 0; i < nAreas; ++i) {
    Rcpp::NumericVector nearestNeighbors = nearestNeighborsList[i];
    nZones += nearestNeighbors.size();
  }
  
  NumericVector allLogLkhd(nZones);
  
  // Loop through single zones and compute log-lkhd
  for (int i = 0; i < nAreas; ++i) {
    cz = 0.0;
    nz = 0.0;
    
    Rcpp::NumericVector nearestNeighbors = nearestNeighborsList[i];
    nNeighbors = nearestNeighbors.size();
    
    // For each area's nearest neighbors
    for(int j = 0; j < nNeighbors; ++j) { 
      // Watch off by 1 vector indexing in C++ as opposed to R
      cz += observedCases[nearestNeighbors[j]-1];
      nz += expectedCases[nearestNeighbors[j]-1];
      
      if (logLkhdType == "poisson") {
        allLogLkhd[incrementor] = poissonLogLkhd(cz, nz, C, N);
      } else if (logLkhdType == "binomial") {
        allLogLkhd[incrementor] = binomialLogLkhd(cz, nz, C, N);
      }
      incrementor++;
    }
  }
  return allLogLkhd;
}

//' Compute log likelihood ratio test statistic for all 
//' 
//' @inheritParams computeAllLogLkhd
//' @param permutedCaseMatrix blah
//' @return Max
//' @export
// [[Rcpp::export]]
NumericVector kulldorffMC(NumericMatrix permutedCaseMatrix, 
                          NumericVector expectedCases, 
                          List nearestNeighborsList, 
                          String logLkhdType) {
  
  // Define variables
  int nAreas = permutedCaseMatrix.nrow(),
    nSims = permutedCaseMatrix.ncol(),
    nZones = 0;
  
  // Compute number of single zones
  // i.e. sapply(nearestNeighborsList, length) %>% sum()
  for (int i = 0; i < nAreas; ++i) {
    Rcpp::NumericVector nearestNeighbors = nearestNeighborsList[i];
    nZones += nearestNeighborsList.size();
  }
  
  NumericVector allLogLkhd(nZones);
  NumericVector permutedCases(nAreas); 
  NumericVector maxLogLkhd(nSims);
  
  /*
   Main loop:  For each simulation, compute all the log-likelihood ratio
   statistics for each zone and return the max
   */
  for(int j=0; j < nSims; ++j) {
    for(int i = 0; i < nAreas; ++i) {
      permutedCases[i] = permutedCaseMatrix(i,j);
    }
    
    // compute max log likelihood
    allLogLkhd = computeAllLogLkhd(permutedCases, expectedCases, 
                                   nearestNeighborsList, logLkhdType);
    maxLogLkhd[j] =  max(allLogLkhd);
  }
  
  // Output results
  return maxLogLkhd;
}