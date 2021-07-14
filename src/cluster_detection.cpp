#include <Rcpp.h>
//[[Rcpp::plugins(cpp11)]]
using namespace Rcpp;


// [[Rcpp::export]]
double binomialLogLkhd(double cz, double nz, double N, double C) {
	double logLkhd = 0;

	if(cz / nz <= (C - cz)/(N - nz)) {
		logLkhd = 0;
	} else {
		logLkhd  =
		N * ( log( N-nz-C+cz ) - log( N-C )  + log( N ) - log( N-nz ) ) +
		nz * ( log( nz-cz ) - log( N-nz-C+cz )  + log( N-nz ) - log( nz ) )+
		cz * ( log( cz ) - log( nz-cz ) + log( N-nz-C+cz ) - log( C-cz ) )+
		C * ( log( C-cz ) - log( C ) + log( N-C ) - log( N-nz-C+cz ) );
	}

	return logLkhd;
}


// [[Rcpp::export]]
double poissonLogLkhd(double cz, double nz, double N, double C) {
	double logLkhd = 0;

	if(cz / nz <= (C - cz)/(N - nz)) {
		logLkhd = 0;
	} else {
		logLkhd =
		cz * log(  (cz / nz) )  +
		cz * log(  ( (N - nz)/( C - cz) ) )  +
		C * log(  ( (C-cz)/(N-nz) )  ) +
		C * log(  ( N/ C )  );
	}
	return logLkhd;
}


// [[Rcpp::export]]
NumericVector computeAllLogLkhd(NumericVector observedCases, NumericVector expectedCases, List nearestNeighborsList, int nZones, String logLkhdType) {

	NumericVector allLogLkhd(nZones);
	int nAreas = expectedCases.size(), C = sum(observedCases),
	N = sum(expectedCases), index = 0;
	int nNeighbors = 0;
	double cz = 0.0, nz = 0.0;

	for (int i = 0; i < nAreas; ++i) {
		cz = 0;
		nz = 0;

		Rcpp::NumericVector nearestNeighbors = nearestNeighborsList[i];
		nNeighbors =  nearestNeighbors.size();

		// For each area's nearest neighbors
		for(int j = 0; j < nNeighbors; ++j) {
			// Watch off by 1 vector indexing in C as opposed to R
			cz += observedCases[nearestNeighbors[j]-1];
			nz += expectedCases[nearestNeighbors[j]-1];

			if (logLkhdType=="poisson") {
				allLogLkhd[index] = poissonLogLkhd(cz, nz, N, C);
			} else if (logLkhdType == "binomial" ) {
				allLogLkhd[index] = binomialLogLkhd(cz, nz, N, C);
			}
			index++;
		}
	}
	return allLogLkhd;
}



// [[Rcpp::export]]
NumericVector kulldorffMC(NumericMatrix permutedCaseMatrix, NumericVector expectedCases, List nearestNeighbors, int nZones, String logLkhdType) {

	// Define variables
	int nAreas = permutedCaseMatrix.nrow();
	int nSimulations = permutedCaseMatrix.ncol();

	NumericVector allLogLkhd(nZones);
	NumericVector permutedCases(nAreas);
	NumericVector maxLogLkhd(nSimulations);

	/*
	Main loop:  For each simulation, compute all the log-likelihood ratio
	statistics for each zone and then return the max
	*/
	for(int j=0; j < nSimulations; ++j) {
		// Load simulatedCases
		for(int i = 0; i < nAreas; ++i) {
			permutedCases[i] = permutedCaseMatrix(i,j);
		}

		// compute max log likelihood
		allLogLkhd = computeAllLogLkhd(permutedCases, expectedCases, nearestNeighbors,
			nZones, logLkhdType);
			maxLogLkhd[j] =  max(allLogLkhd);
		}

		// Output results
		return maxLogLkhd;
	}



	// [[Rcpp::export]]
	List besag_newell_internal(NumericVector observedCases, NumericVector expectedCases, List nearestNeighborsList, int nZones, int k) {

		double sumObserved, sumExpected;
		int nAreas = observedCases.size(), nNeighbors = 0;
		int j = 0;
		NumericVector x(1);
		NumericVector observedPValues(nAreas), observedMValues(nAreas),
		observedKValues(nAreas);

		for(int i=0; i < nAreas; ++i) {
			sumObserved = 0.0;
			sumExpected = 0.0;

			Rcpp::NumericVector nearestNeighbors = nearestNeighborsList[i];
			nNeighbors =  nearestNeighbors.size();

			for(j = 0; j < nNeighbors; ++j) {
				// Watch off by 1 vector indexing in C as opposed to R
				sumObserved += observedCases[nearestNeighbors[j]-1];
				sumExpected += expectedCases[nearestNeighbors[j]-1];

				// Difference from original method:  stop when we see k cases, not k OTHER cases
				if( sumObserved >= k ){
					observedKValues[i] = sumObserved;
					break;
				}
			}

			// Watch off by 1 vector indexing in C as opposed to R
			observedMValues[i] = j + 1;

			// Difference from original method:  use actual observed k, not specified k
			x[0] = sumObserved - 1.0;
			observedPValues[i] = 1 - ppois(x, sumExpected)[0];
		}

		return List::create(
			_["observed.p.values"] = observedPValues,
			_["observed.m.values"] = observedMValues,
			_["observed.k.values"] = observedKValues
		);
	}
