/*
 *  SpatialEpi.h
 *
 *
 *  Created by Albert Y. Kim on 12/20/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */



/*
 ------------------------------------------------------------------------

 Utility Functions

 ------------------------------------------------------------------------
 */

/*
 getListElement:  Get list element named "str"
 */
static SEXP getListElement(SEXP list, char *str){

	SEXP elmt, names;
    PROTECT(elmt = R_NilValue);
	PROTECT(names = getAttrib(list, R_NamesSymbol));
    int i;

    for (i = 0; i < length(list); i++)
		if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
			elmt = VECTOR_ELT(list, i);
			break;
		}
	UNPROTECT(2);
    return elmt;
}


/*
 vecSum:  Computes the sum of the elements of a vector
 */
SEXP vecSum(SEXP Rvec){
	// SEXP coercion of R objects --------------------------
	PROTECT( Rvec = coerceVector(Rvec, REALSXP) );

	int i, n = length(Rvec);

	// Store result here + Initialize to first element
	SEXP value;
	PROTECT(value = allocVector(REALSXP,1) );
	REAL(value)[0] = 0.0;

	// Compute sum
	for (i = 0; i < n; i++)
		REAL(value)[0] = REAL(value)[0] + REAL(Rvec)[i];

	// Output results
	UNPROTECT(2);
	return value;
}


/*
 vecSum:  Computes the sum of the elements of a vector
 */
SEXP vecSumInt(SEXP Rvec){
	// SEXP coercion of R objects --------------------------
	PROTECT( Rvec = coerceVector(Rvec, INTSXP) );

	int i, n = length(Rvec);

	// Store result here + Initialize to first element
	SEXP value;
	PROTECT(value = allocVector(INTSXP,1) );
	INTEGER(value)[0] = 0;

	// Compute sum
	for (i = 0; i < n; i++)
		INTEGER(value)[0] = INTEGER(value)[0] + INTEGER(Rvec)[i];

	// Output results
	UNPROTECT(2);
	return value;
}


/*
 vecSumSq:  Computes the sum of the squares elements of a vector
 */
SEXP vecSumSq(SEXP Rvec){
	// SEXP coercion of R objects --------------------------
	PROTECT( Rvec = coerceVector(Rvec, REALSXP) );
	double *vec = REAL(Rvec);

	int i, n = length(Rvec);

	// Store result here + Initialize to first element
	SEXP value;
	PROTECT(value = allocVector(REALSXP,1) );
	REAL(value)[0] = 0.0;

	// Compute sum
	for (i = 0; i < n; i++)
		REAL(value)[0] = REAL(value)[0] + vec[i]*vec[i];

	// Output results
	UNPROTECT(2);
	return value;
}


/*
 max:	Computes the max of the elements of a vector
 */
SEXP max(SEXP Rvec){
	// SEXP coercion of R objects --------------------------
	PROTECT( Rvec = coerceVector(Rvec, REALSXP) );
	double *vec = REAL(Rvec);

	int i, n = length(Rvec);

	// Store result here + Initialize to first element
	SEXP maxValue;
	PROTECT(maxValue = allocVector(REALSXP,1) );
	REAL(maxValue)[0] = vec[0];

	// Compute Max
	for (i = 1; i < n; i++){
		if(REAL(maxValue)[0] < vec[i]){
			REAL(maxValue)[0] = vec[i];
		}
	}


	// Output results
	UNPROTECT(2);
	return maxValue;
}


/*
 whichMax:	Computes the index of the max of the vector
 If there are ties, it returns the smallest index
 */
SEXP whichMax(SEXP Rvec){
	// SEXP coercion of R objects --------------------------
	PROTECT( Rvec = coerceVector(Rvec, REALSXP) );
	double *vec = REAL(Rvec);

	int i, n = length(Rvec);

	// Store result here + Initialize to first element
	SEXP whichMaxValue, maxValue;
	PROTECT(whichMaxValue = allocVector(INTSXP,1) );
	PROTECT(maxValue = allocVector(REALSXP,1) );

	INTEGER(whichMaxValue)[0] = 0;
	REAL(maxValue)[0] = vec[0];

	// Compute Index of Max
	for (i = 1; i < n; i++){
		if(REAL(maxValue)[0] < vec[i]){
			REAL(maxValue)[0] = vec[i];
			INTEGER(whichMaxValue)[0] = i;
		}
	}

	// Output results
	UNPROTECT(3);
	return whichMaxValue;
}


/*
 min:	Computes the min of the elements of a vector
 */
SEXP min(SEXP Rvec){
	// SEXP coercion of R objects --------------------------
	PROTECT( Rvec = coerceVector(Rvec, REALSXP) );
	double *vec = REAL(Rvec);

	int i, n = length(Rvec);

	// Store result here + Initialize to first element
	SEXP minValue;
	PROTECT(minValue = allocVector(REALSXP,1) );
	REAL(minValue)[0] = vec[0];

	// Compute Min
	for (i = 1; i < n; i++){
		if(REAL(minValue)[0] > vec[i]){
			REAL(minValue)[0] = vec[i];
		}
	}

	// Output results
	UNPROTECT(2);
	return minValue;
}


/*
 whichMin:	Computes the index of the min of the vector
 If there are ties, it returns the smallest index
 */
SEXP whichMin(SEXP Rvec){
	// SEXP coercion of R objects --------------------------
	PROTECT( Rvec = coerceVector(Rvec, REALSXP) );
	double *vec = REAL(Rvec);

	int i,n = length(Rvec);

	// Store result here + Initialize to first element
	SEXP whichMinValue, minValue;
	PROTECT(whichMinValue = allocVector(INTSXP,1) );
	PROTECT(minValue = allocVector(REALSXP,1) );

	INTEGER(whichMinValue)[0] = 0;
	REAL(minValue)[0] = vec[0];

	// Compute Index of Min
	for (i = 1; i < n; i++){

		if(REAL(minValue)[0] > vec[i]){
			REAL(minValue)[0] = vec[i];
			INTEGER(whichMinValue)[0] = i;
		}
	}

	// Output results
	UNPROTECT(3);
	return whichMinValue;
}


/*
 compute log of density of multinomial
 */
SEXP ldmultinom(SEXP x_INPUT, SEXP n_INPUT, SEXP prob_INPUT){
	// SEXP coercion of R objects --------------------------
	PROTECT( x_INPUT = coerceVector(x_INPUT, REALSXP) );
	PROTECT( n_INPUT = coerceVector(n_INPUT, REALSXP) );
	PROTECT( prob_INPUT = coerceVector(prob_INPUT, REALSXP) );

	double n = REAL(n_INPUT)[0], *x = REAL(x_INPUT);
	int i;
	R_len_t len = length(x_INPUT);

	// Store result here + Initialize to first element
	SEXP p, value, sum_p;
	PROTECT(p = allocVector(REALSXP,len) );
	PROTECT(value = allocVector(REALSXP,1) );
	PROTECT(sum_p = allocVector(REALSXP,1) );

	sum_p = vecSum(prob_INPUT);

	// normalize probabilities if needed
	for(i=0; i<len; i++){
		REAL(p)[i] = REAL(prob_INPUT)[i]/REAL(sum_p)[0];
	}

	REAL(value)[0] = lgamma(n + 1);
	for(i=0; i<len; i++){
		if(REAL(p)[i]!=0){
			REAL(value)[0] -= lgamma(x[i] + 1);
			REAL(value)[0] += x[i]*log(REAL(p)[i]);
		}
	}


	// Output results
	UNPROTECT(6);
	return value;
}


/*
 *  Unequal Probability Sampling.
 *
 *  Modelled after Fortran code provided by:
 *    E. S. Venkatraman <venkat@biosta.mskcc.org>
 *  but with significant modifications in the
 *  "with replacement" case.
 */

/* Unequal probability sampling; with-replacement case */

SEXP ProbSampleReplace(SEXP p_input, SEXP nans){
	// SEXP coercion of R objects --------------------------
	PROTECT( p_input = coerceVector(p_input, REALSXP) );
	PROTECT( nans = coerceVector(nans, INTSXP) );

    double rU, sum;
    int i, j, n = length(p_input);
    int nm1 = n - 1;
	int perm[n];
	double p[n];

	SEXP ans;
	PROTECT( ans = allocVector(INTSXP, INTEGER(nans)[0] ) );

	sum = REAL(vecSum(p_input))[0];
	if( sum != 1 ){
		for(i=0; i<length(p_input); i++)
			REAL(p_input)[i] /= sum;
	}


    /* record element identities */
    for (i = 0; i < n; i++){
		perm[i] = i;
		p[i] = REAL(p_input)[i];
	}

    /* sort the probabilities into descending order */
    revsort(p, perm, n);


    /* compute cumulative probabilities */
    for (i = 1 ; i < n; i++){
		p[i] += p[i - 1];
	}

	GetRNGstate();
    /* compute the sample */
    for (i = 0; i < INTEGER(nans)[0]; i++) {
		rU = unif_rand();
		for (j = 0; j < nm1; j++) {
			if (rU <= p[j])
				break;
		}
		INTEGER(ans)[i] = perm[j];
    }
	PutRNGstate();


	UNPROTECT(3);
	return(ans);
}


/*
 Normalize probabilities over only p[include.indices]

 normalize <- function (p, include.indices=NA){
 p.prime <- p
 p.prime[!include.indices] <- 0
 p.prime <- p.prime/sum(p.prime)
 return( p.prime )
 }
 */



SEXP normalize(SEXP values, SEXP include_indices){
	PROTECT( values = coerceVector(values, REALSXP) );
	if(length(include_indices)!=0)
		PROTECT( include_indices = coerceVector(include_indices, INTSXP) );
	else
		PROTECT( include_indices = coerceVector(include_indices, NILSXP) );

	int i; double sum = 0, SXP = 2;
	R_len_t n = length(values);


	// Save Output Here
	SEXP p_prime;
	PROTECT( p_prime = allocVector(REALSXP, n) ); SXP++;

	if( length(include_indices)!=0 ){
		memset( REAL(p_prime), 0, n*sizeof(double) );
		for(i=0; i<n; i++)
			if( INTEGER(include_indices)[i] == 1 )
				REAL(p_prime)[i] = REAL(values)[i];
	}
	else{
		for(i=0; i<n; i++)
			REAL(p_prime)[i] = REAL(values)[i];
	}

	// Renormalize Probabilities
	sum = REAL(vecSum(p_prime))[0];
	for(i=0; i<n; i++)
		REAL(p_prime)[i] /= sum;

	UNPROTECT(3);
	return p_prime;
}



void sortVector(SEXP s, Rboolean decreasing)
{
    int n = LENGTH(s);
    if (n >= 2 && (decreasing || isUnsorted(s, FALSE)))
		switch (TYPEOF(s)) {
			case LGLSXP:
			case INTSXP:
				R_isort(INTEGER(s), n);
				break;
			case REALSXP:
				R_rsort(REAL(s), n);
				break;
		}
}



SEXP sort_matrix(SEXP input){
	PROTECT( input = coerceVector(input, INTSXP) );

	int n_row, n_col, i, j;

	SEXP Rdim;
	PROTECT( Rdim = getAttrib(input, R_DimSymbol));
	n_row = INTEGER(Rdim)[1];
	n_col = INTEGER(Rdim)[0];
	UNPROTECT(1);

	SEXP temp;
	PROTECT( temp = allocVector(INTSXP, n_col) );
	for(i=0; i<n_row; i++){
		for(j=0; j<n_col; j++){
			INTEGER(temp)[j] = INTEGER(input)[j+i*n_col];
		}
		sortVector(temp,0);
		for(j=0; j<n_col; j++){
			INTEGER(input)[j+i*n_col] = INTEGER(temp)[j];
		}
	}

	UNPROTECT(2);
	return input;
}


SEXP sort_matrix_double(SEXP input){
	PROTECT( input = coerceVector(input, REALSXP) );

	int n_row, n_col, i, j;

	SEXP Rdim;
	PROTECT( Rdim = getAttrib(input, R_DimSymbol));
	n_row = INTEGER(Rdim)[1];
	n_col = INTEGER(Rdim)[0];
	UNPROTECT(1);

	SEXP temp;
	PROTECT( temp = allocVector(REALSXP, n_col) );
	for(i=0; i<n_row; i++){
		for(j=0; j<n_col; j++){
			REAL(temp)[j] = REAL(input)[j+i*n_col];
		}
		sortVector(temp,0);
		for(j=0; j<n_col; j++){
			REAL(input)[j+i*n_col] = REAL(temp)[j];
		}
	}

	UNPROTECT(2);
	return input;
}




/*
 compute log of density of negative binomial allowing for non-integer counts
 */
SEXP ldnbinom(SEXP y, SEXP E, SEXP a, SEXP b){
	// SEXP coercion of R objects --------------------------
	PROTECT( y = coerceVector(y, REALSXP) );
	PROTECT( E = coerceVector(E, REALSXP) );
	PROTECT( a = coerceVector(a, REALSXP) );
	PROTECT( b = coerceVector(b, REALSXP) );

	// Store result here
	SEXP value;
	PROTECT(value = allocVector(REALSXP,1) );

	REAL(value)[0] = lgamma( REAL(y)[0]+REAL(a)[0] ) - lgamma( REAL(y)[0]+1 ) - lgamma( REAL(a)[0] ) +
	REAL(y)[0]*( log( REAL(E)[0] )-log( REAL(E)[0]+REAL(b)[0] ) ) +
	REAL(a)[0]*( log( REAL(b)[0] )-log( REAL(E)[0]+REAL(b)[0] ) );

	// Output results
	UNPROTECT(5);
	return value;
}






/*---------------------------------------------------------------------------*/
/*
 index <- match(
 apply(matrix1,1,function(x){paste(x,collapse="-")}),
 apply(matrix2,1,function(x){paste(x,collapse="-")})
 )
 */
SEXP match_matrices(SEXP matrix1,  SEXP matrix2){
	PROTECT( matrix1 = coerceVector(matrix1, INTSXP) );
	PROTECT( matrix2 = coerceVector(matrix2, INTSXP) );
	int SXP = 2;
	int n1, n2, k, i, j;

	// Get Matrix Dimensions
	SEXP Rdim1;
	PROTECT(Rdim1 = getAttrib(matrix1, R_DimSymbol));
	n1 = INTEGER(Rdim1)[1];  k = INTEGER(Rdim1)[0];
	SEXP Rdim2;
	PROTECT(Rdim2 = getAttrib(matrix2, R_DimSymbol));
	n2 = INTEGER(Rdim2)[1];
	UNPROTECT(2);

	// Convert Both Matrices Into Strings
	char value[1000], temp[1000];
	SEXP string1, string2;
	PROTECT( string1 = allocVector(STRSXP, n1) );
	PROTECT( string2 = allocVector(STRSXP, n2) );
	SXP = SXP + 2;

	for(i=0; i<n1; i++){
		strcpy(temp,"");
		strcpy(value,"");
		for(j=0; j<k; j++){
			sprintf(temp, "%i", INTEGER(matrix1)[i*k + j] );
			strcat(value,"&");
			strcat(value,temp);
		}
		SET_STRING_ELT(string1,i,mkChar(value));
	}

	for(i=0; i<n2; i++){
		strcpy(temp,"");
		strcpy(value,"");
		for(j=0; j<k; j++){
			sprintf(temp, "%i", INTEGER(matrix2)[i*k + j] );
			strcat(value,"&");
			strcat(value,temp);
		}
		SET_STRING_ELT(string2,i,mkChar(value));
	}

	SEXP matches;
	PROTECT( matches = match(string2,string1,0L) );
	SXP++;

	UNPROTECT(SXP);
	return matches;
}





/*
 config2vector <- function(config.input, n.zones, k){
 unit <- 10^ceiling(log10(n.zones))
 # when converting configuration matrices to vectors
 if(is.matrix(config.input)){
 config.output <- rep(0, nrow(config.input))
 for(i in 1:k)
 config.output <- config.output + config.input[,i] * unit^(k-i)
 }
 # when converting vectors to configuration matrices
 if(is.vector(config.input)){
 config.output <- matrix(0, nrow = length(config.input), ncol=k)
 for(i in 1:k)
 config.output[,i] <- trunc( (config.input %% unit^(k-i+1))/unit^(k-i) )
 }
 return(config.output)
 }
 ## theta <- sample.theta(1000, rep(TRUE,n.zones), 3, prop.post.1, overlap)$zones
 ## vector <- config2vector(theta,n.zones,k)
 ## sum((theta - config2vector(vector,n.zones,k))^2)
 ## sum((vector - config2vector(theta,n.zones,k))^2)
 */
SEXP config2vector(SEXP config_input, SEXP n_zones, SEXP k_input){
	/*
	 Main problem is that large int's don't seem to load right
	 */
	PROTECT( config_input = coerceVector(config_input, REALSXP) );
	PROTECT( n_zones = coerceVector(n_zones, INTSXP) );
	PROTECT( k_input = coerceVector(k_input, INTSXP) );
	int SXP = 3, k = INTEGER(k_input)[0];
	int i, j, n;
	double unit = pow(10, ceil( log10( (double) INTEGER(n_zones)[0] ) ));

	SEXP Rdim;
	PROTECT( Rdim = getAttrib(config_input, R_DimSymbol));
	SXP++;

	// if matrix, convert to vector
	if(length(Rdim)==2){
		n = INTEGER(Rdim)[1];
		SEXP config_output;
		PROTECT( config_output = allocVector(REALSXP,n) );
		SXP++;

		for(i=0; i<n; i++)
			REAL(config_output)[i] = 0;

		for(i=0; i<n; i++)
			for(j=0; j<k; j++){
				REAL(config_output)[i] += REAL(config_input)[i*k + j] * pow(unit, k-j-1) ;
			}

		UNPROTECT(SXP);
		return config_output;

	// if vector, convert to matrix
	}else{
		n = length(config_input);
		SEXP config_output;
		PROTECT( config_output = allocMatrix(REALSXP, k, n) ) ;
		SXP++;
		double value;


		for(i=0; i<n; i++)
			for(j=0; j<k; j++){
				// problems with size of int.  this simulates modular division
				value = REAL(config_input)[i] / pow(unit,k-j);
				value = value - floor(value);
				value = round(value * pow(unit,k-j));
				value = trunc( value / pow(unit,k-j-1) )  ;
				REAL(config_output)[i*k + j] = value;
			}
		UNPROTECT(SXP);
		return config_output;
	}


}






/*
SEXP unique_zones(SEXP input){
  PROTECT( input = coerceVector(input, INTSXP) );

  SEXP Rdim;
  PROTECT( Rdim = getAttrib(input, R_DimSymbol));
  int i, j, I = INTEGER(Rdim)[0], J = INTEGER(Rdim)[1];
  UNPROTECT(1);

  // find duplicated --------------------------
  char value[1000], temp[1000];

  int *h, n, count=0, counter;
  HashData data;
  SEXP string;
  PROTECT( string = allocVector(STRSXP, J) );

  for(j=0; j<J; j++){
    strcpy(temp,"");
    strcpy(value,"");
    // ignore first index, it is the parent index
    for(i=0; i<I; i++){
      sprintf(temp, "%i", INTEGER(input)[i + j*I] );
      strcat(value,"&");
      strcat(value,temp);
    }
    SET_STRING_ELT(string,j,mkChar(value));
  }

  HashTableSetup(string, &data);
  h = INTEGER(data.HashTable);
    PROTECT(data.HashTable);

  SEXP ans;
  PROTECT(ans = allocVector(LGLSXP, J));

    for(i = 0; i < data.M; i++){
    h[i] = NIL;
  }
  for(i = 0; i < J; i++)
    LOGICAL(ans)[i] = isDuplicated(string, i, &data);

  // return unique --------------------------
  for(i=0; i<J; i++){
    if(LOGICAL(ans)[i] == 0)
      count++;
  }

  SEXP output;
  PROTECT(output = allocMatrix(INTSXP, I, count));
  counter = 0;
  for(i=0; i<J; i++){
    if(LOGICAL(ans)[i] == 0){
      for(j=0; j<I; j++){
        INTEGER(output)[counter*I + j] = INTEGER(input)[i*I + j];
      }
      counter++;
    }
  }


  UNPROTECT(5);
    return output;
}
*/
