#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Random.h>
#include "SpatialEpi.h"
#include "cluster_detection.h"

SEXP return_global_moves(SEXP theta, SEXP overlap, SEXP cluster_coords){
	PROTECT( theta = coerceVector(theta, INTSXP) );
	if(length(overlap)==2)
		PROTECT( overlap = coerceVector(overlap, VECSXP) );
	else
		PROTECT( overlap = coerceVector(overlap, INTSXP) );
	PROTECT( cluster_coords = coerceVector(cluster_coords, INTSXP) );

	int SXP = 3, k = length(theta), n_zones = length(cluster_coords)/2,
	i, j, z, l, m, zone, zone_sub,
	num_birth_moves,
	counter;



	// Check if using list version of overlap (slower) or matrix (look-up table i.e. quicker)
	int large_flag = 0, overlap_vector[n_zones];
	SEXP presence, cluster_list;
	if(length(overlap)==2){
		large_flag = 1;
		PROTECT(presence = getListElement(overlap, "presence"));
		PROTECT(cluster_list = getListElement(overlap, "cluster_list"));
		SXP += 2;
	}

	// Save Moves Here
	SEXP death_moves, birth_moves, death_output, birth_output, move_possible;

	PROTECT(move_possible = allocVector(INTSXP, 2));
	SXP++;


	//---------------------------------------------------------------------
	// Death Moves
	//---------------------------------------------------------------------
	switch(k){
		case 0:
			PROTECT(death_output = R_NilValue);
			SXP++;
			INTEGER(move_possible)[0] = 0;
			break;
		case 1:
			PROTECT(death_output = R_NilValue);
			SXP++;
			INTEGER(move_possible)[0] = 1;
			break;
		default:
			PROTECT(death_moves = allocMatrix(INTSXP, k-1, k));
			SXP++;
			for(j=0; j<k; j++){
				counter=0;
				for(i=0; i<k; i++)
					if(i!=j){
						INTEGER(death_moves)[j*(k-1) + counter] = INTEGER(theta)[i];
						counter++;
					}
			}

			PROTECT(death_output = death_moves);
			SXP++;
			INTEGER(move_possible)[0] = 1;
			break;
	}

	//---------------------------------------------------------------------
	// Birth Moves
	//---------------------------------------------------------------------
	if(k==0){
		PROTECT(birth_output = allocMatrix(INTSXP, 1, n_zones));
		SXP++;
		for(z=0; z<n_zones; z++)
			INTEGER(birth_output)[z] = z + 1;  // Watch R vs C Indexing

		INTEGER(move_possible)[1] = 1;

	}else{
		//---------------------------------------------------------------------
		// Check overlap with existing single zones
		//---------------------------------------------------------------------
		PROTECT(birth_moves = allocVector(INTSXP, n_zones));
		SXP++;

		for(z=0; z<n_zones; z++)
			INTEGER(birth_moves)[z] = 1;

		for(j=0; j<k; j++){
			zone = INTEGER(theta)[j] - 1;  // Watch R vs C Indexing

			//-----------------------------------------------------------------
			// Start split on overlap type:  matrix vs list
			//-----------------------------------------------------------------
			if(large_flag==0){
				for(z=0; z<n_zones; z++)
					INTEGER(birth_moves)[z] *= 1-INTEGER(overlap)[n_zones*zone + z];
			}else{
				memset( overlap_vector, 0, n_zones*sizeof(int) );
				for(m=0; m<length(VECTOR_ELT(cluster_list, zone)); m++){
					zone_sub = INTEGER( VECTOR_ELT(cluster_list, zone) )[m] - 1;
					for(l=0; l<length(VECTOR_ELT(presence, zone_sub)); l++)
						overlap_vector[ INTEGER(VECTOR_ELT(presence, zone_sub))[l] - 1 ] = 1;
				}
				for(z=0; z<n_zones; z++)
					INTEGER(birth_moves)[z] *= 1-overlap_vector[z];
			}
			// End split on overlap type:  matrix vs list
		} // end of j=0..k; checking overlap over all single zones



		//----------------------------------------------------------
		// Make reduced table of all possible configs
		//----------------------------------------------------------
		num_birth_moves = INTEGER(vecSumInt(birth_moves))[0];
		if(num_birth_moves==0){
			PROTECT(birth_output = R_NilValue);
			SXP++;
			INTEGER(move_possible)[1] = 0;
		}else{
			PROTECT(birth_output = allocMatrix(INTSXP, k+1, num_birth_moves));
			SXP++;

			// Fill in original theta
			for(i=0; i<num_birth_moves; i++)
				for(j=0; j<k; j++)
					INTEGER(birth_output)[i*(k+1)+j] = INTEGER(theta)[j];

			// Fill in new single zone
			counter=0;
			for(z=0; z<n_zones; z++){
				if( INTEGER(birth_moves)[z] == 1 ){
					INTEGER(birth_output)[counter*(k+1)+k] = z + 1;  // Watch R vs C Indexing
					counter++;
				}
			}

			birth_output = sort_matrix(birth_output);
			INTEGER(move_possible)[1] = 1;

		}


	} // end of if/else k==0


	//----------------------------------------------------------
	// Store in one list
	//----------------------------------------------------------
	SEXP moves, moves_names;
	PROTECT(moves = allocVector(VECSXP, 2)); SXP++;

	PROTECT(moves_names = allocVector(STRSXP, 2));
	SET_STRING_ELT(moves_names, 0,  mkChar("death"));
	SET_STRING_ELT(moves_names, 1,  mkChar("birth"));
	setAttrib(moves, R_NamesSymbol, moves_names);
	UNPROTECT(1);

	SET_VECTOR_ELT(moves, 0, death_output);
	SET_VECTOR_ELT(moves, 1, birth_output);




	//-------------------------------------------------------------------------
	// Output Results As R List
	//-------------------------------------------------------------------------
	SEXP results, results_names;
	PROTECT(results = allocVector(VECSXP, 2)); SXP++;

	PROTECT(results_names = allocVector(STRSXP, 2));
	SET_STRING_ELT(results_names, 0,  mkChar("moves"));
	SET_STRING_ELT(results_names, 1,  mkChar("move_possible"));
	setAttrib(results, R_NamesSymbol, results_names);
	UNPROTECT(1);

	SET_VECTOR_ELT(results, 0, moves);
	SET_VECTOR_ELT(results, 1, move_possible);

	UNPROTECT(SXP);
	return results;
}





SEXP return_local_moves(SEXP theta, SEXP overlap, SEXP cluster_coords, SEXP move_input){
	PROTECT(theta = coerceVector(theta, INTSXP));
	if(length(overlap)==2)
		PROTECT(overlap = coerceVector(overlap, VECSXP));
	else
		PROTECT(overlap = coerceVector(overlap, INTSXP));
	PROTECT(cluster_coords = coerceVector(cluster_coords, INTSXP));
	PROTECT(move_input = coerceVector(move_input, INTSXP));


	int SXP = 4, k = length(theta), n_zones = length(cluster_coords)/2,
	same_center[n_zones], diff_center[n_zones],
	move = INTEGER(move_input)[0],  // Watch R vs C Indexing: Internal, so we don't adjust by 1
	i, j, z, l, m, zone, zone_sub, zone_sub2, start, end,
	num_moves,
	counter;


	// Check if using list version of overlap (slower) or matrix (look-up table i.e. quicker)
	int large_flag = 0, overlap_vector[n_zones];
	SEXP presence, cluster_list;
	if(length(overlap)==2){
		large_flag = 1;
		PROTECT(presence = getListElement(overlap,"presence"));
		PROTECT(cluster_list = getListElement(overlap,"cluster_list"));
		SXP += 2;
	}

	// Save Moves Here
	SEXP moves, moves_reduced, move_possible;



	//-------------------------------------------------------------------------
	//
	// Get all growth/trim/replacement (local) moves.  The procedure has similar
	// structure for all three moves
	//
	//-------------------------------------------------------------------------
	PROTECT(moves = allocMatrix(INTSXP, n_zones, k));
	SXP++;
	memset(INTEGER(moves), 0, n_zones*k*sizeof(int));

	for(j=0; j<k; j++){
		zone = INTEGER(theta)[j] - 1; // Watch R vs C Indexing
		//---------------------------------------------------------------------
		// Identify indices of smallest (start) & largest (end) unit zones
		// with same center.  i.e. [start, end) zones have same center
		//---------------------------------------------------------------------
		start = zone;
		end = zone + 1;
		while(INTEGER(cluster_coords)[start-1] == INTEGER(cluster_coords)[zone] && start != 0)
			start--;
		while(INTEGER(cluster_coords)[end] == INTEGER(cluster_coords)[zone] && end != n_zones)
			end++;


		//---------------------------------------------------------------------
		// Identify single zones with same center/different center
		//---------------------------------------------------------------------
		for(z=0; z<n_zones; z++)
			diff_center[z] = 1;
		memset( same_center, 0, n_zones*sizeof(int) );

		for(z=start; z<end; z++){
			same_center[z] = 1;
			diff_center[z] = 0;
		}


		//---------------------------------------------------------------------
		// Check overlap with all OTHER single zones i.e. i!=j
		//---------------------------------------------------------------------
		for(i=0; i<k; i++)
			if(i != j){
				zone_sub = INTEGER(theta)[i] - 1; // Watch R vs C Indexing

				// Start split on overlap type:  matrix vs list
				if(large_flag==0){
					// Same center
					for(z=start; z<end; z++)
						same_center[z] *= 1-INTEGER(overlap)[n_zones*zone_sub + z];
					// Different center
					for(z=0; z<n_zones; z++)
						diff_center[z] *= 1-INTEGER(overlap)[n_zones*zone_sub + z];

				}else{
					memset( overlap_vector, 0, n_zones*sizeof(int) );
					for(m=0; m<length(VECTOR_ELT(cluster_list, zone_sub)); m++){
						zone_sub2 = INTEGER( VECTOR_ELT(cluster_list, zone_sub) )[m] - 1;
						for(l=0; l<length(VECTOR_ELT(presence,zone_sub2)); l++)
							overlap_vector[ INTEGER(VECTOR_ELT(presence,zone_sub2))[l] - 1 ] = 1;
					}

					// Same center
					for(z=start; z<end; z++)
						same_center[z] *= 1-overlap_vector[z];
					// Different center
					for(z=0; z<n_zones; z++)
						diff_center[z] *= 1-overlap_vector[z];
				}
				// End split on overlap type:  matrix vs list
			}


		//---------------------------------------------------------------------
		// Store moves
		//---------------------------------------------------------------------
		switch(move){
			case 0:
				// Growth
				for(z=zone+1; z<end; z++){
					INTEGER(moves)[n_zones*j + z] = same_center[z];
				}
				break;
			case 1:
				// Trim
				for(z=start; z<zone; z++){
					INTEGER(moves)[n_zones*j + z] = same_center[z];
				}
				break;
			case 2:
				// Replace
				for(z=0; z<n_zones; z++){
					INTEGER(moves)[n_zones*j + z] = diff_center[z];
				}
				break;
		}
	} // end of j = 1:k; over all component single zones




	//-------------------------------------------------------------------------
	//
	// Make reduced table of all possible configs.  Same for all local moves
	//
	//-------------------------------------------------------------------------
	PROTECT(move_possible = allocVector(INTSXP, 1));
	SXP++;

	num_moves = INTEGER(vecSumInt(moves))[0];
	if(num_moves==0){
		PROTECT(moves_reduced = R_NilValue);
		SXP++;
		INTEGER(move_possible)[0] = 0;
	}else{
		PROTECT(moves_reduced = allocMatrix(INTSXP, k, num_moves));
		SXP++;
		// Fill in old theta
		for(i=0; i<num_moves; i++)
			for(j=0; j<k; j++)
				INTEGER(moves_reduced)[i*k + j] = INTEGER(theta)[j];

		// Adjust configs to include new single zone
		counter=0;
		for(j=0; j<k; j++)
			for(z=0; z<n_zones; z++)
				if(INTEGER(moves)[n_zones*j + z] == 1){
					INTEGER(moves_reduced)[counter*k + j] = z + 1; // Watch R vs C Indexing
					counter++;
				}
		moves_reduced = sort_matrix(moves_reduced);
		INTEGER(move_possible)[0] = 1;
	}


	//-------------------------------------------------------------------------
	// Output Results As R List
	//-------------------------------------------------------------------------
	SEXP results, results_names;
	PROTECT(results = allocVector(VECSXP, 2)); SXP++;

	PROTECT(results_names = allocVector(STRSXP, 2));
	SET_STRING_ELT(results_names, 0,  mkChar("moves"));
	SET_STRING_ELT(results_names, 1,  mkChar("move_possible"));
	setAttrib(results, R_NamesSymbol, results_names);
	UNPROTECT(1);

	SET_VECTOR_ELT(results, 0, moves_reduced);
	SET_VECTOR_ELT(results, 1, move_possible);

	UNPROTECT(SXP);
	return results;
}





SEXP MCMC(SEXP n_sim, SEXP pattern, SEXP theta_init, SEXP overlap, SEXP cluster_coords, SEXP p_moves_orig,
		  SEXP k_max, SEXP lkhd_z, SEXP lambda){
	PROTECT( n_sim = coerceVector(n_sim, INTSXP) );
	PROTECT( pattern = coerceVector(pattern, INTSXP) );
	PROTECT( theta_init = coerceVector(theta_init, INTSXP) );
	// Check If List or Not
	if(length(overlap)==2)
		PROTECT( overlap = coerceVector(overlap, VECSXP) );
	else
		PROTECT( overlap = coerceVector(overlap, INTSXP) );
	PROTECT( cluster_coords = coerceVector(cluster_coords, INTSXP) );
	PROTECT( p_moves_orig = coerceVector(p_moves_orig, REALSXP) );
	PROTECT( k_max = coerceVector(k_max, INTSXP) );
	PROTECT( lkhd_z = coerceVector(lkhd_z, REALSXP) );
	PROTECT( lambda = coerceVector(lambda, REALSXP) );

	int SXP=9, SXP2=0, n = INTEGER(n_sim)[0], n_zones = length(cluster_coords)/2,
	n_moves, n_moves_rev, k, k_rev, K = INTEGER(k_max)[0],
	no_moves_flag, match_flag, reverse[5] = {1,0,2,4,3}, z_star,
	sim, i, j, z, sampling_type, index,
	move_rev_index;
	double p_num, p_denom, sum, u;


	//-------------------------------------------------------------------------
	//
	// Set up SEXP memory objects
	//
	//-------------------------------------------------------------------------
	SEXP theta, theta_star, Rdim, n_samp;
	PROTECT(n_samp = allocVector(INTSXP,1));
	INTEGER(n_samp)[0] = 1;
	SXP++;

	//----------------------------------------------------------
	// Move related
	//----------------------------------------------------------
	SEXP move, move_rev,
	p_moves, p_moves_rev,
	pos_local_moves, pos_local_moves_rev,
	pos_global_moves, pos_global_moves_rev,
	moves, moves_rev;

	PROTECT(p_moves = allocVector(REALSXP, length(p_moves_orig)));
	PROTECT(p_moves_rev = allocVector(REALSXP, length(p_moves_orig)));
	PROTECT(move = allocVector(INTSXP,1));
	PROTECT(move_rev = allocVector(INTSXP,1));
	SXP += 4;

	//----------------------------------------------------------
	// Sampling
	//----------------------------------------------------------
	// Either Uniform or Proportional Sampling
	SEXP unif_p, prop_p, p, probs;
	PROTECT(unif_p = allocVector(REALSXP, n_zones));
	PROTECT(prop_p = allocVector(REALSXP, n_zones));
	PROTECT(p = allocVector(REALSXP, n_zones));
	SXP += 3;

	sum = REAL(vecSum(lkhd_z))[0];
	for(z=0; z<n_zones; z++){
		REAL(unif_p)[z] = (double) 1/n_zones;
		REAL(prop_p)[z] = REAL(lkhd_z)[z]/sum;
	}


	//----------------------------------------------------------
	// Store Outputs Here
	//----------------------------------------------------------
	SEXP sample, move_trace, accpt_trace, ratio;

	PROTECT(sample = allocVector(VECSXP,n));
	PROTECT(move_trace = allocVector(INTSXP,n));
	PROTECT(accpt_trace = allocVector(INTSXP,n));
	PROTECT(ratio = allocVector(REALSXP,n));
	SXP += 4;

	memset( INTEGER(move_trace), 0, n*sizeof(int) );
	memset( INTEGER(accpt_trace), 0, n*sizeof(int) );


	//-------------------------------------------------------------------------
	//
	// Main MCMC Algorithm
	//
	//-------------------------------------------------------------------------
	// First Value
	SET_VECTOR_ELT(sample, 0, theta_init);

	// For loop: over all simulations
	for(sim=1; sim<n; sim++){
		//---------------------------------------------------------------------
		//
		// Intialize stuff
		//
		//---------------------------------------------------------------------
		//------------------------------------------------------
		// Select Sampling Type:  Remove Later
		//------------------------------------------------------
		sampling_type = INTEGER( pattern )[ sim % length(pattern) ];
		if( sampling_type == 0 ){
			for(z=0; z<n_zones; z++)
				REAL(p)[z] = REAL(unif_p)[z];
		}else{
			for(z=0; z<n_zones; z++)
				REAL(p)[z] = REAL(prop_p)[z];
		}

		//------------------------------------------------------
		// Reset Values
		//------------------------------------------------------
		PROTECT(theta = VECTOR_ELT(sample, sim-1));
		SXP2 = 1;
		k = length(theta);
		p_num = 1;
		p_denom = 1;
		// Probs of each move
		for(i=0; i<length(p_moves_orig); i++){
			REAL(p_moves)[i] = REAL(p_moves_orig)[i];
			REAL(p_moves_rev)[i] = REAL(p_moves_orig)[i];
		}
		no_moves_flag = 0;




		//---------------------------------------------------------------------
		//
		// Forward moves
		//
		//---------------------------------------------------------------------
		//------------------------------------------------------
		// Figure out legal forward moves and get their probs
		//------------------------------------------------------
		//---------------------------------------
		// k==0: only birth moves
		//---------------------------------------
		if(k==0){
			memset(REAL(p_moves), 0, length(p_moves_orig)*sizeof(double));
			REAL(p_moves)[4] = 1;
		}

		//---------------------------------------
		// k==K: no birth moves
		//---------------------------------------
		if(k==K){
			REAL(p_moves)[4] = 0;
		}

		//------------------------------------------------------
		// Renormalize prob and sample move
		//------------------------------------------------------
		sum = REAL(vecSum(p_moves))[0];
		for(i=0; i<length(p_moves); i++)
			REAL(p_moves)[i] /= sum;




		INTEGER(move)[0] = INTEGER(ProbSampleReplace(p_moves, n_samp))[0];
		p_num *= REAL(p_moves)[INTEGER(move)[0]];



		//------------------------------------------------------
		// Set up sampling space of configurations
		//------------------------------------------------------
		//---------------------------------------
		// Local moves
		//---------------------------------------
		if(INTEGER(move)[0] == 0 || INTEGER(move)[0] == 1 || INTEGER(move)[0] == 2){
			PROTECT(pos_local_moves = return_local_moves(theta, overlap, cluster_coords, move));
			SXP2++;

			if( INTEGER(getListElement(pos_local_moves,"move_possible"))[0] == 0 ){
				no_moves_flag = 1;
				PROTECT(moves = R_NilValue);
				SXP2++;
			}else{
				PROTECT(moves = getListElement(pos_local_moves,"moves"));
				SXP2++;
			}
		} // end if local move: growth, trim, replace


		//---------------------------------------
		// Global moves
		//---------------------------------------
		if(INTEGER(move)[0] == 3 || INTEGER(move)[0] == 4){
			PROTECT(pos_global_moves = return_global_moves(theta, overlap, cluster_coords));
			SXP2++;

			// Change coding of global moves
			if( INTEGER(getListElement(pos_global_moves,"move_possible"))[ INTEGER(move)[0]-3 ] == 0 ){
				no_moves_flag = 1;
				PROTECT(moves = R_NilValue);
				SXP2++;
			}else{
				PROTECT(moves = VECTOR_ELT(getListElement(pos_global_moves,"moves"),INTEGER(move)[0]-3 ));
				SXP2++;
			}
		} // end if global move: death, birth



		//---------------------------------------------------------------------
		//
		// Continue only if move exists.  If doesn't, set ratio = 0 i.e. always reject
		//
		//---------------------------------------------------------------------
		if(no_moves_flag != 1){
			//--------------------------------------------------
			// Sample configs: handle null configuration differently
			// i.e. if k==1 and death move
			//--------------------------------------------------
			if(INTEGER(move)[0]==3 && k==1){
				PROTECT( theta_star = R_NilValue );
				SXP2++;
				p_num *= 1;

			}else{
				// Get number of possible forward proposal configurations and k*
				PROTECT( Rdim = getAttrib(moves, R_DimSymbol));
				k_rev = INTEGER(Rdim)[0];
				n_moves = INTEGER(Rdim)[1];
				UNPROTECT(1);


				// Compute probability of each possible forward proposed configuration
				PROTECT(probs = allocVector(REALSXP, n_moves));
				for(i=0; i<n_moves; i++){
					REAL(probs)[i] = 1.0;
					for(j=0; j<k_rev; j++){
						index = INTEGER(moves)[i*k_rev + j] - 1; // Watch R vs C Indexing
						REAL(probs)[i] *= REAL(p)[ index ];
					}
				}

				//Rprintf("Forward probs:\n");
				sum = REAL(vecSum(probs))[0];
				for(i=0; i<n_moves; i++){
					REAL(probs)[i] /= sum;
					//Rprintf("%4.6f, ", REAL(probs)[i]);
				}
				//Rprintf("\n");


				// Sample a new configuration theta_star and sort it
				z_star = INTEGER(ProbSampleReplace(probs, n_samp))[0];
				p_num *= REAL(probs)[z_star];
				UNPROTECT(1);


				PROTECT( theta_star = allocVector(INTSXP, k_rev) );
				SXP2++;
				for(j=0; j<k_rev; j++){
					INTEGER(theta_star)[j] = INTEGER(moves)[z_star*k_rev + j]; // Watch R vs C Indexing
				}

				if(length(theta_star)!=0)
					sortVector(theta_star,0);
			} // end if/else null configuration





			//-----------------------------------------------------------------
			//
			// Backward moves
			//
			//-----------------------------------------------------------------
			k_rev = length(theta_star);
			INTEGER(move_rev)[0] = reverse[INTEGER(move)[0]];


			//------------------------------------------------------
			// Figure out legal backward moves and get their probs
			//------------------------------------------------------
			//-----------------------------------
			// k_rev==0: only birth moves
			//-----------------------------------
			if(k_rev==0){
				memset(REAL(p_moves_rev), 0, length(p_moves_orig)*sizeof(double));
				REAL(p_moves_rev)[4] = 1;
			}

			//-----------------------------------
			// k_rev==K: no death moves
			//-----------------------------------
			if(k_rev==K){
				REAL(p_moves_rev)[4] = 0;
			}

			//--------------------------------------------------
			// Renormalize prob and sample move
			//--------------------------------------------------
			sum = REAL(vecSum(p_moves_rev))[0];
			for(i=0; i<length(p_moves_rev); i++)
				REAL(p_moves_rev)[i] /= sum;

			p_denom *= REAL(p_moves_rev)[INTEGER(move_rev)[0]];


			//--------------------------------------------------
			// Set up sampling space of configurations
			//--------------------------------------------------
			//-----------------------------------
			// Local moves
			//-----------------------------------
			if(INTEGER(move_rev)[0] == 0 || INTEGER(move_rev)[0] == 1 || INTEGER(move_rev)[0] == 2){
				PROTECT(pos_local_moves_rev = return_local_moves(theta_star, overlap, cluster_coords, move_rev));
				SXP2++;

				if( INTEGER(getListElement(pos_local_moves_rev,"move_possible"))[0] == 0 ){
					// superfluous, as existence of theta guarantees move exists
				}else{
					PROTECT(moves_rev = getListElement(pos_local_moves_rev,"moves"));
					SXP2++;
				}
			} // end if local move: growth, trim, replace


			//-----------------------------------
			// Global moves
			//-----------------------------------
			if( INTEGER(move_rev)[0] == 3 || INTEGER(move_rev)[0] == 4 ){
				PROTECT(pos_global_moves_rev = return_global_moves(theta_star, overlap, cluster_coords));
				SXP2++;

				move_rev_index = INTEGER(move_rev)[0]-3;

				if( INTEGER(getListElement(pos_global_moves_rev,"move_possible"))[ move_rev_index ] == 0 ){
				}else{
					PROTECT(moves_rev = VECTOR_ELT(getListElement(pos_global_moves_rev,"moves"), move_rev_index));
					SXP2++;
				}
			} // end if global move: death, birth




			//Rprintf("Theta Star:\n");
			//for(i=0; i<length(theta_star); i++)
			//	Rprintf("%i ", INTEGER(theta_star)[i]);
			//Rprintf("\n");
			//--------------------------------------------------
			// Sample configs: handle null configuration differently
			// i.e. if k_rev==1 and death is reverse move
			//--------------------------------------------------
			if(INTEGER(move_rev)[0]==3 && k_rev==1){
				p_denom *= 1;
			}else{
				// Get number of possible backward proposal configurations
				PROTECT( Rdim = getAttrib(moves_rev, R_DimSymbol));
				n_moves_rev = INTEGER(Rdim)[1];
				UNPROTECT(1);

				// Compute probability of each possible forward proposed configuration
				PROTECT(probs = allocVector(REALSXP, n_moves_rev));
				for(i=0; i<n_moves_rev; i++){
					REAL(probs)[i] = 1.0;
					for(j=0; j<k; j++){
						z = INTEGER(moves_rev)[i*k + j] - 1; // Watch R vs C Indexing
						REAL(probs)[i] *= REAL(p)[z];
					}
				}


				//Rprintf("Reverse probs:\n");
				sum = REAL(vecSum(probs))[0];
				for(i=0; i<n_moves_rev; i++){
					REAL(probs)[i] /= sum;
					//Rprintf("%4.6f, ", REAL(probs)[i]);
				}
				//Rprintf("\n");

				// match index
				for(i=0; i<n_moves_rev; i++){
					match_flag = 1;
					for(j=0; j<k; j++){
						if( INTEGER(moves_rev)[i*k + j] == INTEGER(theta)[j] ){
							match_flag *= 1;
						}else{
							match_flag *= 0;
							break;
						}
					}
					if(match_flag == 1)
						break;
				}
				z = i;
				p_denom *= REAL(probs)[z];
				UNPROTECT(1);
				// Rprintf("%i %i %4.6f %4.6f\n", z_star+1, z+1, p_num, p_denom);

			} // end if/else null configuration

			//Rprintf("Num and Denom:\n");
			REAL(ratio)[sim-1] = p_denom/p_num;

			REAL(ratio)[sim-1] *= REAL(lambda)[k_rev]; // Watch R vs C Indexing
			if( k_rev != 0 ){
				for(j=0; j<k_rev; j++)
					REAL(ratio)[sim-1] *= REAL(lkhd_z)[ INTEGER(theta_star)[j]-1 ]; // Watch R vs C Indexing
			}

			REAL(ratio)[sim-1] /= REAL(lambda)[k]; // Watch R vs C Indexing
			if( k != 0 ){
				for(j=0; j<k;j++)
					REAL(ratio)[sim-1] /= REAL(lkhd_z)[ INTEGER(theta)[j]-1 ]; // Watch R vs C Indexing
			}

		}else{
			REAL(ratio)[sim-1] = 0;
		} // end if/else move exists


		//-------------------------------------------------------------------------
		//
		// Accept/Reject:  As ratio goes up, more chance of accepting
		//
		//-------------------------------------------------------------------------
		GetRNGstate();
		u = unif_rand();
		PutRNGstate();

		if( u < REAL(ratio)[sim-1] ){
			SET_VECTOR_ELT(sample, sim, theta_star);
			INTEGER(accpt_trace)[sim-1] = 1;
		}else{
			SET_VECTOR_ELT(sample, sim, theta);
		}

		INTEGER(move_trace)[sim-1] = INTEGER(move)[0] + 1; // Watch R vs C Indexing
		UNPROTECT(SXP2);
	} // end for loop over all simulations


	//-------------------------------------------------------------------------
	//
	// Output Results As R List
	//
	//-------------------------------------------------------------------------
	SEXP output, output_names;
	PROTECT(output = allocVector(VECSXP, 4)); SXP++;

	PROTECT(output_names = allocVector(STRSXP, 4));
	SET_STRING_ELT(output_names, 0,  mkChar("sample"));
	SET_STRING_ELT(output_names, 1,  mkChar("move_trace"));
	SET_STRING_ELT(output_names, 2,  mkChar("accpt_trace"));
	SET_STRING_ELT(output_names, 3,  mkChar("ratio_trace"));
	setAttrib(output, R_NamesSymbol, output_names);
	UNPROTECT(1);

	SET_VECTOR_ELT(output, 0, sample);
	SET_VECTOR_ELT(output, 1, move_trace);
	SET_VECTOR_ELT(output, 2, accpt_trace);
	SET_VECTOR_ELT(output, 3, ratio);

	UNPROTECT(SXP);
	return output;
}



/*
 Computes the internally standardized expected number of counts

 IMPORTANT:  We assume that the population and cases vectors are balanced:
 all counts are sorted by area first, and then within each area the
 counts for all strata are listed (even if 0 count) in the same order.
 That way we can loop over indices easily.

 INPUT:
 -population_INPUT:  population for each strata in each area
 -cases_INPUT:  cases for each strata in each ares
 -nStrata_INPUT:  number of strata considered

 OUTPUT:
 -E:  expected counts
 -q:  strata-specific rates
 */
SEXP expected( SEXP population_INPUT, SEXP cases_INPUT, SEXP nStrata_INPUT ){
	// SEXP coercion of R objects --------------------------
	PROTECT( population_INPUT = coerceVector(population_INPUT, REALSXP) );
	PROTECT( cases_INPUT = coerceVector(cases_INPUT, REALSXP) );
	PROTECT( nStrata_INPUT = coerceVector(nStrata_INPUT, INTSXP) );

	double	*population = REAL(population_INPUT),
	*cases = REAL(cases_INPUT);
	int		nStrata = INTEGER(nStrata_INPUT)[0],
	n = length(population_INPUT)/nStrata;
	int i, j;


	// Creation of SEXP objects --------------------------
	SEXP qNum, qDenom, q, E;
	PROTECT(q = allocVector(REALSXP,nStrata) );
	PROTECT(E = allocVector(REALSXP,n) );
	PROTECT(qNum = allocVector(REALSXP,nStrata) );
	PROTECT(qDenom = allocVector(REALSXP,nStrata) );


	// Initialize vectors to 0
	for(i = 0; i < nStrata; i++){
		REAL(qNum)[i] = 0.0; REAL(qDenom)[i] = 0.0; REAL(q)[i]=0.0;
	}
	for(i = 0; i < n; i++){
		REAL(E)[i] = 0.0;
	}


	// Compute q:  strata-specific rates.  We do the numerator and denominator separately
	for(i = 0; i < nStrata; i++){
		for(j = 0; j < n; j++){
			REAL(qNum)[i] = REAL(qNum)[i] +  cases[j*nStrata + i];
			REAL(qDenom)[i] = REAL(qDenom)[i] +  population[j*nStrata + i];
		}
	}
	for(i = 0; i < nStrata; i++){
		REAL(q)[i] = REAL(qNum)[i]/REAL(qDenom)[i];
	}
	// unprotect num and denominator
	UNPROTECT(2);


	// Compute E expected counts
	for(i=0; i < n; i++)
		for(j=0; j < nStrata; j++){
			REAL(E)[i] = REAL(E)[i] + REAL(q)[j] * population[i*nStrata + j];
		}


	// output results as R list
	SEXP list, list_names;
	PROTECT(list_names = allocVector(STRSXP, 2));
	SET_STRING_ELT(list_names, 0,  mkChar("E"));
	SET_STRING_ELT(list_names, 1,  mkChar("q"));

	PROTECT(list = allocVector(VECSXP, 2));	// Creating list w/2 vector elements
	SET_VECTOR_ELT(list, 0, E);				// attaching vector to list
	SET_VECTOR_ELT(list, 1, q);

	setAttrib(list, R_NamesSymbol, list_names); //and attaching the vector names


	// Output results
	UNPROTECT(7);
	return list;
}



SEXP coeff(SEXP y_INPUT, SEXP E_INPUT, SEXP a_INPUT, SEXP b_INPUT, SEXP cluster_list_INPUT){
	// SEXP coercion of R objects --------------------------
	PROTECT( y_INPUT = coerceVector(y_INPUT, REALSXP) );
	PROTECT( E_INPUT = coerceVector(E_INPUT, REALSXP) );
	PROTECT( a_INPUT = coerceVector(a_INPUT, REALSXP) );
	PROTECT( b_INPUT = coerceVector(b_INPUT, REALSXP) );
	PROTECT( cluster_list_INPUT = coerceVector(cluster_list_INPUT, VECSXP) );
	int SXP = 5;

	int i,j,k,index;

	R_len_t n = length(E_INPUT), nZones = length(cluster_list_INPUT);
	// double *y = REAL(y_INPUT), *E = REAL(E_INPUT), *a = REAL(a_INPUT), *b = REAL(b_INPUT), Ez, yz;
	SEXP y, E, a, b, yz, Ez;
	PROTECT( y = allocVector(REALSXP,1) );
	PROTECT( E = allocVector(REALSXP,1) );
	PROTECT( a = allocVector(REALSXP,1) );
	PROTECT( b = allocVector(REALSXP,1) );
	PROTECT( yz = allocVector(REALSXP,1) );
	PROTECT( Ez = allocVector(REALSXP,1) );
	SXP = SXP + 6;

	SEXP log_fy_n, sum_log_fy_inside, coeff, y_temp, E_temp, inside;
	PROTECT( log_fy_n = allocVector(REALSXP,n) );
	PROTECT( sum_log_fy_inside = allocVector(REALSXP,1) );
	PROTECT( coeff = allocVector(REALSXP,nZones) );
	PROTECT( y_temp = allocVector(REALSXP,n) );
	PROTECT( E_temp = allocVector(REALSXP,n) );
	PROTECT( inside = allocVector(REALSXP,n) );
	SXP = SXP + 6;

	// Compute n normalization terms (log negbinom)
	REAL(a)[0] = REAL(a_INPUT)[0];  REAL(b)[0] = REAL(b_INPUT)[0];
	for(i=0; i<n; i++){
		REAL(y)[0] = REAL(y_INPUT)[i];  REAL(E)[0] = REAL(E_INPUT)[i];
		REAL(log_fy_n)[i] = REAL(ldnbinom(y,E,a,b))[0];
	}

	// Main loop
	REAL(a)[0] = REAL(a_INPUT)[1];  REAL(b)[0] = REAL(b_INPUT)[1];
	for(j=0; j<nZones; j++){

		// Reset Variables
		REAL(yz)[0] = 0;  REAL(Ez)[0] = 0;
		for(i=0; i<n; i++){
			REAL(y_temp)[i] = 0;
			REAL(E_temp)[i] = 0;
			REAL(inside)[i] = 0;
		}

		for(k=0; k<length( VECTOR_ELT(cluster_list_INPUT,j) ); k++){
			// Watch off by 1 vector indexing in C as opposed to R
			index = INTEGER(VECTOR_ELT(cluster_list_INPUT,j))[k] - 1;
			REAL(yz)[0] = REAL(yz)[0] + REAL(y_INPUT)[ index ];
			REAL(Ez)[0] = REAL(Ez)[0] + REAL(E_INPUT)[ index ];
			REAL(y_temp)[ index ] = REAL(y_INPUT)[ index ];
			REAL(E_temp)[ index ] = REAL(E_INPUT)[ index ];
			REAL(inside)[ index ] = 1;
		}

		// Compute normalization terms inside the zone
		REAL(sum_log_fy_inside)[0] = 0;
		for(i=0; i<n; i++){
			if( REAL(inside)[i] == 1 ){
				REAL(sum_log_fy_inside)[0] += REAL(log_fy_n)[i];
			}
		}

		REAL(coeff)[j] = REAL(ldmultinom(y_temp,vecSum(y_temp),E_temp))[0] +
		REAL(ldnbinom(yz,Ez,a,b))[0] - REAL(sum_log_fy_inside)[0] ;
	}

	// output results as R list
	UNPROTECT(SXP);
	return(coeff);
}







/*
SEXP compute_coeff(SEXP clusterlist, SEXP coeff){
	PROTECT( clusterlist = coerceVector(clusterlist, INTSXP) );
	PROTECT( coeff = coerceVector(coeff, REALSXP) );

	int i,n,k,K;
	double sum;
	SEXP Rdim;
	PROTECT( Rdim = getAttrib(clusterlist, R_DimSymbol));
	n = INTEGER(Rdim)[1];
	K = INTEGER(Rdim)[0];
	UNPROTECT(1);

	SEXP output;
	PROTECT( output = allocVector(REALSXP, n) );

	for(i=0; i<n; i++){
		sum = 0;
		for(k=0; k<K; k++){
			sum = sum + REAL(coeff)[ INTEGER(clusterlist)[k + i*K] - 1];
		}
		REAL(output)[i] = sum;
	}


	UNPROTECT(3);
	return output;
}
*/



SEXP check_overlap(SEXP configs, SEXP overlap){
	PROTECT( configs = coerceVector(configs, INTSXP) );
	PROTECT( overlap = coerceVector(overlap, VECSXP) );

	int SXP = 2, k, n, i, j, m, l, index, index2;

	SEXP Rdim;
	PROTECT( Rdim = getAttrib(configs, R_DimSymbol));
	k = INTEGER(Rdim)[0];
	n = INTEGER(Rdim)[1];
	UNPROTECT(1);

	SEXP presence, cluster_list, indicator;
	PROTECT( presence = getListElement(overlap,"presence") );
	PROTECT( cluster_list = getListElement(overlap,"cluster_list") );
	PROTECT( indicator = allocVector(INTSXP, n) );
	for(i=0; i<n; i++)
		INTEGER(indicator)[i] = 1;
	SXP+=3;

	int n_zones = length(cluster_list), overlap_vector[n_zones];

	for(i=0; i<n; i++){
		memset( overlap_vector, 0, n_zones*sizeof(int) );

		for(j=0; j<k; j++){
			index = INTEGER(configs)[i*k + j] - 1;  // Watch R vs C Indexing

			if(overlap_vector[index]==1){
				INTEGER(indicator)[i] = 0;
				break;
			}


			for(m=0; m<length(VECTOR_ELT(cluster_list, index)); m++){
				index2 = INTEGER( VECTOR_ELT(cluster_list, index) )[m] - 1;  // Watch R vs C Indexing
				for(l=0; l<length(VECTOR_ELT(presence,index2)); l++)
					overlap_vector[ INTEGER(VECTOR_ELT(presence,index2))[l] - 1 ] = 1;
			}
		}
	}

	UNPROTECT(SXP);
	return indicator;
}




/*
SEXP return_moves(SEXP theta, SEXP overlap, SEXP cluster_coords, SEXP k_max, SEXP move){
	PROTECT( theta = coerceVector(theta, INTSXP) );
	if(length(overlap)==2)
		PROTECT( overlap = coerceVector(overlap, VECSXP) );
	else
		PROTECT( overlap = coerceVector(overlap, INTSXP) );
	PROTECT( cluster_coords = coerceVector(cluster_coords, INTSXP) );
	PROTECT( k_max = coerceVector(k_max, INTSXP) );
	PROTECT( move = coerceVector(move, INTSXP) );

	int SXP = 5, k = length(theta), n_zones = length(cluster_coords)/2,
	same_center[n_zones], diff_center[n_zones],
	i, j, z, l, m, index, index_sub, index_sub2, start, end, num_moves, counter;

	// Save Moves Here
	SEXP moves, output, possible_moves;

	PROTECT( possible_moves = allocVector(INTSXP, 1) );
	INTEGER(possible_moves)[0] = 1;
	SXP++;

	// Check If Use Large Version of Overlap
	int large_flag = 0, overlap_vector[n_zones];
	SEXP presence, cluster_list;
	if(length(overlap)==2){
		large_flag = 1;
		PROTECT( presence = getListElement(overlap,"presence") );
		PROTECT( cluster_list = getListElement(overlap,"cluster_list") );
		SXP += 2;
	}

	// Note:  R_NilValue == allocVector(NILSXP, 1)


	//-------------------------------------------------------------------------
	//
	// Moves 1-3: Growth/Trim/Recenter Moves Have Similar Structure:  Split on k=0 or not
	//
	//-------------------------------------------------------------------------
	// If k=0, no such moves
	if(( INTEGER(move)[0]==0 || INTEGER(move)[0]==1 || INTEGER(move)[0]==2 ) && k==0){
		PROTECT(output = R_NilValue);
		SXP++;
		INTEGER(possible_moves)[0] = 0;
	}

	// If k!=0, loop over all k unit zones in theta
	if(( INTEGER(move)[0]==0 || INTEGER(move)[0]==1 || INTEGER(move)[0]==2 ) && k!=0){
		PROTECT(moves = allocMatrix(INTSXP, n_zones, k));
		SXP++;
		memset( INTEGER(moves), 0, n_zones*k*sizeof(int) );

		for(j=0; j<k; j++){
			index = INTEGER(theta)[j] - 1; // Watch R vs C Indexing

			// Identify indices of smallest (start) & largest (end) unit zones with same center
			start = index;
			end = index + 1;
			while(INTEGER(cluster_coords)[start-1] == INTEGER(cluster_coords)[index] && start != 0)
				start--;
			while(INTEGER(cluster_coords)[end] == INTEGER(cluster_coords)[index] && end != n_zones)
				end++;

			// Identify unit zones with same/different center depending on move
			switch(INTEGER(move)[0]){
				case 0:
				case 1:
					memset( same_center, 0, n_zones*sizeof(int) );
					for(z=start; z<end; z++)
						same_center[z] = 1;

					// Check overlap with all OTHER unit zones
					for(i=0; i<k; i++)
						if(i != j){
							index_sub = INTEGER(theta)[i] - 1; // Watch R vs C Indexing

							if(large_flag==0){
								for(z=start; z<end; z++)
									same_center[z] *= 1-INTEGER(overlap)[n_zones*index_sub + z];
							}else{
								memset( overlap_vector, 0, n_zones*sizeof(int) );
								for(m=0; m<length(VECTOR_ELT(cluster_list, index_sub)); m++){
									index_sub2 = INTEGER( VECTOR_ELT(cluster_list, index_sub) )[m] - 1;
									for(l=0; l<length(VECTOR_ELT(presence,index_sub2)); l++)
										overlap_vector[ INTEGER(VECTOR_ELT(presence,index_sub2))[l] - 1 ] = 1;
								}
								for(z=start; z<end; z++)
									same_center[z] *= 1-overlap_vector[z];
							}

						}
					break;
				case 2:
					for(z=0; z<n_zones; z++)
						diff_center[z] = 1;

					for(z=start; z<end; z++)
						diff_center[z] = 0;

					// Check overlap with all OTHER unit zones
					for(i=0; i<k; i++)
						if(i != j){
							index_sub = INTEGER(theta)[i] - 1; // Watch R vs C Indexing

							if(large_flag==0){
								for(z=0; z<n_zones; z++)
									diff_center[z] *= 1-INTEGER(overlap)[n_zones*index_sub + z];
							}else{
								memset( overlap_vector, 0, n_zones*sizeof(int) );
								for(m=0; m<length(VECTOR_ELT(cluster_list, index_sub)); m++){
									index_sub2 = INTEGER( VECTOR_ELT(cluster_list, index_sub) )[m] - 1;
									for(l=0; l<length(VECTOR_ELT(presence,index_sub2)); l++)
										overlap_vector[ INTEGER(VECTOR_ELT(presence,index_sub2))[l] - 1 ] = 1;
								}
								for(z=0; z<n_zones; z++)
									diff_center[z] *= 1-overlap_vector[z];
							}

						}
					break;
			}

			// Store moves
			switch(INTEGER(move)[0]){
				case 0:
					for(z=index+1; z<end; z++)
						INTEGER(moves)[n_zones*j + z] = same_center[z];
					break;
				case 1:
					for(z=start; z<index; z++)
						INTEGER(moves)[n_zones*j + z] = same_center[z];
					break;
				case 2:
					for(z=0; z<n_zones; z++)
						INTEGER(moves)[n_zones*j + z] = diff_center[z];
					break;
			}
		} // end for k


		// Make table of all possible new theta's
		num_moves = INTEGER(vecSumInt(moves))[0];
		if(num_moves==0){
			PROTECT(output = R_NilValue);
			SXP++;
			INTEGER(possible_moves)[0] = 0;
		}else{
			PROTECT(output = allocMatrix(INTSXP, k, num_moves));
			SXP++;
			// Fill in old theta
			for(i=0; i<num_moves; i++)
				for(j=0; j<k; j++)
					INTEGER(output)[i*k + j] = INTEGER(theta)[j];

			// Adjust configs to include new single zone
			counter=0;
			for(j=0; j<k; j++)
				for(z=0; z<n_zones; z++)
					if(INTEGER(moves)[n_zones*j + z] == 1){
						INTEGER(output)[counter*k + j] = z + 1; // Watch R vs C Indexing
						counter++;
					}
			output = sort_matrix(output);

		}

	} // end if k!=0



	//---------------------------------------------------------------------
	// 4. Death Move
	//---------------------------------------------------------------------
	if(INTEGER(move)[0]==3){
		switch(k){
			case 0:
				PROTECT(output = R_NilValue);
				SXP++;
				INTEGER(possible_moves)[0] = 0;
				break;
			case 1:
				PROTECT(output = R_NilValue);
				SXP++;
				INTEGER(possible_moves)[0] = 1;
				break;
			default:
				PROTECT(output = allocMatrix(INTSXP, k-1, k));
				SXP++;
				for(j=0; j<k; j++){
					counter=0;
					for(i=0; i<k; i++)
						if(i!=j){
							INTEGER(output)[j*(k-1) + counter] = INTEGER(theta)[i];
							counter++;
						}
				}
				break;
		}
	}



	//---------------------------------------------------------------------
	// 5. Birth Move:  Split on k==0, k in middle, k==k_max
	//---------------------------------------------------------------------
	if(INTEGER(move)[0]==4){
		if(k == INTEGER(k_max)[0]){
			PROTECT(output = R_NilValue);
			SXP++;
			INTEGER(possible_moves)[0] = 0;
		}

		if(k == 0){
			PROTECT(output = allocMatrix(INTSXP, 1, n_zones));
			SXP++;
			for(z=0; z<n_zones; z++)
				INTEGER(output)[z] = z + 1;  // Watch R vs C Indexing
		}

		if(k != INTEGER(k_max)[0] && k !=0 ){
			PROTECT(moves = allocVector(INTSXP, n_zones));
			for(z=0; z<n_zones; z++)
				INTEGER(moves)[z] = 1;
			// Check Overlap
			for(j=0; j<k; j++){
				index = INTEGER(theta)[j] - 1;  // Watch R vs C Indexing

				if(large_flag==0){
					for(z=0; z<n_zones; z++)
						INTEGER(moves)[z] *= 1-INTEGER(overlap)[n_zones*index + z];
				}else{
					memset( overlap_vector, 0, n_zones*sizeof(int) );
					for(m=0; m<length(VECTOR_ELT(cluster_list, index)); m++){
						index_sub2 = INTEGER( VECTOR_ELT(cluster_list, index) )[m] - 1;
						for(l=0; l<length(VECTOR_ELT(presence,index_sub2)); l++)
							overlap_vector[ INTEGER(VECTOR_ELT(presence,index_sub2))[l] - 1 ] = 1;
					}
					for(z=0; z<n_zones; z++)
						INTEGER(moves)[z] *= 1-overlap_vector[z];
				}

			}
			SXP++;

			num_moves = INTEGER(vecSumInt(moves))[0];
			if(num_moves==0){
				PROTECT(output = R_NilValue);
				SXP++;
				INTEGER(possible_moves)[0] = 0;
			}else{
				PROTECT(output = allocMatrix(INTSXP, k+1, num_moves));
				SXP++;
				for(i=0; i<num_moves; i++)
					for(j=0; j<k; j++)
						INTEGER(output)[i*(k+1)+j] = INTEGER(theta)[j];

				counter=0;
				for(z=0; z<n_zones; z++){
					if( INTEGER(moves)[z]==1 ){
						INTEGER(output)[counter*(k+1)+k] = z + 1;  // Watch R vs C Indexing
						counter++;
					}
				}
				output = sort_matrix(output);
			}
		}
	}



	//-------------------------------------------------------------------------
	// Output Results
	//-------------------------------------------------------------------------
	//-------------------------------------------------------------------------
	// Output Results As R List
	//-------------------------------------------------------------------------
	SEXP results, results_names;
	PROTECT(results = allocVector(VECSXP, 2)); SXP++;

	PROTECT(results_names = allocVector(STRSXP, 2));
	SET_STRING_ELT(results_names, 0,  mkChar("output"));
	SET_STRING_ELT(results_names, 1,  mkChar("possible_moves"));
	setAttrib(results, R_NamesSymbol, results_names);
	UNPROTECT(1);

	SET_VECTOR_ELT(results, 0, output);
	SET_VECTOR_ELT(results, 1, possible_moves);

	UNPROTECT(SXP);
	return results;
}
*/



/*end of spatepi functions*/
