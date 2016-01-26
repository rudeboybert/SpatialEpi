#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;



// [[Rcpp::export]]
NumericVector normalize(NumericVector x) {
  int n = x.size();
  NumericVector y(n);
  double sum_x = sum(x);
  for(int i=0; i<n; ++i){
    y[i] = x[i]/sum_x;
	}
  return y;
}



// [[Rcpp::export]]
int NumericVectorEquality(NumericVector x, NumericVector y){  
  if(x.size() != y.size()) {
    return 0;
  }

  int vector_equal = 1;
  for(int i=0; i<x.size(); ++i) {
    if(x[i] != y[i]) {
      vector_equal = 0;
      break;
    }
  }
 
 return vector_equal;
}



// [[Rcpp::export]]
double ldnbinom(double y, double E, double a, double b) {
  double value = lgamma(y+a) - lgamma(y+1) - lgamma(a) + y*(log(E)-log(E+b)) +
  a*(log(b)-log(E+b));

  return value;
}



// [[Rcpp::export]]
double ldmultinom(NumericVector x, NumericVector prob) {
  int len = x.size();
	NumericVector p(len);
	double value = 0, sum_p = sum(prob);

	// normalize probabilities
	for(int i=0; i<len; ++i) {
		p[i] = prob[i]/sum_p;
	}

	value = lgamma(sum(x) + 1);
	for(int i=0; i<len; ++i) {
		if(p[i] != 0) {
			value -= lgamma(x[i] + 1);
			value += x[i] * log(p[i]);
		}
	}

	return value;
}



// [[Rcpp::export]]
NumericVector coeff(NumericVector y_vector, NumericVector E_vector, 
	NumericVector a_values, NumericVector b_values, List cluster_list) {
	
	double yz, Ez, sum_log_fy_inside;
	int n = y_vector.size(), n_zones = cluster_list.size(), index;
	NumericVector log_fy_n(n), coeff(n_zones), y_temp(n), E_temp(n), inside(n);

	// Compute n normalization terms (log negbinom)
	for(int i=0; i<n; ++i) {
		log_fy_n[i] = ldnbinom(y_vector[i], E_vector[i], a_values[0], b_values[0]);
	}

	// Main loop
	for(int j=0; j<n_zones; ++j) {
		// Reset Variables to 0
		yz = 0;  
		Ez = 0;
		for(int i=0; i<n; ++i){
			y_temp[i] = 0;
			E_temp[i] = 0;
			inside[i] = 0;
		}

		Rcpp::NumericVector areas = cluster_list[j];
		for(int k=0; k<areas.size(); ++k){
			// Watch off by 1 vector indexing in C as opposed to R
			index = areas[k] - 1;
			yz = yz + y_vector[ index ];
			Ez = Ez + E_vector[ index ];
			y_temp[ index ] = y_vector[ index ];
			E_temp[ index ] = E_vector[ index ];
			inside[ index ] = 1;
		}

		// Compute normalization terms inside the zone
		sum_log_fy_inside = 0;
		for(int i=0; i<n; ++i) {
			if( inside[i] == 1 ) {
				sum_log_fy_inside += log_fy_n[i];
			}
		}

		coeff[j] = ldmultinom(y_temp, E_temp) + 
		ldnbinom(yz, Ez, a_values[1], b_values[1]) - sum_log_fy_inside ;
	}

	return coeff;
}



// [[Rcpp::export]]
int ProbSampleReplace(NumericVector prob) {
	int n = prob.size();
	NumericVector x(n);
	
	for(int i=0; i<n; ++i) {
		x[i] = i;
	}
	
	RNGScope scope;
	int sample = RcppArmadillo::sample(x, 1, FALSE, prob)[0];  
	return sample;
}



// [[Rcpp::export]]
NumericVector check_overlap(NumericMatrix config, List overlap) {
  List presence = as<List>(overlap["presence"]);
  List cluster_list = as<List>(overlap["cluster.list"]);
  int nSim = config.nrow();
  int k = config.ncol(), index, index2;
  int nZones = cluster_list.size();

  // Set all to true
  NumericVector indicator(nSim), overlap_vector(nZones);
  for(int i=0; i<nSim; ++i) {
    indicator[i] = 1;
  }
  
  for(int i=0; i<nSim; ++i){
    // reset overlap vector
    for(int j=0; j<nZones; ++j){
      overlap_vector[j] = 0;
    }
    
    for(int j=0; j<k; ++j) {
      index = config(i, j) - 1;  // Watch R vs C Indexing
      if(overlap_vector[index]==1){
        indicator[i] = 0;
        break;
      }
      
      Rcpp::NumericVector areas = cluster_list[index];
      for(int m=0; m<areas.size(); ++m){
       // Watch R vs C Indexing
        index2 = areas[m] - 1;
        
        Rcpp::NumericVector stuff = presence[index2];
        for(int l=0; l<stuff.size(); ++l)
         overlap_vector[ stuff[l] - 1 ] = 1;
     }
   } 
 }
 return indicator;
}



// [[Rcpp::export]]
NumericMatrix clean_moves_matrix(NumericVector theta, NumericMatrix moves, int n_zones) {
  int n_moves = sum(moves), k = theta.size(), counter;
  NumericMatrix moves_reduced(k, n_moves);
  
  // Fill in old theta
  for(int i=0; i<n_moves; ++i) {
    for(int j=0; j<k; ++j) {
      moves_reduced(j, i) = theta[j];
    }
  }

  // Adjust configs to include new single zone
  counter=0;
  for(int j=0; j<k; ++j) {
    for(int z=0; z<n_zones; ++z) {
      if(moves(j, z) == 1){
        // Watch R vs C indexing
        moves_reduced(j, counter) = z + 1; 
        counter++;
      }
    }
  }
  
  /*
   * Sort values in each column
   */
   NumericVector values(k);
   for(int i=0; i<n_moves; ++i) {
     for(int j=0; j<k; ++j) {
     	values[j] = moves_reduced(j, i);
   	}
   	std::sort(values.begin(), values.end());
   	for(int j=0; j<k; ++j) {
   		moves_reduced(j, i) = values[j];
   	}
   }

  return moves_reduced;
}    



// [[Rcpp::export]]
NumericMatrix return_death_moves(NumericVector theta) {
  // Note this assumes that theta is not the null configuration.
  // i.e. != 0
	int counter, k = theta.size();

	if (k==1) {
    // If single zone config, return null config
		NumericMatrix death_moves(0, 0);
		return death_moves;
	} else {
		// If multiple zone config, identify all possible new configs with a single
		// zone dropped
		NumericMatrix death_moves(k-1, k);

		for(int j=0; j<k; ++j) {
			counter = 0;
			for(int i=0; i<k; ++i){
				if(i != j) {
					death_moves(counter, j) = theta[i];
					counter++;
				}
			}
		}
		return death_moves;
	}
}



// [[Rcpp::export]]
NumericMatrix return_birth_moves(NumericVector theta, List overlap) {
	List presence = as<List>(overlap["presence"]);
	List cluster_list = as<List>(overlap["cluster.list"]);

	int k = theta.size(), n_zones = cluster_list.size(), counter, zone, zone_area, 
	n_birth_moves;//, move_possible;
	
	NumericVector overlap_vector(n_zones);
  NumericMatrix birth_moves(1, n_zones);

	// Set all zones as possible birth moves
	for(int z=0; z<n_zones; ++z) {
		birth_moves(0, z) = 1;
	}

	// Check overlap with existing single zones
	for(int j=0; j<k; ++j) {
		// Identify zone within config, and all areas in the zone
    // Watch R vs C indexing
		zone = theta[j] - 1;  
		Rcpp::NumericVector zone_areas = cluster_list[zone];

      // Reset the vector that tracks if other zones overlap the current zone
		for(int z=0; z<n_zones; ++z) {
			overlap_vector[z] = 0;
		}
		
      // For all areas in zone, identify which other zones it is present in
		for(int m=0; m<zone_areas.size(); ++m) {
				// Watch R vs C indexing
			zone_area = zone_areas[m] - 1;
			
			Rcpp::NumericVector other_zones = presence[zone_area];
			for(int l=0; l<other_zones.size(); ++l) {
				overlap_vector[other_zones[l] - 1 ] = 1;
			}
			
        // Update which moves are possible using boolean algebra 
			for(int z=0; z<n_zones; ++z) {
				birth_moves(0, z) *= 1-overlap_vector[z];
			}
		}
	}
	 // end of check overlap


  //Make reduced table of all possible configs of birth moves
  n_birth_moves = sum(birth_moves);
  //move_possible = n_birth_moves > 0;
  NumericMatrix birth_moves_reduced(k+1, n_birth_moves);
  
  // Fill in old theta
  for(int i=0; i<n_birth_moves; ++i) {
    for(int j=0; j<k; ++j) {
      birth_moves_reduced(j, i) = theta[j];
    }
  }
  
  // Adjust configs to include new single zone
  counter=0;
  for(int z=0; z<n_zones; ++z) {
    if(birth_moves(0, z) == 1){
      // Watch R vs C indexing
      birth_moves_reduced(k, counter) = z + 1; 
      counter++;
    }
  }
    
    
  /*
   * Sort values in each column
   */
   NumericVector values(k+1);
   for(int i=0; i<n_birth_moves; ++i) {
     for(int j=0; j<k+1; ++j) {
       values[j] = birth_moves_reduced(j, i);
   	}
   	std::sort(values.begin(), values.end());
   	for(int j=0; j<k+1; ++j) {
   		birth_moves_reduced(j, i) = values[j];
   	}
   }

	/*
   * Output results
   */
   return birth_moves_reduced;
//   List results = List::create(_["birth_moves"] = birth_moves_reduced, 
//   	_["move_possible"] = move_possible);
//   return results;
}
    
    
    
// [[Rcpp::export]]
List return_local_moves(NumericVector theta, List overlap, NumericMatrix 
	cluster_coords) {

	List presence = as<List>(overlap["presence"]);
	List cluster_list = as<List>(overlap["cluster.list"]);

	int k = theta.size(), n_zones = cluster_list.size(), zone, start, end, zone_sub,
	zone_sub2, center;

	NumericVector overlap_vector(n_zones), diff_center(n_zones), 
  same_center(n_zones), 
  move_possible(3);

  NumericMatrix growth_moves(k, n_zones), trim_moves(k, n_zones),
  replace_moves(k, n_zones), 
  growth_moves_reduced, trim_moves_reduced, replace_moves_reduced;


  /*
   * Loop over all component single zones in config 
   */
	for(int j=0; j<k; ++j){
		// Identify single zone in question and its center
		// Watch R vs C indexing
		zone = theta[j] - 1; 
   	center = cluster_coords(zone, 0);  

		// Identify all single zones that the same/different centering area.  
		// Note: [start, end) are the indices of all single zones with the same 
		// center
		for(int i=0; i<n_zones; ++i) {
			if(cluster_coords(i, 0) == center) {
				start = i;
				break;
   		}
   	}
   	for(int i=n_zones-1; i>=0; --i) {
   		if(cluster_coords(i, 0) == center) {
   			end = i+1;
   			break;
   		}
   	}
   	for(int z=0; z<n_zones; ++z) {
   		diff_center[z] = 1;
   		same_center[z] = 0;
   	}
   	for(int z=start; z<end; z++){
   		same_center[z] = 1;
   		diff_center[z] = 0;
   	}

    // Find overlap of single zone with all other single zones
   	for(int i=0; i<k; ++i) {
   		if(i != j){
				// Watch R vs C indexing
   			zone_sub = theta[i] - 1;

				// Reset the vector that tracks if other zones overlap the current zone
   			for(int z=0; z<n_zones; ++z) {
   				overlap_vector[z] = 0;   
   			}

   			Rcpp::NumericVector areas = cluster_list[zone_sub];
   			for(int m=0; m<areas.size(); ++m){
   				zone_sub2 = areas[m] - 1;
   				Rcpp::NumericVector flags = presence[zone_sub2];
   				for(int l=0; l<flags.size(); ++l) {
   					overlap_vector[flags[l] - 1] = 1;
   				}
   			}

				// Update same center
   			for(int z=start; z<end; ++z) {
   				same_center[z] *= 1-overlap_vector[z];
   			}
				// Update different center
   			for(int z=0; z<n_zones; ++z) {
   				diff_center[z] *= 1-overlap_vector[z];
   			}
   		} // end if(i!=j)

			// Store moves
      // Growth
   		for(int z=zone+1; z<end; ++z) {
   			growth_moves(j, z) = same_center[z];
   		}
			// Trim
   		for(int z=start; z<zone; ++z) {
   			trim_moves(j, z) = same_center[z];
   		}
      // Replace
   		for(int z=0; z< n_zones; ++z) {
   				replace_moves(j, z) = diff_center[z];
   		}
   	} // end for(int i=0; i<k; ++i)
  } // end of loop over component single zones


  //Make reduced table of all possible configs.  Same for all local moves
  move_possible[0] = sum(growth_moves) > 0;
  move_possible[1] = sum(trim_moves) > 0;
  move_possible[2] = sum(replace_moves) > 0;
   
  growth_moves_reduced = clean_moves_matrix(theta, growth_moves, n_zones);
  trim_moves_reduced = clean_moves_matrix(theta, trim_moves, n_zones);
  replace_moves_reduced = clean_moves_matrix(theta, replace_moves, n_zones);

  /*
   * Output results
   */
   List results = List::create(
   	_["growth_moves"] = growth_moves_reduced, 
    _["trim_moves"] = trim_moves_reduced, 
    _["replace_moves"] = replace_moves_reduced, 
    _["move_possible"] = move_possible
   	);

   return results;  
}



// [[Rcpp::export]]
List MCMC_simulation(int n_sim, NumericVector pattern, NumericVector theta_init, 
	List overlap, NumericMatrix cluster_coords, NumericVector p_moves_orig, 
	int J, NumericVector lkhd_z, NumericVector lambda){
	
	int n_zones = cluster_coords.nrow(), z_star, z,
  k, k_rev, 
  move, move_rev,
  n_moves, n_moves_rev;

	NumericVector reverse = NumericVector::create(1, 0, 2, 4, 3);  
	double p_num, p_denom, u;

	// Store Outputs Here
	NumericVector move_trace(n_sim), accpt_trace(n_sim), ratio(n_sim);
	List sample(n_sim);

	//----------------------------------------------------------
	// Sampling
	//----------------------------------------------------------
	// Either Uniform or Proportional Sampling
	NumericVector unif_p(n_zones), prop_p(n_zones);

  prop_p = normalize(lkhd_z);
	for(int i=0; i<n_zones; ++i){
		unif_p[i] = 1.0/n_zones;
	}

	//----------------------------------------------------------
	// Main Loop
	//----------------------------------------------------------
  // initial value
	sample[0] = theta_init;
  
	for(int sim=1; sim<n_sim; ++sim) {
 
		// Select Sampling Type:  Remove Later?
    NumericVector p(n_zones);
		int sampling_type = pattern[ sim % pattern.size() ];
		if(sampling_type == 0) {
      for(int i=0; i<n_zones; ++i)
        p[i] = unif_p[i];
		} else {
      for(int i=0; i<n_zones; ++i)
        p[i] = prop_p[i];
		}

    // Get theta info 
    NumericVector theta = sample[sim-1];
    k = theta.size();


    //---------------------------------------------------------------------
    // Forward moves
    //--------------------------------------------------------------------- 
    // Reset forward values
    p_num = 1;
    NumericVector p_moves(p_moves_orig.size());
    for(int i=0; i<p_moves_orig.size(); ++i)
  	  p_moves[i] = p_moves_orig[i];
    
    //-----------------------------------------------------    
		// Figure out legal forward moves, get their probs, sample
    //-----------------------------------------------------    
    // k==0: only birth move
		if(k==0){
      for(int i=0; i<p_moves.size(); ++i)
        p_moves[i] = 0;
			p_moves[4] = 1;
		}
    // k==J: no birth move
		if(k==J){
			p_moves[4] = 0;
		}  
    // Figure out if local moves are possible
    List local_moves = return_local_moves(theta, overlap, cluster_coords);    
    Rcpp::NumericVector local_move_possible = local_moves["move_possible"];
    for(int i=0; i<local_move_possible.size(); ++i) {
      if(local_move_possible[i] == 0) {
        p_moves[i] = 0;
      }
    }  
    // Figure out if birth moves are possible
    NumericMatrix birth_moves = return_birth_moves(theta, overlap);
    if(birth_moves.ncol() == 0) {
      p_moves[4] = 0;
    }

		// Renormalize prob and sample move
    p_moves = normalize(p_moves); 
    
		move = ProbSampleReplace(p_moves);
		p_num *= p_moves[move];

    
    //-----------------------------------------------------    
		// Set up sampling space of configurations
    NumericMatrix moves;
		// Local moves
		if(move == 0 || move == 1 || move == 2){
      moves = as<NumericMatrix>(local_moves[move]);
		}
    // Death moves
		if(move == 3) {
      moves = return_death_moves(theta);
		}
		// Birth moves
		if(move == 4) {
      moves = birth_moves;        
		}


    //--------------------------------------------------
		// Sample configs: 
		//--------------------------------------------------
    k_rev = moves.nrow();
    NumericVector theta_star(k_rev);
    // If k==1 and we have death move, treat differently
		if(move==3 && k==1){
			p_num *= 1;
		} else {
			// Get number of possible forward proposal configurations and k*
			n_moves = moves.ncol();

			// Compute probability of each possible forward proposed configuration
			NumericVector probs(n_moves);
			for(int i=0; i<n_moves; ++i) {
				probs[i] = 1;
				for(int j=0; j<k_rev; ++j) {
					// Watch R vs C indexing
					probs[i] *= p[moves(j, i) - 1];
				}
			}
			probs = normalize(probs);
      
			// Sample a new configuration theta_star and sort it
      z_star = ProbSampleReplace(probs);
			for(int j=0; j<k_rev; ++j){
				theta_star[j] = moves(j, z_star); 
			}
      
      // Update prob
      p_num *= probs[z_star];
      
		} // end if/else if(move==3 && k==1)



		//-----------------------------------------------------------------
		// Backward moves
		//-----------------------------------------------------------------
    // Reset forward values
    p_denom = 1;
    NumericVector p_moves_rev(p_moves_orig.size());
    for(int i=0; i<p_moves_orig.size(); ++i)
      p_moves_rev[i] = p_moves_orig[i];
    
  	//------------------------------------------------------
		// Figure out legal backward moves and get their probs
		//------------------------------------------------------
		// k_rev==0: only birth moves
    if(k_rev==0){
      for(int i=0; i<p_moves_rev.size(); ++i)
        p_moves_rev[i] = 0;
			p_moves_rev[4] = 1;
		}
		// k_rev==J: no death moves
		if(k_rev==J){
			p_moves_rev[4] = 0;
		}
    // Figure out if local moves are possible
    List local_rev_moves = return_local_moves(theta_star, overlap, cluster_coords);    
    Rcpp::NumericVector local_rev_move_possible = local_rev_moves["move_possible"];
    for(int i=0; i<3; ++i) {
      if(local_rev_move_possible[i] == 0) {
        p_moves_rev[i] = 0;
      }
    }  
    // Figure out if birth moves are possible
    NumericMatrix birth_rev_moves = return_birth_moves(theta_star, overlap);
    if(birth_rev_moves.ncol() == 0){
      p_moves_rev[4] = 0;
    }
    
		// Renormalize prob, set determined reverse move
    p_moves_rev = normalize(p_moves_rev);  
		move_rev = reverse[move];
		p_denom *= p_moves_rev[move_rev];
    
		//--------------------------------------------------
		// Set up sampling space of configurations
		//--------------------------------------------------
    NumericMatrix moves_rev;
    // Local moves
		if(move_rev == 0 || move_rev == 1 || move_rev == 2){
      moves_rev = as<NumericMatrix>(local_rev_moves[move_rev]);
		}
    // Death moves
		if(move_rev == 3) {
      moves_rev = return_death_moves(theta_star);
		}
		// Birth moves
		if(move_rev == 4) {
      moves_rev = birth_rev_moves;        
		}


    //--------------------------------------------------
  	// Compute probability of reverse move
    //--------------------------------------------------
    // If k_rev==1 and we have reverse death move, treat differently
		if(move_rev==3 && k_rev==1){
			p_denom *= 1;
		} else {
			// Get number of possible forward proposal configurations and k*
			n_moves_rev = moves_rev.ncol();

			// Compute probability of each possible reverse proposed configuration
			NumericVector probs(n_moves_rev);
			for(int i=0; i<n_moves_rev; ++i) {
				probs[i] = 1;
				for(int j=0; j<k; ++j) {
					// Watch R vs C indexing
					probs[i] *= p[moves_rev(j, i) - 1];
				}
			}
			probs = normalize(probs);

			// Find original configuration theta FIX THIS
  	  for(z=0; z<n_moves_rev; ++z){
        NumericVector config = moves_rev( _, z);    
        if(NumericVectorEquality(config, theta))
          break;
		  }
      
      // Update prob
      p_denom *= probs[z]; 

		} // end if/else if(move==3 && k==1)

//    printf("%u %u %u \n%4.6f %4.6f %4.6f %4.6f %4.6f \n%4.6f %4.6f %4.6f %4.6f %4.6f \n%4.6f %4.6f\n", 
//    move+1, z_star, move_rev+1,
//    p_moves[0], p_moves[1], p_moves[2], p_moves[3], p_moves[4],
//    p_moves_rev[0], p_moves_rev[1], p_moves_rev[2], p_moves_rev[3], p_moves_rev[4],
//    p_num, p_denom
//    );
    
  	//-----------------------------------------------------------------
    // M-H Ratio
  	//-----------------------------------------------------------------    
    ratio[sim-1] = p_denom/p_num;
		ratio[sim-1] *= lambda[k_rev];
  	ratio[sim-1] /= lambda[k]; 
        
		if(k_rev != 0) {
			for(int j=0; j<k_rev; ++j) {
				// Watch R vs C indexing
				ratio[sim-1] *= lkhd_z[theta_star[j]-1]; 
			}
		}
		if(k != 0) {
			for(int j=0; j<k; ++j) {
				// Watch R vs C indexing
				ratio[sim-1] /= lkhd_z[theta[j]-1]; 
			}
		}

    // Accept/Reject:  As ratio goes up, more chance of accepting
		u = runif(1)[0];
		if( u < ratio[sim-1] ){
			sample[sim] = theta_star;
			accpt_trace[sim-1] = 1;
		} else {
			sample[sim] = theta;
		}
		
		// Watch R vs C indexing
		move_trace[sim-1] = move + 1;

  
  } // end overall n_sim
  
  return List::create(
    _["sample"] = sample, _["move_trace"] = move_trace,
    _["accpt_trace"] = accpt_trace, _["ratio"] = ratio
    ); 
}

