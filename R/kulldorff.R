#' Obtain all areas included in a single zone
#'
#' @param cluster_index integer index of single zone
#' @param zone_info output of \code{\link{define_single_zones}} function
#'
#' @return vector of areas
#' @export
#'
#' @examples
#' 1+1
return_single_zone_areas <- function(cluster_index, zone_info){
  # Obtain single zone center and radial areas
  center <- zone_info$cluster_coords[cluster_index, 1]
  radial <- zone_info$cluster_coords[cluster_index, 2]
 
  # Obtain all areas (in order of distance) from center to radial area
  cluster <- zone_info$nearest_neighbors[[center]]
  cluster <- cluster[1:which(cluster == radial)]
  
  return(cluster)
}


get_cluster_report <- function(cluster_index, 
                             population, cases, expected_cases,
                             zone_info,
                             all_log_lkhd, lambdas){
  cluster <- return_single_zone_areas(cluster_index, zone_info)

  output <- list(
    location_IDs_included = cluster,
    population = sum(population[cluster]),	
    number_of_cases = sum(cases[cluster]),
    expected_cases = sum(expected_cases[cluster]),
    SMR = sum(cases[cluster])/sum(expected_cases[cluster]),
    log_lkhd_ratio = all_log_lkhd[cluster_index], 
    monte_carlo_rank = sum(lambdas >= all_log_lkhd[cluster_index]),
    p_value = 1-mean(lambdas < all_log_lkhd[cluster_index])
  )
}



kulldorff <- function(centroids, cases, population, expected_cases = NULL, 
                      pop_upper_bound, alpha_level, n_sim = 9999, plot = TRUE){
  # Get geographic information
  zone_info <- define_single_zones(centroids, population, pop_upper_bound)
  
  # Depending on whether expected_cases were specified, set Kulldorff method
  # to either binomial or poisson and appropriate denominator
  if(is.null(expected_cases)){
    type <- "binomial"
    denominator <- population
    expected_cases <- sum(cases) * (denominator/sum(denominator))
    # poisson case	
  }else{
    type <- "poisson"
    denominator <- expected_cases
  }
  

  
  #-------------------------------------------------------------------------------
  # Compute Monte Carlo randomized p-value
  #-------------------------------------------------------------------------------
  # Observed statistic computation
  all_log_lkhd <- computeAllLogLkhd(cases, denominator, 
                                    zone_info$nearest_neighbors, type)
  
  # Simulate cases under null hypothesis of no area effects Compute simulated 
  # lambda's:  max log-lkhd in region
  lambdas <- 
    rmultinom(n=n_sim, size=round(sum(cases)), prob=denominator) %>% 
    kulldorffMC(., denominator, zone_info$nearest_neighbors, type) %>% 
    c(., max(all_log_lkhd))
  p_value <- 1-mean(lambdas < max(all_log_lkhd))
  
  #-------------------------------------------------------------------------------
  # Create Most Likely Cluster Object
  #-------------------------------------------------------------------------------
  cluster_index <- which.max(all_log_lkhd)
  most_likely_cluster <- 
    get_cluster_report(cluster_index, population, cases, expected_cases, 
                       zone_info, all_log_lkhd, lambdas)

  # Plot histogram
  if(plot){
    hist(lambdas, main="", xlab=expression(log(lambda)))
    title(expression(paste("Monte Carlo Distribution of ", log(lambda))))
    abline(v=max(all_log_lkhd), col="red")
    obs_log_lkhd <- max(all_log_lkhd) %>% round(3)
    p_value <- most_likely_cluster$p_value %>% round(log10(n_sim + 1))
    leg_lab <- paste(c("Observed value =", "p-value ="), c(obs_log_lkhd, p_value))
    legend("top", legend=leg_lab, lty=c(1, 1), col=c("red", "white"), 
           bty="n")
  }
  
  
  #-------------------------------------------------------------------------------
  # Investigate Secondary Clusters
  #-------------------------------------------------------------------------------
  # running list of areas already covered to test overlap
  current_clusters <- most_likely_cluster$location_IDs_included
  secondary_clusters <- NULL
  
  # go through log-likelihoods in decreasing order, skipping largest which 
  # corresponds to most likely cluster
  indices <- order(all_log_lkhd, decreasing=TRUE)
  
  for(i in 2:length(indices)){
    # Get areas included in new potential cluster
    cluster_index <- indices[i]
    new_cluster <- return_single_zone_areas(cluster_index, zone_info)
    
    # if there is no overlap between existing clusters and the new potential cluster
    if(length(intersect(new_cluster, current_clusters)) == 0){		
      new_secondary_cluster <- 
        get_cluster_report(cluster_index, population, cases, expected_cases, 
                           zone_info, all_log_lkhd, lambdas)
      
      # If no longer significant, end loop.
      if(new_secondary_cluster$p_value > alpha_level){
        break
      }
      
      # Otherwise add to existing clusters
      secondary_clusters <- c(secondary_clusters, list(new_secondary_cluster))
      current_clusters <- unique(c(current_clusters, new_cluster))
    }		
  }
  
  
  #-------------------------------------------------------------------------------
  # Output results
  #-------------------------------------------------------------------------------
  results <- list(
    most_likely_cluster = most_likely_cluster,
    secondary_clusters = secondary_clusters,
    type = type,
    log_lkhd = all_log_lkhd,
    simulated_log_lkhd = lambdas
  )
  return(results)
}




