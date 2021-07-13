#' Generate cluster reports
#'
#' Return a report for a given cluster containing the same information as SatScan
#'
#' @inheritParams define_single_zones
#' @inheritParams return_single_zone_areas
#' @inheritParams kulldorff
#' @param all_log_lkhd vector of the log likelihood of each single zone
#' @param lambdas vector of simulated max log likelihood
#'
#' @return For a particular clusters a report including \describe{ 
#'   \item{\code{location_IDs_included}}{ID's of areas in cluster, in order of distance} 
#'   \item{\code{population}}{population of cluster} 
#'   \item{\code{number_of_cases}}{number of cases in cluster}
#'   \item{\code{expected_cases}}{expected number of cases in cluster}
#'   \item{\code{SMR}}{estimated SMR of cluster}
#'   \item{\code{log_lkhd_ratio}}{log-likelihood of cluster}
#'   \item{\code{monte_carlo_rank}}{rank of lkhd of cluster within Monte Carlo simulated values}
#'   \item{\code{p_value}}{Monte Carlo \eqn{p}-value}
#'   }
#' 
#' @export
#'
#' @examples
#' 1+1
get_cluster_report <- function(single_zone_index, zone_info, cases, population,  
                               expected_cases, all_log_lkhd, lambdas){
  # Get areas included in single zone
  cluster_areas <- return_single_zone_areas(single_zone_index, zone_info)
  
  # Return information
  output <- list(
    location_IDs_included = cluster_areas,
    population = sum(population[cluster_areas]),	
    number_of_cases = sum(cases[cluster_areas]),
    expected_cases = sum(expected_cases[cluster_areas]),
    SMR = sum(cases[cluster_areas])/sum(expected_cases[cluster_areas]),
    log_lkhd_ratio = all_log_lkhd[single_zone_index], 
    monte_carlo_rank = sum(lambdas >= all_log_lkhd[single_zone_index]),
    p_value = 1-mean(lambdas < all_log_lkhd[single_zone_index])
  )
}



#' Kulldorff method for cluster detection
#'
#' @inheritParams define_single_zones
#' @param cases Vector of number of cases for each area
#' @param expected_cases Vector of number of expected cases for each area
#' @param alpha_level alpha-level threshold used to declare significance
#' @param n_sim number of Monte Carlo samples used for significance measures
#' @param plot logical variable of whether to plot Monte Carlo distribution
#'
#' @return A list containing \describe{ 
#'   \item{\code{most_likely_cluster}}{\code{\link{get_cluster_report}} report on most likely cluster} 
#'   \item{\code{secondary_clusters}}{list of \code{\link{get_cluster_report}} reports on secondary clusters; if none \code{NULL} is returned} 
#'   \item{\code{type}}{type of likelihood: binomial or poisson}
#'   \item{\code{log_lkhd}}{log-likelihood of each zone considered}
#'   \item{\code{simulated_log_lkhd}}{\code{n.simulations} Monte Carlo samples of the log-likelihood of the most likely cluster}
#'   }
#' @export
#'
#' @examples
#' data("NYleukemia")
#' centroids <- sp::coordinates(NYleukemia)
#' output <- kulldorff(centroids = centroids, 
#'    cases = NYleukemia$cases, population = NYleukemia$population, 
#'    expected_cases = NULL, pop_upper_bound = 0.15, alpha_level=0.05,
#'    n_sim=999)
#' output$most_likely_cluster
kulldorff <- function(centroids, cases, population, expected_cases = NULL, 
                      pop_upper_bound, alpha_level, n_sim = 9999, plot = FALSE){
  # Get geographic information
  zone_info <- define_single_zones(centroids, population, pop_upper_bound)
  
  # Depending on whether expected_cases were specified, set Kulldorff method
  # to use either binomial or Poisson likelihood and appropriate denominator
  if(is.null(expected_cases)){
    type <- "binomial"
    denominator <- population
    expected_cases <- sum(cases) * (denominator/sum(denominator))
    # poisson case	
  }else{
    type <- "poisson"
    denominator <- expected_cases
  }
  
  # Compute Monte Carlo p-value: After computing the log likelihood for each 
  # single zone, simulate n_sim cases under the null, and construct monte carlo
  # distribution
  all_log_lkhd <- computeAllLogLkhd(cases, denominator, 
                                    zone_info$nearest_neighbors, type)
  lambdas <- 
    rmultinom(n=n_sim, size=round(sum(cases)), prob=denominator) %>% 
    kulldorffMC(., denominator, zone_info$nearest_neighbors, type) %>% 
    c(., max(all_log_lkhd))
  p_value <- 1-mean(lambdas < max(all_log_lkhd))

  # Create Most Likely Cluster Object
  cluster_index <- which.max(all_log_lkhd)
  most_likely_cluster <- 
    get_cluster_report(cluster_index, zone_info, cases, population, 
                       expected_cases, 
                       all_log_lkhd, lambdas)

  # Plot histogram
  if(plot){
    hist(lambdas, main="", xlab=expression(log(lambda)))
    title(expression(paste("Monte Carlo Distribution of ", log(lambda))))
    abline(v=max(all_log_lkhd), col="red")
    obs_log_lkhd <- max(all_log_lkhd) %>% round(3)
    p_value <- most_likely_cluster$p_value %>% round(log10(n_sim + 1))
    leg_lab <- paste(c("Observed value =", "p-value ="), c(obs_log_lkhd, p_value))
    legend("top", legend=leg_lab, lty=c(1, 1), col=c("red", "white"), bty="n")
  }
  
  
  # Investigate Secondary Clusters
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
        get_cluster_report(cluster_index, zone_info, cases, population, expected_cases, 
                           all_log_lkhd, lambdas)
      
      # If no longer significant, end loop.
      if(new_secondary_cluster$p_value > alpha_level){
        break
      }
      
      # Otherwise add to existing clusters
      secondary_clusters <- c(secondary_clusters, list(new_secondary_cluster))
      current_clusters <- unique(c(current_clusters, new_cluster))
    }		
  }
  
  # Output results
  results <- list(
    most_likely_cluster = most_likely_cluster,
    secondary_clusters = secondary_clusters,
    type = type,
    log_lkhd = all_log_lkhd,
    simulated_log_lkhd = lambdas
  )
  return(results)
}



#' Plot results of \code{kulldorff} method
#' 
#' Wrapper function for \code{\link{map_variable}} to plot the results of a kulldorff analysis.  
#'
#' @param kulldorff_output List output of \code{\link{kulldorff}}
#' @param ... other arguments passed on to  \code{\link{map_variable}}
#' @inheritParams map_variable
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' data("NYleukemia")
#' centroids <- sp::coordinates(NYleukemia)
#' output <- kulldorff(centroids = centroids, 
#'    cases = NYleukemia$cases, population = NYleukemia$population, 
#'    expected_cases = NULL, pop_upper_bound = 0.15, alpha_level=0.05,
#'    n_sim=999)
#' map_kulldorff(NYleukemia, output)
map_kulldorff <- function(sp_obj, kulldorff_output, ...){
  # Get cluster information
  clusters <- c(
    list(kulldorff_output$most_likely_cluster),
    kulldorff_output$secondary_clusters
  )
  p_values <- sapply(clusters, function(x){x$p_value})
  
  # Insert p-values
  sp_obj$p_values <- rep(NA, length(sp_obj))
  for(i in 1:length(p_values)){
    areas <- clusters[[i]]$location_IDs_included
    sp_obj$p_values[areas] <- p_values[i]
  }
  map_ggplot <- map_variable(sp_obj, variable_name = "p_values", type="discrete", ...)
  
  return(map_ggplot)
}







