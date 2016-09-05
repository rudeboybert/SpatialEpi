#' Besag-Newell Cluster Detection Method
#' 
#' Besag-Newell cluster detection method.  There are differences with the original paper and our implementation:
#'
#' @param geo an \code{n x 2} table of the (x,y)-coordinates of the area centroids
#' @param population aggregated case counts for all \code{n} areas
#' @param cases aggregated population counts for all \code{n} areas
#' @param expected.cases expected numbers of disease for all \code{n} areas
#' @param k number of cases to consider
#' @param alpha.level alpha-level threshold used to declare significanc
#'
#' @details For the \code{population} and \code{cases} tables, the rows are 
#' bunched by areas first, and then for each area, the counts for each strata 
#' are listed.  It is important that the tables are balanced:  the strata 
#' information are in the same order for each area, and counts for each 
#' area/strata combination appear exactly once (even if zero
#' @return List containing
#' \item{clusters}{information on all clusters that are \eqn{\alpha}-level 
#' significant, in decreasing order of the \eqn{p}-value}
#' \item{p.values}{for each of the \eqn{n} areas, \eqn{p}-values of each cluster 
#' of size at least \eqn{k}}
#' \item{m.values}{for each of the \eqn{n} areas, the number of areas need to 
#' observe at least \eqn{k} cases}
#' \item{observed.k.values}{based on \code{m.values}, the actual number of cases 
#' used to compute the \eqn{p}-values}
#' @export
#' @seealso \code{\link{pennLC}}
#' @seealso \code{\link{expected}}
#' @references Besag J. and Newell J. (1991) The Detection of Clusters in Rare 
#' Diseases \emph{Journal of the Royal Statistical Society. Series A (Statistics 
#' in Society)}, \bold{154}, 143--155
#'
#' @examples
#' # Load Pennsylvania Lung Cancer Data
#' data(pennLC)
#' data <- pennLC$data
#' 
#' # Process geographical information and convert to grid
#' geo <- pennLC$geo[,2:3]
#' geo <- latlong2grid(geo)
#' 
#' # Get aggregated counts of population and cases for each county
#' population <- tapply(data$population,data$county,sum)
#' cases <- tapply(data$cases,data$county,sum)
#' 
#' # Based on the 16 strata levels, computed expected numbers of disease
#' n.strata <- 16
#' expected.cases <- expected(data$population, data$cases, n.strata)
#' 
#' # Set Parameters
#' k <- 1250
#' alpha.level <- 0.05
#' 
#' # not controlling for stratas
#' results <- besag_newell(geo, population, cases, expected.cases=NULL, k, alpha.level)
#' # controlling for stratas
#' results <- besag_newell(geo, population, cases, expected.cases, k, alpha.level)
besag_newell <-
  function(geo, population, cases, expected.cases=NULL, k, alpha.level){
    
    #-------------------------------------------------------------------------------
    # Initialization 
    #-------------------------------------------------------------------------------
    # If no expected.cases provided, set them if there are no expected counts
    if(is.null(expected.cases)){
      p <- sum(cases)/sum(population)
      expected.cases <- population*p
    }
    
    # geographical information computation
    geo.results <- zones(geo, population, 1)
    nearest.neighbors <- geo.results$nearest.neighbors
    distance <- geo.results$dist
    n.zones <- length(unlist(nearest.neighbors))
    
    
    #-------------------------------------------------------------------------------
    # Observed statistic computation
    #-------------------------------------------------------------------------------
    results <- besag_newell_internal(cases, expected.cases, nearest.neighbors, 
                                     n.zones, k)
    
    # observed p.values for each areas
    p.values <- results$observed.p.values	
    # observed number of neighbors needed to observe k cases
    m.values <- results$observed.m.values	
    # actual observed number of cases
    k.values <- results$observed.k.values	
    
    # pick out areas that were significant and order them by p-value
    signif.indices <- order(p.values)[1:sum(p.values <= alpha.level)]
    
    # order remaining values
    signif.p.values <- p.values[signif.indices]
    signif.m.values <- m.values[signif.indices]
    signif.k.values <- k.values[signif.indices]
    
    
    # Create object to output
    # If none are significant, return NULL
    if(length(signif.indices) == 0){
      clusters <- NULL
    } else {
      clusters <- vector("list", length=length(signif.indices))
      
      for( i in 1:length(clusters) ){	
        # find areas included in cluster
        cluster <- order(distance[signif.indices[i],])[1:signif.m.values[i]]
        
        new.cluster <- list(
          location.IDs.included = cluster,
          population = sum(population[cluster]),
          number.of.cases = sum(cases[cluster]),
          expected.cases = sum(expected.cases[cluster]),
          SMR = sum(cases[cluster])/sum(expected.cases[cluster]),
          p.value = signif.p.values[i]	
        )
        clusters[[i]] <- new.cluster
      }
    }
    
    
    #-------------------------------------------------------------------------------
    # Output results
    #-------------------------------------------------------------------------------
    results <- list(
      clusters=clusters,
      p.values=p.values,
      m.values=m.values,
      observed.k.values=k.values
    )	
    return(results)
  }




#' Compute Expected Numbers of Disease
#' 
#' Compute the internally indirect standardized expected numbers of disease. 
#' The \code{population} and \code{cases} vectors must be \emph{balanced}: all 
#' counts are sorted by area first, and then within each area the counts for all 
#' strata are listed (even if 0 count) in the same order.
#' 
#' @param population a vector of population counts for each strata in each area
#' @param cases a vector of the corresponding number of cases
#' @param n.strata number of strata considered
#'
#' @return A vector of expected counts
#' \item{expected.cases}{a vector of the expected numbers of disease for each area}
#' @export
#' @references Elliot, P. et al. (2000) \emph{Spatial Epidemiology:  
#' Methods and Applications}.  Oxford Medical Publications.
#' @examples
#' data(pennLC)
#' population <- pennLC$data$population
#' cases <- pennLC$data$cases
#' # In each county in Pennsylvania, there are 2 races, gender and 4 age bands 
#' # considered = 16 strata levels
#' pennLC$data[1:16,]
#' expected(population, cases, 16)
expected <-
  function(population, cases, n.strata){
    
    n <- length(population)/n.strata
    E <- rep(0, n)
    qNum <- rep(0, n.strata)
    qDenom <- rep(0, n.strata)
    q <- rep(0, n.strata)
    
    
    # Compute q: strata-specific rates. Numerator and denominator separately
    for(i in 1:n.strata){
      indices <- rep(i, n) + seq(0, n-1)*n.strata
      qNum[i] <- qNum[i] + sum(cases[indices])
      qDenom[i] <- qDenom[i] + sum(population[indices])
    }
    q <- qNum/qDenom
    
    
    # Compute E expected counts
    for(i in 1:n) {
      indices <- 1:n.strata + (i-1)*n.strata
      E[i] <- sum(population[indices]*q)
    }
    
    return(E)
  }



#' Title
#'
#' @param geo 
#' @param cases 
#' @param population 
#' @param expected.cases 
#' @param pop.upper.bound 
#' @param n.simulations 
#' @param alpha.level 
#' @param plot 
#'
#' @return
#' @export
#'
#' @examples
kulldorff <-
  function(geo, cases, population, expected.cases=NULL, pop.upper.bound, 
           n.simulations, alpha.level, plot=TRUE){
    
    #-------------------------------------------------------------------------------
    # Initialization 
    #-------------------------------------------------------------------------------
    # Determine likelihood type: binomial or poisson
    if(is.null(expected.cases)){
      type <- "binomial"
      denominator <- population
      expected.cases <- sum(cases) * (denominator/sum(denominator))
      # poisson case	
    }else{
      type <- "poisson"
      denominator <- expected.cases
    }
    
    # Get geographic information
    geo.results <- zones(geo, population, pop.upper.bound)
    nearest.neighbors <- geo.results$nearest.neighbors
    cluster.coords <- geo.results$cluster.coords
    n.zones <- nrow(cluster.coords)
    
    
    #-------------------------------------------------------------------------------
    # Observed statistic computation
    #-------------------------------------------------------------------------------
    lkhd <- computeAllLogLkhd(cases, denominator, nearest.neighbors, n.zones, type)
    
    # Get areas included in most likely cluster
    cluster.index <- which.max(lkhd)
    
    # cluster center and radial area
    center <- cluster.coords[cluster.index,1]
    end <- cluster.coords[cluster.index,2]
    
    # list of all areas included in cluster	
    cluster <- nearest.neighbors[[center]]
    cluster <- cluster[1:which(cluster == end)]
    
    
    #-------------------------------------------------------------------------------
    # Compute Monte Carlo randomized p-value
    #-------------------------------------------------------------------------------
    # Simulate cases under null hypothesis of no area effects i.e. conditioned on E
    perm <- rmultinom(n.simulations, round(sum(cases)), prob=denominator)
    
    # Ensure that # of cases is less than population size of each area. If not,
    # resample
    # for(i in 1:ncol(perm)){
    #   while(any(perm[,i] > population))
    #     perm[,i] <- rmultinom(1, round(sum(cases)), prob=denominator)
    # }
    
    
    # Compute simulated lambda's:  max log-lkhd in region
    sim.lambda <- kulldorffMC(perm, denominator, nearest.neighbors, n.zones, type)
    
    # Compute Monte Carlo p-value
    combined.lambda <- c(sim.lambda, max(lkhd))
    p.value <- 1-mean(combined.lambda < max(lkhd))
    
    # Plot histogram
    if(plot){
      hist(combined.lambda, 
           main="Monte Carlo Distribution of Lambda",
           xlab=expression(log(lambda)))
      abline(v=max(lkhd), col="red")
      legend("top",
             c(paste("Obs. log(Lambda) = ",round(max(lkhd),3),sep=""), 
               paste("p-value = ", round(p.value,log10(n.simulations + 1)),sep="")),
             lty=c(1, 1), 
             col=c("red","white"), 
             bty="n"
      )
    }
    
    
    #-------------------------------------------------------------------------------
    # Create Most Likely Cluster Object
    #-------------------------------------------------------------------------------
    most.likely.cluster = list(
      location.IDs.included = cluster,
      population = sum(population[cluster]),	
      number.of.cases = sum(cases[cluster]),
      expected.cases = sum(expected.cases[cluster]),
      SMR = sum(cases[cluster])/sum(expected.cases[cluster]),
      log.likelihood.ratio = lkhd[cluster.index],
      monte.carlo.rank = sum(combined.lambda >= lkhd[cluster.index]),
      p.value = p.value
    )
    
    
    #-------------------------------------------------------------------------------
    # Investigate Secondary Clusters
    #-------------------------------------------------------------------------------
    # running list of areas already covered to test overlap
    current.cluster <- cluster
    secondary.clusters <- NULL
    
    # go through log-likelihoods in decreasing order, skipping largest which 
    # corresponds to most likely cluster
    indices <- order(lkhd, decreasing=TRUE)
    
    for(i in 2:length(indices)){
      #---------------------------------------------------
      # Get areas included in new potential cluster
      new.cluster.index <- indices[i]
      new.center <- cluster.coords[new.cluster.index, 1]
      new.end <- cluster.coords[new.cluster.index, 2]
      
      new.cluster <- nearest.neighbors[[new.center]]
      new.cluster <- new.cluster[1:which(new.cluster == new.end)]	
      #---------------------------------------------------
      # if there is no overlap between existing clusters and the new potential cluster
      if(length(intersect(new.cluster,current.cluster)) == 0){		
        new.secondary.cluster <- list(
          location.IDs.included = new.cluster,
          population = sum(population[new.cluster]),
          number.of.cases = sum(cases[new.cluster]),
          expected.cases = sum(expected.cases[new.cluster]),
          SMR = sum(cases[new.cluster])/sum(expected.cases[new.cluster]),
          log.likelihood.ratio = lkhd[new.cluster.index],
          monte.carlo.rank = sum(combined.lambda >= lkhd[new.cluster.index]),
          p.value =  1 - mean(combined.lambda < lkhd[new.cluster.index])
        )
        
        # If no longer significant, end loop.
        if(new.secondary.cluster$p.value > alpha.level){
          break
        }
        # Otherwise add to existing clusters
        secondary.clusters <- c(secondary.clusters,list(new.secondary.cluster))
        current.cluster <- unique(c(current.cluster, new.cluster))
      }		
    }
    
    
    #-------------------------------------------------------------------------------
    # Output results
    #-------------------------------------------------------------------------------
    results <- list(
      most.likely.cluster=most.likely.cluster,
      secondary.clusters=secondary.clusters,
      type = type,
      log.lkhd=lkhd,
      simulated.log.lkhd=sim.lambda
    )
    return(results)
  }
