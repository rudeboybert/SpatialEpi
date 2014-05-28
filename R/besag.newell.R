`besag.newell` <- function(geo, population, cases, expected.cases=NULL, k, alpha.level){

#-------------------------------------------------------------------------
# Initialization 
#-------------------------------------------------------------------------
#---------------------------------------------------
# If no expected.cases provided, set them

# if there are no expected counts
if(is.null(expected.cases)){
	p <- sum(cases)/sum(population)
	expected.cases <- population*p
}

#---------------------------------------------------
# geographical information computation
geo.results <- zones(geo, population, 1)

# list of all nearest neighbors for each area
nearest.neighbors <- geo.results$nearest.neighbors
# interpoint distance matrix
distance <- geo.results$dist





#-------------------------------------------------------------------------
# Observed statistic computation
#-------------------------------------------------------------------------
results <- .Call("besagNewell", as.double(cases), as.double(expected.cases), nearest.neighbors, as.integer(k), PACKAGE="SpatialEpi") 





#-------------------------------------------------------------------------
# Process results
#-------------------------------------------------------------------------
#---------------------------------------------------
# significant areas
p.values <- results$observed.p.values	# observed p.values for each ares
m.values <- results$observed.m.values	# observed number of neighbors needed to observe k cases
k.values <- results$observed.k.values	# actual observed number of cases

# pick out areas that were significant and order them by p-value
signif.indices <- order(p.values)[1:sum(p.values <= alpha.level)]

# order remaining values
signif.p.values <- p.values[signif.indices]
signif.m.values <- m.values[signif.indices]
signif.k.values <- k.values[signif.indices]


#---------------------------------------------------
# Create object to output

# If none are significant, return NULL
if( length(signif.indices) == 0){
	clusters <- NULL
}else{
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





#-------------------------------------------------------------------------
# Output results
#-------------------------------------------------------------------------
results <- list(
	clusters=clusters,
	p.values=p.values,
	m.values=m.values,
	observed.k.values=k.values
)	
return(results)
}
	