ranking <-
function(variable, cluster.list){
	indices <- order(variable, decreasing = TRUE)
	current.cluster <- cluster.list[[indices[1]]]
	nonoverlap <- 1
	
	for (i in 2:length(indices)) {
		new.cluster <- cluster.list[[indices[i]]]
	
		if(!any(is.element(current.cluster,new.cluster))){
			current.cluster <- union(current.cluster,new.cluster)
			nonoverlap <- c(nonoverlap,i)
			next			
		}
	}	
	output <-	list( nonoverlap, indices[ nonoverlap ] )
	attributes(output)$names <- c("ranking","indices")
	output
}

