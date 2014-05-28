`expected` <-
function(population, cases, n.strata){
	expected.cases <- .Call("expected", population, cases, n.strata, PACKAGE = "SpatialEpi")$E
	expected.cases
}

