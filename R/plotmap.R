globalVariables(c(
  "grey","leglabs"
))
#' Plot Levels of a Variable in a Colour-Coded Map
#'
#' @description Plot levels of a variable in a colour-coded map.
#' @param values variable to plot
#' @param map an object of class SpatialPolygons (See [SpatialPolygons-class][sp::SpatialPolygons-class])
#' @param log boolean of whether to plot values on log scale
#' @param nclr number of colour-levels to use
#' @param include.legend boolean of whether to include legend
#' @param lwd line width of borders of areas
#' @param round number of digits to round to in legend
#' @param brks if desired, pre-specified breaks for legend
#' @param legend if desired, a pre-specified legend
#' @param location location of legend
#' @param rev boolean of whether to reverse colour scheme (darker colours for smaller values)
#'
#' @author Albert Y. Kim
#' 
#' @return
#' A map colour-coded to indicate the different levels of `values`.
#' @export
#'
#' @examples
#' ## Load data
#' data(scotland)
#' map <- scotland$spatial.polygon
#' y <- scotland$data$cases
#' E <- scotland$data$expected
#' SMR <- y/E
#' ## Plot SMR
#' plotmap(SMR, map, nclr=9, location="topleft")
plotmap <-
function(values, map, log=FALSE, nclr=7, include.legend=TRUE, lwd=0.5, round=3,
         brks=NULL, legend=NULL, location='topright', rev=FALSE){

# create colors, each based on quantiles of data
plotclr <- grey(1-seq(0, 1, by=1/(nclr-1)))

nclr <- nclr + 1
if(is.null(brks)){
  if(log){
    brks <- exp(
      seq(from=min(log(values)), to=max(log(values)), length.out=nclr)
      )
  }else{
    brks <- seq(from=min(values), to=max(values), length.out=nclr)
  }
}
nclr <- nclr - 1


# Assign colors based on intervals
if(rev){
  plotclr <- rev(plotclr)
}
colornum <- findInterval(values, brks, all.inside=T);
colcode <- plotclr[colornum]


# Modify legend
if(is.null(legend)){
  legend <- leglabs(signif(brks,digits=round))
  # First element of legend
  legend[1] <- paste(
    signif(min(values), digits=round), "-", 
    substr(legend[1], 7, nchar(legend[1]))
  )
  # Last element of legend
  legend[nclr]<- paste(
    substr(legend[nclr], 6, nchar(legend[nclr])), "-",
    signif(max(values), digits=round)
  )
}


plot(map, axes=TRUE, col=colcode, lwd=lwd)
if(include.legend) {
  legend(location, legend=legend, fill=plotclr, bty="n")
}

}
