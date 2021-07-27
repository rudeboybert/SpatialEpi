globalVariables(c(
  "palette","gray","par","layout","axis","rect"
))

#' Plot Levels of a Variable in a Colour-Coded Map
#' 
#' @description Plot levels of a variable in a colour-coded map along with a legend.
#' 
#' @param y variable to plot
#' @param spatial.polygon an object of class SpatialPolygons (See [SpatialPolygons-class][sp::SpatialPolygons-class])
#' @param ncut number of cuts in colour levels to plot
#' @param nlevels number of levels to include in legend
#' @param lower lower bound of levels
#' @param upper upper bound of levels
#' @param main an overall title for the plot
#' @param xlab a title for the x axis
#' @param ylab a title for the y axis
#'
#' @references Bivand, R. S., Pebesma E. J., and Gomez-Rubio V. (2008) *Applied Spatial Data Analysis with R*.  Springer Series in Statistics.
#' E. J. Pebesma and R. S. Bivand. (2005) Classes and methods for spatial data in R. *R News*, **5**, 9--13. 
#' 
#' @author Jon Wakefield, Nicky Best, Sebastien Haneuse, and Albert Y. Kim
#' 
#' @return A map colour-coded to indicate the different levels of `y`
#' 
#' @export
#' @import sp
#'
#' @examples
#' data(scotland)
#' map <- scotland$spatial.polygon
#' y <- scotland$data$cases
#' E <- scotland$data$expected
#' SMR <- y/E
#' mapvariable(SMR,map,main="Scotland",xlab="Eastings (km)",ylab="Northings (km)")
mapvariable <-
function(y, spatial.polygon, ncut=1000, nlevels=10, lower=NULL, upper=NULL, 
         main=NULL, xlab=NULL, ylab=NULL){

#-------------------------------------------------------------------------------
# Create id indicators for coloring scheme
#-------------------------------------------------------------------------------
if (is.null(lower)) 
	lower <- min(y) #- 0.0001
if (is.null(upper)) 
	upper <- max(y) #+ 0.0001

id <- cut(y, breaks=seq(from=lower, to=upper, length=(ncut+1)))
id <- as.numeric(id)
id[is.na(id)] <- 0
id <- id + 1

# Set colours to grey scale
palette(gray(seq(1,0, len=(ncut+1))))


#-------------------------------------------------------------------------------
# Make the scale of the two axes the same
#-------------------------------------------------------------------------------
xrnge <- spatial.polygon@bbox[1,]
yrnge <- spatial.polygon@bbox[2,]	

xd <- xrnge[2] - xrnge[1]
yd <- yrnge[2] - yrnge[1]

if(xd > yd){		
	xplot <- xrnge
	yplot <- NULL
	yplot[1] <- ((yrnge[2] + yrnge[1])/2) - xd/2
	yplot[2] <- ((yrnge[2] + yrnge[1])/2) + xd/2
}

if(xd <= yd){
	yplot <- yrnge
	xplot <- NULL
	xplot[1] <- ((xrnge[2] + xrnge[1])/2) - yd/2
	xplot[2] <- ((xrnge[2] + xrnge[1])/2) + yd/2
}


#-------------------------------------------------------------------------------
# Plots
#-------------------------------------------------------------------------------
def.par <- par(no.readonly = TRUE)
layout(matrix(c(1,2),ncol=2,nrow=1), heights=c(.3,.3), widths=c(.4,.1))

# plot variable
plot(spatial.polygon, axes=TRUE, col=id)

# plot title and axis labels
if(!is.null(main)){
	title(main=main)
}
if(!is.null(xlab)){
	title(xlab=xlab)
}
if(!is.null(ylab)){
	title(ylab=ylab)
}

# plot legend
plot(c(0,1), c(0,1), type="n", axes=FALSE, xlab="", ylab="")
xlims <- rep(0, nlevels)
atpts <- rep(0, nlevels)
for(i in 1:nlevels){
	# xlims[i] <- format(min(x)+(i-1)*(max(x)-min(x))/(nlevels-1),digits=2)
	xlims[i] <- format(lower+(i-1)*(upper-lower)/(nlevels-1),digits=2)
	atpts[i] <- (i-1)/(nlevels-1)
}

axis(2, at=c(atpts[1:nlevels]), labels=c(xlims[1:nlevels]))
yb <- seq(0, (nlevels-2)/(nlevels-1), 1/(nlevels-1))
yt <- seq(1/(nlevels-1), 1, 1/(nlevels-1))
xl <- rep(0, nlevels-1)
xr <- rep(1, nlevels-1)
gr <- seq(0, 1, 1/nlevels)
gr <- max(gr) - gr
rect(xleft=xl, ybottom=yb, xright=xr, ytop=yt, col=gray(gr), border=TRUE)

# Reset colours to default
palette("default")
par(def.par)
}
