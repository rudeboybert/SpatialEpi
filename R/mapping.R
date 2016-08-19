#' Title
#'
#' @param y 
#' @param spatial.polygon 
#' @param ncut 
#' @param nlevels 
#' @param lower 
#' @param upper 
#' @param main 
#' @param xlab 
#' @param ylab 
#'
#' @return
#' @export
#'
#' @examples
mapvariable <- function(y, spatial.polygon, ncut=1000, nlevels=10, lower=NULL, 
                        upper=NULL, main=NULL, xlab=NULL, ylab=NULL){
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



#' Title
#'
#' @param values 
#' @param map 
#' @param log 
#' @param nclr 
#' @param include.legend 
#' @param lwd 
#' @param round 
#' @param brks 
#' @param legend 
#' @param location 
#' @param rev 
#' @param leg.cex 
#' 
#' @importFrom maptools leglabs
#'
#' @return
#' @export
#'
#' @examples
#' 
plotmap <- function(values, map, log=FALSE, nclr=7, include.legend=TRUE, 
                    lwd=0.5, round=3, brks=NULL, legend=NULL, 
                    location='topright', rev=FALSE, leg.cex=1){
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
    legend <- maptools::leglabs(signif(brks,digits=round))
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
  
  
  sp::plot(map, axes=TRUE, col=colcode, lwd=lwd)
  if(include.legend) {
    legend(location, legend=legend, fill=plotclr, bty="n", cex=leg.cex)
  }
  
}
