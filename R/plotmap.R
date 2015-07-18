plotmap <-
function(values, map, log=FALSE, nclr=7, include.legend=TRUE, lwd=0.5, round=3,
         brks=NULL, legend=NULL, location='topright', rev=FALSE, leg.cex){

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
  legend <- paste(
    signif(brks, digits=round)[1:nclr],
    signif(brks, digits=round)[2:(nclr+1)],
    sep=" - "
  )
}


plot(map, axes=TRUE, col=colcode, lwd=lwd)
if(include.legend) {
  legend(location, legend = legend, fill = plotclr, bty = "n", cex=leg.cex, 
         pt.cex = 1)}

}
