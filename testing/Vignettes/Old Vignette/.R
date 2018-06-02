###################################################
### chunk number 1: 
###################################################
library(SpatialEpi)
Rcode <- system.file("doc","SpatialEpi.Rnw",package="SpatialEpi")
options(device.ask.default=FALSE)
Stangle(Rcode)


###################################################
### chunk number 2: 
###################################################
data(scotland)


###################################################
### chunk number 3: 
###################################################
data(scotland)
polygon <- scotland$polygon$polygon; nrepeats <- scotland$polygon$nrepeats
names <- scotland$data$county.names

spatial.polygon <- polygon2spatial.polygon(polygon,coordinate.system='+proj=utm',names,nrepeats)

par(mfrow=c(1,2))
plot(polygon,type='n',xlab="Eastings (km)",ylab="Northings (km)",main="Polygon File")
polygon(polygon)

plot(spatial.polygon,axes=TRUE)
title(xlab="Eastings (km)",ylab="Northings (km)",main="Spatial Polygon")
plot(spatial.polygon[23],add=TRUE,col="red")


###################################################
### chunk number 4: 
###################################################
par(mfrow=c(1,2))
plot(polygon,type='n',xlab="Eastings (km)",ylab="Northings (km)",main="Polygon File")
polygon(polygon)

plot(spatial.polygon,axes=TRUE)
title(xlab="Eastings (km)",ylab="Northings (km)",main="Spatial Polygon")
plot(spatial.polygon[23],add=TRUE,col="red")


###################################################
### chunk number 5: 
###################################################
library(maps)
county.map <- map('county',c('pennsylvania','vermont'),fill=TRUE,plot=FALSE)
county.names <- as.character(county.map$names)
county <- map2SpatialPolygons(county.map,IDs=county.names,proj4string=CRS("+proj=longlat"))
state.map <- map('state',c(),fill=TRUE,plot=FALSE)
state.names <- as.character(state.map$names)
state <- map2SpatialPolygons(state.map,IDs=state.names,proj4string=CRS("+proj=longlat"))
plot(county,axes=TRUE,border="red")
plot(state,add=TRUE,lwd=2)


###################################################
### chunk number 6: 
###################################################
plot(county,axes=TRUE,border="red")
plot(state,add=TRUE,lwd=2)


###################################################
### chunk number 7: 
###################################################
county.grid <- latlong2grid(county)
state.grid <-  latlong2grid(state)
plot(county.grid,axes=TRUE,border="red")
plot(state.grid,add=TRUE,lwd=2)


###################################################
### chunk number 8: 
###################################################
plot(county.grid,axes=TRUE,border="red")
plot(state.grid,add=TRUE,lwd=2)


###################################################
### chunk number 9: 
###################################################
coord <- rbind(
c(-73.7500, 45.4667),
c(-122.6042, 45.6605)
)
latlong2grid(coord)


###################################################
### chunk number 10: 
###################################################
data(scotland)
scotland.map <- scotland$spatial.polygon
y <- runif(nrow(scotland$data))
mapvariable(y, scotland.map)


###################################################
### chunk number 11: 
###################################################
file="scotland.pdf"
pdf(file=file)  
mapvariable(y, scotland.map, xlab="Eastings (km)", ylab="Northings (km)")
graphics.off()  
cat("\\includegraphics{", file, "}\n\n", sep="")   


###################################################
### chunk number 12: 
###################################################
data(pennLC)
penn.map <- pennLC$spatial.polygon
penn.map <- latlong2grid(penn.map)

population <- tapply(pennLC$data$population,pennLC$data$county,sum)
cases <- tapply(pennLC$data$cases,pennLC$data$county,sum) 
geo <- latlong2grid(pennLC$geo[,2:3])

incidence <- (cases/population)*1000

mapvariable(incidence,penn.map)
#par(new=TRUE)
#layout(matrix(c(1, 2), ncol = 2, nrow = 1), heights = c(0.3, 0.3), widths = c(0.1, 0.4))
#plot.new()
#plot(penn.map)
#plot(SpatialPoints(geo, penn.map@proj4string),pch="x",col="red",add=T)


###################################################
### chunk number 13: 
###################################################
file="pennLCIncidence.pdf"
pdf(file=file)  
mapvariable(incidence,penn.map)
graphics.off()  
cat("\\includegraphics{", file, "}\n\n", sep="")   


###################################################
### chunk number 14: 
###################################################
data(scotland)
scotland.map <- scotland$spatial.polygon
y <- scotland$data$cases
E <- scotland$data$expected
SMR <- y/E
mapvariable(SMR, scotland.map)


###################################################
### chunk number 15: 
###################################################
file="scotlandSMR.pdf"
pdf(file=file)  
mapvariable(SMR, scotland.map, xlab="Eastings (km)", ylab="Northings (km)")
graphics.off()  
cat("\\includegraphics{", file, "}\n\n", sep="")   


###################################################
### chunk number 16: 
###################################################
data(pennLC)
n.strata <- 16
population <- tapply(pennLC$data$population,pennLC$data$county,sum)
cases <- tapply(pennLC$data$cases,pennLC$data$county,sum)
expected.cases <- expected(pennLC$data$population,pennLC$data$cases,n.strata)


###################################################
### chunk number 17: 
###################################################
data(scotland)
data <- scotland$data
x <- data$AFF
Xmat <- cbind(x,x^2)
results <- eBayes(data$cases,data$expected,Xmat)
scotland.map <- scotland$spatial.polygon
mapvariable(results$RR, scotland.map)


###################################################
### chunk number 18: 
###################################################
file="scotlandEB.pdf"
pdf(file=file)  
mapvariable(results$RR, scotland.map)
graphics.off()  
cat("\\includegraphics{", file, "}\n\n", sep="")   


###################################################
### chunk number 19: 
###################################################
## Load Pennsylvania Lung Cancer Data
data <- pennLC$data

## Process geographical information and convert to grid
geo <- latlong2grid( pennLC$geo[,2:3] )

## Get aggregated counts of population and cases for each county
population <- tapply(data$population,data$county,sum)
cases <- tapply(data$cases,data$county,sum)

## Based on the 16 strata levels, computed expected numbers of disease
expected.cases <- expected(data$population, data$cases, 16)

## Set Parameters
pop.upper.bound <- 0.5
n.simulations <- 999
alpha.level <- 0.05
plot <- TRUE


###################################################
### chunk number 20: 
###################################################
## Kulldorff using Binomial likelihoods
binomial <- kulldorff(geo, cases, population, NULL, pop.upper.bound, n.simulations, alpha.level, plot)
cluster <- binomial$most.likely.cluster$location.IDs.included

## plot
plot(pennLC$spatial.polygon,axes=TRUE)
plot(pennLC$spatial.polygon[cluster],add=TRUE,col="red")
title("Most Likely Cluster")


###################################################
### chunk number 21: 
###################################################
plot(pennLC$spatial.polygon,axes=TRUE)
plot(pennLC$spatial.polygon[cluster],add=TRUE,col="red") 


###################################################
### chunk number 22: 
###################################################
## Kulldorff using Poisson likelihoods
poisson <- kulldorff(geo, cases, population, expected.cases, pop.upper.bound, n.simulations, alpha.level, plot)
cluster <- poisson$most.likely.cluster$location.IDs.included

## plot
plot(pennLC$spatial.polygon,axes=TRUE)
plot(pennLC$spatial.polygon[cluster],add=TRUE,col="red")
title("Most Likely Cluster Controlling for Strata")


###################################################
### chunk number 23: 
###################################################
plot(pennLC$spatial.polygon,axes=TRUE)
plot(pennLC$spatial.polygon[cluster],add=TRUE,col="red")
title("Most Likely Cluster")


###################################################
### chunk number 24: 
###################################################
k <- 1250
alpha.level <- 0.05

# not controlling for stratas
results <- besag.newell(geo, population, cases, expected.cases=NULL, k, alpha.level)

# controlling for stratas
results <- besag.newell(geo, population, cases, expected.cases, k, alpha.level)


###################################################
### chunk number 25: 
###################################################
library(INLA)
# map file taken from the graph file in INLA
g <- system.file("demodata/scotland.graph", package="INLA")
m <- numNeighbors(g)
adj <- listNeighbors(g)
# Create, simulate, estimate omega2 from sample
Q <- make.Q(m, adj,  omega.sq=1)
tmp.x1 <- sim.Q(Q)


###################################################
### chunk number 26: 
###################################################
file="simulated_u_scotland.pdf"
pdf(file=file)
data(scotland) 
mapvariable(tmp.x1, scotland$spatial.polygon)
graphics.off()  
cat("\\includegraphics{", file, "}\n\n", sep="")  


###################################################
### chunk number 27: 
###################################################
# Function to choose the parameters of a lognormal, by specifying two
# quantiles with associated probs.
param <- LogNormalPriorCh(theta1=1, theta2=3, prob1=0.5, prob2=0.95)
param


###################################################
### chunk number 28: 
###################################################
x<-rnorm(1000,param$mu,param$sigma)
plot( function(y) dnorm(y, mean(x), sd(x)), from=min(x), to=max(x), xlab=expression(beta[1]), ylab="Density")
title("Prior Distribution")


###################################################
### chunk number 29: 
###################################################
# degrees of freedom is 1
param1 <- GammaPriorCh(log(10), 0.975,1)
# degrees of freedom is 1
param2 <- GammaPriorCh(log(5), 0.975,2)


###################################################
### chunk number 30: 
###################################################
x<-rnorm(1000,param$mu,param$sigma)
# plot of the prior distribution
par(mfrow=c(1,2))
curve(dgamma(x,shape=param1$a,rate=param1$b),from=0,to=6,n=1000,xlab=expression(tau^{-1}),ylab="density",main="d=1")
curve(dgamma(x,shape=param2$a,rate=param2$b),from=0,to=6,n=1000,xlab=expression(tau^{-1}),ylab="density", main="d=2")


