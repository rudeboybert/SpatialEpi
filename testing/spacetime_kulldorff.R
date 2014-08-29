library(dplyr)
library(SpatialEpi)


#-------------------------------------------------------------------------------
#
# Import data, compute expected and obseved counts, get map
#
#-------------------------------------------------------------------------------
#---------------------------------------------------------------
# Import data
#---------------------------------------------------------------
geo <- 
  read.table("./testing/nm.geo", header=FALSE, col.names=c("county", "lat", "long")) %>%
  tbl_df()
levels(geo$county)[10] <- "Guadalupe"
levels(geo$county) <- tolower(levels(geo$county))

pop <- 
  read.table("./testing/nm.pop", header=FALSE, col.names=c("county", "year", "pop", "age", "gender")) %>%
  tbl_df() %>%
  mutate(year=year+1900) %>%
  # mutate(county=tolower(county)) %>%
  arrange(county, year, age, gender) %>%
  select(county, year, age, gender, pop)
levels(pop$county)[10] <- "Guadalupe"
levels(pop$county) <- tolower(levels(pop$county))

cases <- 
  read.table("./testing/nm.cas", header=FALSE, col.names=c("county", "cases", "year", "age", "gender")) %>%
  tbl_df() %>%
  # mutate(county=tolower(county)) %>% 
  group_by(county, year, age, gender) %>%
  summarize(cases=sum(cases)) %>%
  arrange(county, year, age, gender) %>%
  ungroup()
levels(cases$county)[10] <- "Guadalupe"
levels(cases$county) <- tolower(levels(cases$county))


#---------------------------------------------------------------
# Associate population years to case years
#---------------------------------------------------------------
pop.years <- unique(pop$year)
case.years <- unique(cases$year)
case.ref.years <- rep(0, length(case.years))

# Find closest year in population data set
for (i in 1:length(case.ref.years)) {
  closest.pop.year <- which.min(abs(pop.years-case.years[i]))
  case.ref.years[i] <- pop.years[closest.pop.year]
}

# Number of times to use each population count
year.counts <- table(case.ref.years)

# Match case year with pop year
ref.years <- data.frame(orig.year=case.years, ref.year=case.ref.years)
get.year <- function(orig.year, ref.years) {
  return(ref.years$ref.year[which(ref.years$orig.year==orig.year)])
}


#---------------------------------------------------------------
# Compute reference rates q
#---------------------------------------------------------------
q <- group_by(cases, age, gender) %>% 
  summarise(cases=sum(cases)) %>% 
  mutate(pop=0) 

for(i in 1:nrow(q)) {
  pop.sub <- group_by(pop, year, age, gender) %>% 
    summarise(pop=sum(pop)) %>%
    filter(age==q[i, "age"], gender==q[i, "gender"]) 
  
  # Associate years
  q$pop[i] <- sum(rep(pop.sub$pop, times=year.counts))
}
q <- mutate(q, q=cases/pop) %>% select(-c(cases:pop))


#---------------------------------------------------------------
# Compute expected values E
#---------------------------------------------------------------
counts <- tbl_df(expand.grid(geo$county, unique(cases$year))) %>%
  mutate(county=Var1, year=Var2, E=0) %>% select(county:E)

# Add observed counts
cases.agg <- group_by(cases, county, year) %>% 
  summarise(y=sum(cases)) 
counts <- left_join(counts, cases.agg, by=c("county", "year")) %>%
  mutate(y = ifelse(is.na(y), 0, y))


for (i in 1:nrow(counts)) {
  ref.year <- get.year(counts[i, "year"], ref.years)  
  counts[i, "E"] <- filter(pop, county==counts[i, "county"], year==ref.year) %>%
    inner_join(q, by=c("age", "gender")) %>%
    mutate(E=pop*q) %>%
    summarize(E=sum(E)) %>%
    ungroup() %>%
    select(E)
}


#---------------------------------------------------------------
# Load sp map
#---------------------------------------------------------------
library(maps)
library(maptools)

# Load Geographic Centroid Data + Map
nmTemp <- map('county','new.mexico',fill=TRUE,plot=FALSE)
nmIDs <- substr(nmTemp$names,1+nchar("new.mexico,"),nchar(nmTemp$names) )
nm <- map2SpatialPolygons(nmTemp,IDs=nmIDs,proj4string=CRS("+proj=longlat"))

# Get current county names
nm.id <- NULL
for (i in 1:length(nm)) {
  nm.id <- c(nm.id, nm@polygons[[i]]@ID)
}
dput(nm.id)

# Reassign cibola county to valenica
nm.id <- c("bernalillo", "catron", "chaves", "valencia", "colfax", "curry", 
           "de baca", "dona ana", "eddy", "grant", "guadalupe", "harding", 
           "hidalgo", "lea", "lincoln", "los alamos", "luna", "mckinley", 
           "mora", "otero", "quay", "rio arriba", "roosevelt", "san juan", 
           "san miguel", "sandoval", "santa fe", "sierra", "socorro", "taos", 
           "torrance", "union", "valencia")
nm <- unionSpatialPolygons(nm, nm.id)
plot(nm, axes=TRUE, border="black", lwd=2)

# Get county names
nm.id <- NULL
for (i in 1:length(nm)) {
  nm.id <- c(nm.id, nm@polygons[[i]]@ID)
}
nm.id <- str_replace_all(nm.id, " ", "")

# Re-arrange sp objects
index <- match(geo$county, nm.id)
nm <- unionSpatialPolygons(nm, nm.id[index])
plot(nm, axes=TRUE, border="black", lwd=2)
text(coordinates(nm), nm.id)

values <- group_by(counts, county) %>% summarise(y=sum(y), E=sum(E))
values[11, "y"] <- 1
plotmap(log(values$y/values$E), nm)



#-------------------------------------------------------------------------------
#
# Import data, compute expected and obseved counts, get map
#
#-------------------------------------------------------------------------------
#---------------------------------------------------------------
# Compare with most likely cluster from SatScan
#---------------------------------------------------------------
cluster <- tolower(c("Torrance", "Bernalillo", "Valencia", "SantaFe", "Guadalupe", "Socorro", "Sandoval", "SanMiguel", "LosAlamos"))
count <- filter(counts, county %in% cluster & 1985 <= year & year <= 1989) %>%
  summarize(y=sum(y), E=sum(E))

cz <- count$y
nz <- count$E
C <- sum(counts$y)
N <- sum(counts$E)
(count$y/count$E) / ((sum(counts$y)-count$y)/(sum(counts$E)-count$E))

log.lkhd(cz, nz, C, N)


#---------------------------------------------------------------
# Set up SatScan
#---------------------------------------------------------------
log.lkhd <- function (cz, nz, C, N) {  
  if(cz / nz <= (C - cz)/(N - nz)) {
    log.lkhd <- 0
  } else {
    log.lkhd <- cz * log(cz / nz) + 
      cz * log((N - nz)/( C - cz))  +
      C * log((C-cz)/(N-nz)) +
      C * log( N/ C )
  }  
  return(log.lkhd)
}


# Get geographic zones
population <- group_by(pop, county) %>% summarise(pop=round(sum(pop)/3)) %>% select(pop) %>% as.data.frame()
population <- population[,1]
geo.results <- zones(as.data.frame(geo[,2:3]), population, 0.50)
geo.objects <- create_geo_objects(0.5, population, as.data.frame(geo[,2:3]), nm)
cluster.list <- geo.objects$overlap$cluster.list
n.zones <- length(cluster.list)


# Get temporal zones
years <- unique(cases$year)
lengths <- c(1:5)

windows <- as.vector(years, mode="list")
# two years
for(i in 2:5) {
  blah <- matrix(0, nrow=i, ncol=length(years)-i+1)
  for(j in 1:i) {
    blah[j, ] <- c(years[j]:years[length(years)-i+j])
  }
  temp <- vector(mode="list", length=ncol(blah))
  for (j in 1:ncol(blah))
    temp[[j]] <- blah[, j]
  
  windows <- c(windows, temp)
}


#---------------------------------------------------------------
# Observed statistic computation
#---------------------------------------------------------------
lkhd <- matrix(0, nrow=n.zones, ncol=length(windows))
C <- sum(counts$y)
N <- sum(counts$E)
for (i in 1:n.zones)
  for (j in 1:length(windows)) {
    areas <- geo$county[cluster.list[[i]]]
    year.min <- min(windows[[j]])
    year.max <- max(windows[[j]])
    count <- filter(counts, county %in% areas & year.min <= year & year <= year.max) %>%
      summarize(y=sum(y), E=sum(E))    
    cz <- count$y
    nz <- count$E
    lkhd[i,j] <- log.lkhd(cz, nz, C, N)
  }

which(lkhd == max(lkhd), arr.ind=TRUE)






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

# Compute simulated lambda's:  max log-lkhd in region
sim.lambda <- kulldorffMC(perm, denominator, nearest.neighbors, n.zones, type)

# Compute Monte Carlo p-value
combined.lambda <- c(sim.lambda, max(lkhd))
p.value <- 1-mean(combined.lambda < max(lkhd))




