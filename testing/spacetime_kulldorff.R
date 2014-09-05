library(dplyr)
library(SpatialEpi)
library(maps)
library(maptools)


#-------------------------------------------------------------------------------
#
# Create expected counts and maps
#
#-------------------------------------------------------------------------------
#---------------------------------------------------------------
# Import data
#---------------------------------------------------------------
geo <- 
  read.table("./testing/nm.geo", header=FALSE, col.names=c("county", "lat", "long")) %>%
  tbl_df()

pop <- 
  read.table("./testing/nm.pop", header=FALSE, col.names=c("county", "year", "pop", "age", "gender")) %>%
  tbl_df() %>%
  mutate(year=year+1900) %>%
  arrange(county, year, age, gender) %>%
  select(county, year, age, gender, pop)

cases <- 
  read.table("./testing/nm.cas", header=FALSE, col.names=c("county", "cases", "year", "age", "gender")) %>%
  tbl_df() %>%
  group_by(county, year, age, gender) %>%
  summarize(cases=sum(cases)) %>%
  arrange(county, year, age, gender) %>%
  ungroup()


#---------------------------------------------------------------
# sp object
#---------------------------------------------------------------
nmTemp <- map('county','new.mexico',fill=TRUE,plot=FALSE)
nmIDs <- substr(nmTemp$names,1+nchar("new.mexico,"),nchar(nmTemp$names) )
nm <- map2SpatialPolygons(nmTemp,IDs=nmIDs,proj4string=CRS("+proj=longlat"))

nm.id <- NULL
for (i in 1:length(nm)) {
  nm.id <- c(nm.id, nm@polygons[[i]]@ID)
}
nm.id <- c("bernalillo", "catron", "chaves", "valencia", "colfax", "curry", 
          "debaca", "donaana", "eddy", "grant", "guadalupe", "harding", 
          "hidalgo", "lea", "lincoln", "losalamos", "luna", "mckinley", 
          "mora", "otero", "quay", "rioarriba", "roosevelt", "sanjuan", 
          "sanmiguel", "sandoval", "santafe", "sierra", "socorro", "taos", 
          "torrance", "union", "valencia")
nm <- unionSpatialPolygons(nm, nm.id)
plot(nm, axes=TRUE)
nm.id <- NULL
for (i in 1:length(nm)) {
  nm.id <- c(nm.id, nm@polygons[[i]]@ID)
}
text(coordinates(nm), labels=nm.id, cex=0.75)
text(geo[,3:2], labels=geo$county, cex=0.75, col="red")

tolower(as.character(geo$county)) == rownames(coordinates(nm))




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
ref.years <- data.frame(orig.year=case.years, ref.year=case.ref.years)
year.counts <- table(case.ref.years)

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
q <- mutate(q, q=cases/pop) %>% 
  select(-c(cases:pop))


#---------------------------------------------------------------
# Compute expected values E, y, pop for each year/county
#---------------------------------------------------------------
counts <- tbl_df(expand.grid(geo$county, unique(cases$year))) %>%
  mutate(county=Var1, year=Var2, E=0, pop=0) %>% 
  select(county:E)

# Add observed counts
cases.agg <- group_by(cases, county, year) %>% 
  summarise(y=sum(cases)) 
pop.agg <- group_by(pop, county, year) %>% 
  summarise(pop=sum(pop)) 

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
  
  counts[i, "pop"] <- filter(pop, county==counts[i, "county"], year==ref.year) %>%
    inner_join(q, by=c("age", "gender")) %>%
    summarize(pop=sum(pop)) 
}





#-------------------------------------------------------------------------------
#
# Compute SatScan
#
#-------------------------------------------------------------------------------
#---------------------------------------------------------------
# Compare E's with most likely cluster from SatScan
#---------------------------------------------------------------
cluster <- c("Torrance", "Bernalillo", "Valencia", "SantaFe", "Guadelupe", "Socorro", "Sandoval", "SanMiguel", "LosAlamos")
count <- filter(counts, county %in% cluster & 1985 <= year & year <= 1989) %>%
  summarize(y=sum(y), E=sum(E))
(count$y/count$E) / ((sum(counts$y)-count$y)/(sum(counts$E)-count$E))


#---------------------------------------------------------------
# Geographic zones.  Based on average population over the 19 years
#---------------------------------------------------------------
population <- group_by(counts, county) %>% 
  summarise(pop=mean(pop))
population <- population[,2]

new.geo <- latlong2grid(coordinates(nm)[,2:1])
geo.objects <- create_geo_objects(0.5, population, new.geo, nm)
cluster.list <- geo.objects$overlap$cluster.list
n.zones <- length(cluster.list)


#---------------------------------------------------------------
# Temporal zones
#---------------------------------------------------------------
years <- unique(cases$year)
n.year.span <- 7
time.windows <- as.vector(years, mode="list")

for(i in 2:n.year.span) {
  # define all n-tuples
  blah <- matrix(0, nrow=i, ncol=length(years)-i+1)
  for(j in 1:i) {
    blah[j, ] <- c(years[j]:years[length(years)-i+j])
  }
  
  # Append to list here
  temp.list <- vector(mode="list", length=ncol(blah))
  for (j in 1:ncol(blah))
    temp.list[[j]] <- blah[, j]
  
  time.windows <- c(time.windows, temp.list)
}


#---------------------------------------------------------------
# Observed statistic computation
#---------------------------------------------------------------
# log lkhd function
log.lkhd <- function (cz, nz, C, N) {  
  if(cz / nz <= (C - cz)/(N - nz)) {
    log.lkhd <- 0
  } else {
    log.lkhd <-
      cz * log((cz / nz)) + cz * log(((N - nz)/(C - cz))) + C * log(((C-cz)/(N-nz))) + C * log((N/ C))
  }  
  return(log.lkhd)
}

# compute log lkhd function over all possible zones
compute.lkhds <- function(counts){
  lkhd <- matrix(0, nrow=n.zones, ncol=length(time.windows))
  C <- sum(counts$y)
  N <- sum(counts$E)
  for (i in 1:n.zones) {
    for (j in 1:length(time.windows)) {
      areas <- geo$county[cluster.list[[i]]]
      year.min <- min(time.windows[[j]])
      year.max <- max(time.windows[[j]])
      
      count <- 
        filter(counts, county %in% areas & year.min <= year & year <= year.max) %>%
        summarize(y=sum(y), E=sum(E))
      lkhd[i, j] <- log.lkhd(count$y, count$E, C, N)
    }
  }
  return(lkhd)
}

# Monte Carlo simulation
n.sim <- 100
max.lkhds <- rep(0, n.sim)
sim.counts <- counts
for(i in 1:100){
  sim.counts$y <- rmultinom(1, sum(sim.counts$E), normalize(sim.counts$E))
  max.lkhds[i] <- max(compute.lkhds(sim.counts))
  if(i %% 10 == 0) print(i)
}

save.image(file="satscan.RData")


# Most likely cluster
lkhd <- compute.lkhds(counts)
max(lkhd)
which(lkhd == max(lkhd), arr.ind=TRUE)

# Secondary non-overlapping cluster
non.overlap <- 
  lapply(cluster.list, function(x){!any(x %in% mlc)}) %>%
  unlist()
second.max <- max(lkhd[non.overlap,])
which(lkhd == second.max, arr.ind=TRUE)

# Secondary Cluster
plot(nm, axes=TRUE)
plot(nm[cluster.list[[388]]], add=TRUE, col="red")
plot(nm[3], add=TRUE, col="blue")

# SatScan most likely
plot(nm[match(cluster, geo$county)], add=TRUE, col="blue")

