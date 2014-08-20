library(dplyr)
library(SpatialEpi)


#-------------------------------------------------------------------------------
#
# Recreate Kulldorff spatial results
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

cluster <- c("Torrance", "Bernalillo", "Valencia", "SantaFe", "Guadelupe", "Socorro", "Sandoval", "SanMiguel", "LosAlamos")
filter(counts, county %in% cluster & 1985 <= year & year <= 1989) %>%
  summarize(y=sum(y), E=sum(E))






cases <- 
  mutate(cases, year=as.integer(substr(cases$year, 3, 4))) %>%
  group_by(county, age, gender) %>%
  summarise(cases=sum(cases)) 

pop <- group_by(pop, county, age, gender) %>%
  summarise(pop=sum(pop))

cases.full <- select(pop, county, age, gender) %>%
  left_join(cases) 
cases.full$cases <- ifelse(is.na(cases.full$cases), 0, cases.full$cases)

E <- expected(pop$pop, cases.full$cases, 18*2)
y <- group_by(cases.full, county) %>%
  summarise(cases=sum(cases)) %>%
  select(cases)
y <- y[,1]




filter(cases, county %in% counties & 1985 <= year & year <= 1989) %>%
  summarize(counts=sum(cases))

sum(E[geo$county %in% counties])/3
