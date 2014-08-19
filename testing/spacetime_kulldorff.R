library(dplyr)
library(SpatialEpi)

#
# Recreate Kulldorff spatial results
#
cases <- tbl_df(read.table("./testing/nm.cas", header=FALSE))
names(cases) <- c("county", "cases", "year", "age", "gender")
cases <- arrange(cases, county, year, age, gender)

pop <- tbl_df(read.table("./testing/nm.pop", header=FALSE))
names(pop) <- c("county", "year", "pop", "age", "gender")
pop <- arrange(pop, county, year, age, gender)

geo <- tbl_df(read.table("./testing/nm.geo", header=FALSE))
names(geo) <- c("county", "lat", "long")



# Compute expected numbers
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




counties <- c("Torrance", "Bernalillo", "Valencia", "SantaFe", "Guadelupe", "Socorro", "Sandoval", "SanMiguel", "LosAlamos")
filter(cases, county %in% counties & 1985 <= year & year <= 1989) %>%
  summarize(counts=sum(cases))

sum(E[geo$county %in% counties])/3
