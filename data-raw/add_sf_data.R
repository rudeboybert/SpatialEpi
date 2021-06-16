library(sf)
library(plyr)
library(dplyr)
library(Hmisc)
library(usethis)
library(readr)
library(SpatialEpi)

# Adding pennLC Pennsylvania lung cancer data as shapefile
pennLC_sf <- read_csv("data-raw/pennLC_sf.csv")
pennLC_sf <- st_as_sf(pennLC$spatial.polygon)

# add all of the cases together by county
pennLC_df <- pennLC$data %>%
  group_by(county) %>%
  summarise(cases = sum(cases))

# add quotes to each county 
county <- Hmisc::Cs(adams, allegheny,
                    armstrong,
                    beaver,
                    bedford,
                    berks,
                    blair,
                    bradford,
                    bucks,
                    butler,
                    cambria,
                    cameron,
                    carbon,
                    centre,
                    chester,
                    clarion,
                    clearfield,
                    clinton,
                    columbia,
                    crawford,
                    cumberland,
                    dauphin,
                    delaware,
                    elk,
                    erie,
                    fayette,
                    forest,
                    franklin,
                    fulton,
                    greene,
                    huntingdon,
                    indiana,
                    jefferson,
                    juniata,
                    lackawanna,
                    lancaster,
                    lawrence,
                    lebanon,
                    lehigh,
                    luzerne,
                    lycoming,
                    mckean,
                    mercer,
                    mifflin,
                    monroe,
                    montgomery,
                    montour,
                    northampton,
                    northumberland,
                    perry,
                    philadelphia,
                    pike,
                    potter,
                    schuylkill,
                    snyder,
                    somerset,
                    sullivan,
                    susquehanna,
                    tioga,
                    union,
                    venango,
                    warren,
                    washington,
                    wayne,
                    westmoreland,
                    wyoming,
                    york)

# add county to sf object manually
pennLC_sf$county <- county

# now that we have a common value we can join
pennLC_sf <- left_join(pennLC_sf, pennLC_df)

# reorder columns so that geometry is first
pennLC_sf <- pennLC_sf[c(3,2,1)]


usethis::use_data(pennLC_sf,overwrite = TRUE)
