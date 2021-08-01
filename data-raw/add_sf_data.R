library(sf)
library(dplyr)
library(Hmisc)
library(usethis)
library(SpatialEpi)

# Adding pennLC Pennsylvania lung cancer data as shapefile
pennLC_sf <- st_as_sf(pennLC$spatial.polygon)

#load pennLC$data 
pennLC_df <- pennLC$data 

# add quotes to each county 
county <- Hmisc::Cs(adams, 
                    allegheny,
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
pennLC_sf <- left_join(pennLC_sf, pennLC_df, by = "county")

# Adding smoking variable
pennLC_sf <- left_join(pennLC_sf, tibble(pennLC$smoking), by = "county")

# Set projection: 
# https://stackoverflow.com/questions/61286108/error-in-cpl-transformx-crs-aoi-pipeline-reverse-ogrcreatecoordinatetrans
pennLC_sf <- st_transform(pennLC_sf, 4326)

usethis::use_data(pennLC_sf,overwrite = TRUE)

