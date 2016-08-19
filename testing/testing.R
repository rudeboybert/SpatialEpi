# TODO
# * Write better package description
# * Add URL to Journal of Statistical Software!
# * Add county names to NYleukemia
# * Cite SIDS data
# * Add Scotland data
# * NYleukemia, scotland source
# * mapproj, gpclib packages
# * Fix gpclibPermit()










# R CMD CHECK returns NOTE that there is "no visible binding
# for global variable" for variables used by dplyr/ggplot2, hence this hack fix.
# utils::globalVariables(c("lat", "long", "group", ".", "variable"))

#' Plot ggplot compliant map
#'
#' As described in
#' \url{https://github.com/hadley/ggplot2/wiki/plotting-polygon-shapefiles}
#'
#' @param sp_obj SpatialPolygonsDataFrame object
#' @param variable_name Variable within data frame \code{sp_obj@@data}
#' @param type string indicating type of plot. Must be either
#'   \code{"continuous"} for continuous scale, \code{"discrete"} for categorical
#'   variables, \code{"quantile"} for quantile-based bins, or
#'   \code{"equidistant"} for equally-sized bins.
#' @param n_int For non-continuous scales, number of bins.
#' @param x_label Label for eastings
#' @param lwd line width between areas
#' @param y_label Label for northings
#' @param title Plot title
#' @param legend_title label for legend. Defaults to \code{variable_name}.
#'
#' @import maptools
#'
#' @return a ggplot object
#' @export
map_variable <- function(sp_obj, variable_name = NULL, type = "continuous", n_int=5,
                         x_label = "eastings", y_label = "northings", lwd=0.5,
                         title = NULL, legend_title = NULL){

  if(is.null(legend_title)){
    legend_title <- variable_name
  }

  # Create variables to plot
  if(!is.null(variable_name)) {
    sp_obj@data <- sp_obj@data %>% mutate_(variable = variable_name)

    if(type == "continuous"){
      # Do nothing
    } else if(type == "discrete"){
      sp_obj$variable <- as.factor(sp_obj$variable)
    } else if(type == "quantile") {
      sp_obj$variable <- cut_number(sp_obj$variable, n=n_int)
    } else if (type == "equidistant") {
      sp_obj$variable <- cut_interval(sp_obj$variable, n=n_int)
    }
  }




  # Create data frame from SpatialPolygonsDataFrame object for ggplot'ing
  data("NYleukemia")
  sp_obj <- NYleukemia
  x_label <- "eastings"
  y_label <- "northings"
  lwd <- 0.5
  title <- "blah"

  library(maptools)
  library(dplyr)

  sp_obj$id <- rownames(sp_obj@data)
  map_points <- fortify(sp_obj, region="id")
  map_df <- inner_join(map_points, sp_obj@data, by="id")



  # Set up base plot
  ggplot(map_df, aes(x=long, y=lat, group=group)) +
    coord_map() +
    #theme_bw() +
    xlab(x_label) +
    ylab(y_label) +
    ggtitle(title) +
    geom_polygon(aes(fill=population)) +
    scale_fill_continuous(low="white", high="black") +
    geom_path(col="black", size=lwd)





  # Depending on plot type, fill polygons and set legend
  if(!is.null(variable_name)){
    map_ggplot <- map_ggplot + geom_polygon(aes(fill=variable))

    if(type == "continuous") {
      map_ggplot <- map_ggplot + scale_fill_continuous(low="white", high="black", name=legend_title)
    } else if (type == "discrete"){
      n_col <- nlevels(sp_obj$variable)
      if (n_col == 1) {
        colours <- "black"
      } else {
        colours <- brewer.pal(name="Greys", n=n_col+1)[-1]
      }
      names(colours) <- rev(levels(sp_obj$variable))
      map_ggplot <- map_ggplot + scale_fill_manual(values=colours, name=legend_title)
    } else if (type %in% c("quantile", "equidistant", "jenks")) {
      map_ggplot <- map_ggplot + scale_fill_brewer(palette="Greys", name=legend_title)
    }
  }

  # Trace outlines of polygons
  map_ggplot <- map_ggplot + geom_path(col="black", size=lwd)

  return(map_ggplot)
}










# Save scotland data as SpatialPolygonsDataFrame
library(rgdal)
shapefile_dir <- "./inst/extdata/scotland/"
scotland <- rgdal::readOGR(dsn=shapefile_dir, layer="scot", p4s = "+proj=longlat +ellps=clrk66",
                           verbose = FALSE
)

data <- read.table("./inst/extdata/scotland/scotland.dat", skip=1)
colnames(data) <- c("District", "Observed", "Expected", "AFF", "Latitude", "Longitude")

scotland@data <- inner_join(scotland@data, data, by=c("ID"="District"))
scotland$SMR <- scotland$Observed/scotland$Expected
map_variable(scotland, "SMR", type="jenks")

library(stringr)
load("/Users/aykim/Dropbox/SpatialEpi/data/scotland.rda")
sp <- scotland$spatial.polygon

getSpPPolygonsIDSlots(sp)
names <- sapply(slot(sp, "polygons"), function(x) slot(x, "ID")) %>%
  str_replace_all("\\.", "_")

for(i in 1:length(sp)){
  slot(slot(sp, "polygons")[[i]], "ID") <- as.character(i)
}

data <- scotland$data %>%
  mutate(county_names=str_replace_all(county.names, "\\.", "_")) %>%
  select(county_names, cases, expected, AFF)

scotland <- SpatialPolygonsDataFrame(sp, data)
scotland$SMR <- scotland$cases/scotland$expected
map_variable(scotland, "SMR", type="equidistant")
# devtools::use_data(scotland, scotland, overwrite=TRUE)










# Save NYleukemia data as SpatialPolygonsDataFrame
library(rgeos)
load("/Users/rudeboybert/Dropbox/SpatialEpi/data/NYleukemia.rda")
sp <- NYleukemia$spatial.polygon
sp <- gSimplify(sp, tol = 0.00001)
sp <- gBuffer(sp, byid=TRUE, width=0)
# this is a well known R / GEOS hack (usually combined with the above) to
# deal with "bad" polygons
# http://gis.stackexchange.com/questions/163445/
data <- NYleukemia$data %>%
  rename(FIPS = censustract.FIPS) %>%
  mutate(ID = 1:nrow(NYleukemia$data)) %>%
  select(ID, FIPS, population, cases)
NYleukemia <- SpatialPolygonsDataFrame(sp, data)
# devtools::use_data(NYleukemia, NYleukemia, overwrite=TRUE)

library(rgdal)
shapefile_dir <- "/Users/rudeboybert/Google Drive/Documents/Research/SC/Datasets/NY/"
NYleukemia <- rgdal::readOGR(dsn=shapefile_dir, layer="NY", verbose = FALSE)

devtools::use_data(NYleukemia)










# Two different ways of loading in shapefiles
options(stringsAsFactors=FALSE)
library(maptools)
shapefile_name <- system.file("extdata/sids/sids.shp", package="SpatialEpi")
nc_sids <- maptools::readShapePoly(fn = shapefile_name, ID = "FIPSNO",
  proj4string = CRS("+proj=longlat +ellps=clrk66"))

library(rgdal)
shapefile_dir <- "/Library/Frameworks/R.framework/Versions/3.2/Resources/library/spdep/etc/shapes/"
nc_sids <- rgdal::readOGR(dsn=shapefile_dir, layer="sids", p4s = "+proj=longlat +ellps=clrk66",
  verbose = FALSE)

library(rgdal)
shapefile_dir <- "/Users/aykim/Downloads/"
scotland <- rgdal::readOGR(dsn=shapefile_dir, layer="scot", p4s = "+proj=longlat +ellps=clrk66",
                          verbose = FALSE)




