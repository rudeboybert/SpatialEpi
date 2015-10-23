# Save NYleukemia data as SpatialPolygonsDataFrame
load("/Users/rudeboybert/Dropbox/SpatialEpi/data/NYleukemia.rda")
sp <- NYleukemia$spatial.polygon
data <- NYleukemia$data %>% 
  rename(FIPS = censustract.FIPS) %>% 
  mutate(ID = 1:nrow(NYleukemia$data)) %>% 
  select(ID, FIPS, population, cases)

NYleukemia <- SpatialPolygonsDataFrame(sp, data)
devtools::use_data(NYleukemia, NYleukemia, overwrite=TRUE)


# Testing plotting function
options(stringsAsFactors=FALSE)


library(maptools)
shapefile_name <- system.file("extdata/sids/sids.shp", package="SpatialEpi")
nc_sids <- maptools::readShapePoly(
  fn = shapefile_name,
  ID = "FIPSNO",
  proj4string = CRS("+proj=longlat +ellps=clrk66")
)


library(rgdal)
shapefile_dir <- "/Library/Frameworks/R.framework/Versions/3.2/Resources/library/spdep/etc/shapes/"
nc_sids <- readOGR(
  dsn=shapefile_dir, 
  layer="sids", 
  p4s = "+proj=longlat +ellps=clrk66",
  verbose = FALSE
)



sp_obj <- nc_sids
variable <- "SID74"
x_label <- "longitude"
y_label <- "latitude"
title <- "North Carolina SIDS Deaths"
legend_label <- "SIDS Deaths"

plot <- map_variable(sp_obj, "SID74", x_label, y_label)

map_variable <- function(sp_obj, variable=NULL, x_label = NULL, y_label = NULL, 
                        title = NULL, legend_label = NULL){
  if(is.null(x_label)){
    x_label <- "eastings"
  }
  if(is.null(y_label)){
    y_label <- "northings"
  }
  if(is.null(legend_label)){
    legend_label <- variable
  }
  
  # Create data frame from SpatialPolygonsDataFrame object for ggplot'ing
  sp_obj@data$id <- rownames(sp_obj@data)
  map_points <- fortify(sp_obj, region="id")
  map_df <- inner_join(map_points, sp_obj@data, by="id")
  
  # Set up base plot
  map_ggplot <- 
    ggplot(map_df, aes(x=long, y=lat, group=group)) +
    coord_map() +
    theme_bw() +
    xlab(x_label) + 
    ylab(y_label) + 
    ggtitle(title) + 
    geom_polygon(fill="white") + 
    geom_path(col="black", size=0.5)
  
  # If specified, plot variable
  if(!is.null(variable)) {
    map_ggplot <- map_ggplot + 
      scale_fill_continuous(low="white", high="black", name=legend_label) +
      geom_polygon(aes_string(fill=variable)) + 
      geom_path(col="black", size=0.5)
  }

  # Return ggplot object
  return(map_ggplot)
}