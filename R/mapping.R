# R CMD CHECK returns NOTE that there is "no visible binding 
# for global variable" for variables used by dplyr/ggplot2, hence this hack fix.
utils::globalVariables(c("lat", "long", "group", "."))


#' Plot ggplot compliant map
#' 
#' As described in \url{https://github.com/hadley/ggplot2/wiki/plotting-polygon-shapefiles}
#'
#' @param sp_obj SpatialPolygonsDataFrame object
#' @param variable Variable within data frame \code{sp_obj@@data}
#' @param x_label Label for eastings
#' @param y_label Label for northings
#' @param title Plot title
#' @param legend_label label for legend. Defaults to variable name.
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' shapefile_name <- system.file("extdata/sids/sids.shp", package="SpatialEpi")
#' proj4string <- sp::CRS("+proj=longlat +ellps=clrk66")
#' nc_sids <- maptools::readShapePoly(
#'    fn = shapefile_name,
#'    ID = "FIPSNO",
#'    proj4string = proj4string
#'    )
#' map_variable(nc_sids)
#' names(nc_sids)
#' nc_sids$RATE74 <- nc_sids$SID74/nc_sids$BIR74
#' nc_map <- map_variable(nc_sids, "RATE74")
#' nc_map
map_variable <- function(sp_obj, variable = NULL, x_label = "eastings", 
                         y_label = "northings", title = NULL, 
                         legend_label = NULL){
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