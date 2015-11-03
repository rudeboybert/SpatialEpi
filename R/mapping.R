# R CMD CHECK returns NOTE that there is "no visible binding 
# for global variable" for variables used by dplyr/ggplot2, hence this hack fix.
utils::globalVariables(c("lat", "long", "group", ".", "variable"))


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
#' @importFrom RColorBrewer brewer.pal
#' @importFrom gpclib gpclibPermit
#'
#' @return a ggplot object
#' @export
#' 
#' @examples
#' # Read in shapefile
#' shapefile_name <- system.file("extdata/sids/sids.shp", package="SpatialEpi")
#' proj4string <- sp::CRS("+proj=longlat +ellps=clrk66")
#' nc_sids <- maptools::readShapePoly(
#'    fn = shapefile_name, ID = "FIPSNO", proj4string = proj4string
#'    )
#'    
#' # Blank map
#' map_variable(nc_sids)
#' 
#' # Map SIDS rate
#' nc_sids$RATE74 <- nc_sids$SID74/nc_sids$BIR74
#' nc_map <- map_variable(nc_sids, "RATE74")
#' nc_map
map_variable <- function(sp_obj, variable_name = NULL, type = "continuous", n_int=5,
                         x_label = "eastings", y_label = "northings", lwd=0.5,
                         title = NULL, legend_title = NULL){
  gpclibPermit()
  
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
    } else if (type == "jenks"){
      breaks <- classIntervals(sp_obj$variable, n_int, style="jenks")$brks
      sp_obj$variable <- cut(sp_obj$variable, breaks, include.lowest = TRUE)
    }
  }
  
  # Create data frame from SpatialPolygonsDataFrame object for ggplot'ing
  sp_obj$id <- rownames(sp_obj@data)
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
    geom_polygon(fill="white")
  
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