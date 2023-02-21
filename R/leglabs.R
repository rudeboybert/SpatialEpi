#' Make legend labels
#'
#' leglabs makes character strings from the same break points. This function was copied from the soon-to-be 
#' deprecated `maptools` package with permission from author Roger Bivand
#'
#' @param vec vector of break values
#' @param under character value for under
#' @param over character value for over
#' @param between character value for between
#' @param reverse flag to reverse order of values, you will also need to reorder colours, see example
#'
#' @author Roger Bivand, Nick Bearman, Nicholas Lewin-Koh 
leglabs <- function(vec, under = "under", over = "over", between = "-", 
                    reverse = FALSE) {
  x <- vec
  lx <- length(x)
  if (lx < 3) 
    stop("vector too short")
  if (reverse) {
    x <- rev(x)
    under <- "over"
    over <- "under"
  }
  res <- character(lx - 1)
  res[1] <- paste(under, x[2])
  for (i in 2:(lx - 2)) res[i] <- paste(x[i], between, x[i + 
                                                           1])
  res[lx - 1] <- paste(over, x[lx - 1])
  res
}