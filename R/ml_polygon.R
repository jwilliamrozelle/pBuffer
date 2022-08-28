#' @name ml_polygon
#' @rdname ml_polygon
#' @title  ml_polygon
#'
#' @description  Determines the maximum likely polygon for a given point
#' 
#' @returns returns a dataframe with one row containing the data from most likely health facility joined, and its estimated probability.
#' 
#' @details Replaces the weightedClose function
#'
#'
#' @param bufferPoints.sf Buffer points of a given displaced location containing weights
#' @param polygons.sf Polygons to join the displaced points to. For now, the function assumes there is a column with unique ID for each polygon called "SPAID"
#' @param polygons.id Name of column of polygons.id sf object that contains a unique ID, default is "SPAID".
#' @param weightsCol Name of the column in the bufferPoints.sf object that contains the weights to use, defaults to "layer".
#'
#' @author J.W. Rozelle
#'
#'
#' @export ml_polygon
#' @examples
#'
#' # coming soon!
#' 
#'



ml_polygon <- function(bufferPoints.sf, polygons.sf, polygons.id = "SPAID", weightsCol = "layer") {
  
  intersection <- st_intersection(polygons.sf, bufferPoints.sf) |> suppressWarnings()
  
  # Rename the columns to the columns required by the function
  polygons.sf$SPAID <- polygons.sf[[polygons.id]]
  bufferPoints.sf$layer <- bufferPoints.sf[[weightsCol]]
  
  # reweight to sum of 1
  bufferPoints.sf$layer <- bufferPoints.sf$layer / sum(bufferPoints.sf$layer, na.rm = T)
  
  # non_missing <- length(intersection$SPAID)
  
  # group all intersections by the number of possible polygons
  int_result <- intersection %>% st_drop_geometry() %>% group_by(SPAID) %>% summarize(weight = sum(layer, na.rm = TRUE))
  
  # merge the grouped possible polygons and probabilities with the data in the polygons. Keep only the possible polygons
  includedScores <- merge(int_result, 
                          st_drop_geometry(polygons.sf[c(polygons.id)]), 
                          all.x = TRUE, all.y = FALSE, 
                          by = polygons.id)
  
  # Filter to the maximum likely polygon
  mostLikelyHF.df <- filter(int_result, weight == max(weight, na.rm = T)) |> as.data.frame()
  
  # throw a warning if linked to more than one health facility
  if (nrow(mostLikelyHF) > 1) {
    warning(paste0("Warning: this was linked to more than one feature. The first feature of ", nrow(mostLikelyHF)," is used"))
    # select only the first one linked
    mostLikelyHF.df <- mostLikelyHF[1,]
  }
  
  
  return(mostLikelyHF.df)
}



