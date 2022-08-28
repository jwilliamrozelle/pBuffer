#' @name trimProbBuff
#' @rdname trimProbBuff
#' @title  trimProbBuff
#'
#' @description Trims and reweights probability buffer (from a raster object) by a polygon that it is limited to
#' 
#' @returns Returns a trimmed object of class raster, with values summed to 1.
#' 
#' @param probRaster Probability raster, object of class raster
#' @param adminBound The administrative / polygon boundary of class SF that the raster should be masked by
#'
#' @author J.W. Rozelle
#'
#'
#' @export trimProbBuff
#' @examples
#'
#' # coming soon!
#' 
#'




trimProbBuff <- function(probRaster = NULL, adminBound = NULL) {
  
  # mask the raster
  trimmedProbBuff <- raster::mask(probRaster, mask = adminBound)
  
  # reweight
  trimmedProbBuff[] <- trimmedProbBuff[] / sum(trimmedProbBuff[], na.rm = T)
  
  return(trimmedProbBuff)
  
}

