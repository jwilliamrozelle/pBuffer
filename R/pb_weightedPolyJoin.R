#' @name pb_weightedPolyJoin
#' @rdname pb_weightedPolyJoin
#' @title  pb_weightedPolyJoin
#'
#' @description For each given displaced community, estimates the raster value based on the probability weighted mean of possible true locations
#' 
#' @returns Returns a dataframe of length nrow(displaced.sf), with a DHSID column containing the unique ID of displaced.sf and SPAID containing the unique ID of the most likely polygon that a given point falls in.
#'
#'
#' @param displaced.sf Displaced coordinates in an sf object type.
#' @param polygons.sf Polygons to join the displaced points to. For now, the function assumes there is a column with unique ID for each polygon called "SPAID"
#' @param metrics Character vector of column names containing the numeric metrics you wish to calculate.
#' @param densityBuffer A density buffer list object created with the pb_integratedDens function.
#' @param adminBound (Optional) the administrative boundary that circumscribes the displacement.
#' @param adminID The unique ID for the adminBound features, defaults to "ID_2".
#' @param n.cores This allows for parallelization using the futures package. This can be slightly unstable, but normally functions well and dramatically speeds compute time when it does. n.cores specifies the number of cores to use. The default is 1, and does not parallelize.
#'
#' @author J.W. Rozelle
#'
#'
#' @export pb_weightedPolyJoin
#' @examples
#'
#' # coming soon!
#' 

pb_weightedPolyJoin <- function(displaced.sf, 
                                polygons.sf, 
                                displaced.id = "DHSID", 
                                metrics = c("sri_score", 
                                            "sri_basicamenities", 
                                            "sri_basicequip", 
                                            "sri_diagcapacity", 
                                            "sri_infprev", 
                                            "sri_med"
                                ), 
                                densityBuffer = densitybuffer, 
                                adminBound = haiti_adm2, 
                                adminID = "ID_2",
                                n.cores = 1) {
  
  displaced.sf$DHSID <- displaced.sf[[displaced.id]]
  
  
  if (n.cores > 1) {
    
    tryCatch({
      
      # set up the workers
      cl <- parallel::makeCluster(n.cores)
      plan(cluster, workers = cl)
      
      closestHFs <- future_apply(st_drop_geometry(displaced.sf), 1, function(x) {
        
        rowDHSID <- x[displaced.id]
        
        # extract an sf object for each row of the data frame
        singleComm <- filter(displaced.sf, DHSID == rowDHSID)
        crs2use <- raster::crs(singleComm)
        
        if (st_coordinates(singleComm)[2] != 0) {
          
          # rasterize the displacement buffer around the single community
          singleDens.raster <- rasterizeDisplacement(
            densitybuffer, 
            initialCoords = st_coordinates(singleComm), 
            inputCRS = crs2use
          )
          
          # Trim the weighted raster
          singleDens.raster <- trimProbBuff(singleDens.raster, adminBound = adminBound)
          
          # turn these into a point object
          singleDens.sf <- st_rasterAsPoint(singleDens.raster)
          rm(singleDens.raster)
          
          # Get the most probable nearby health facility
          probableMetrics.result <- weightedMetrics(singleDens.sf, polygons.sf, metrics = metrics)
          probableMetrics.df <- probableMetrics.result$weightedMetrics.df
          probableMetrics.df$DHSID <- rowDHSID
          probableMetrics.df$n_linked <- nrow(probableMetrics.result$hfWeights)
          
          hfWeights.df <- probableMetrics.result$hfWeights
          hfWeights.df$DHSID <- rowDHSID
          
          
          
        } else {
          probableMetrics.df <- NA
          hfWeights.df <- NA
        }
        
        return(list(probableMetrics.df = probableMetrics.df, hfWeights.df = hfWeights.df))
        
      })
      
    }, 
    finally = {parallel::stopCluster(cl)}
    )
    
    
  } else {
    
    closestHFs <- apply(st_drop_geometry(displaced.sf), 1, function(x) {
      
      rowDHSID <- x[polygons.sf]
      
      # extract an sf object for each row of the data frame
      singleComm <- filter(displaced.sf, DHSID == rowDHSID)
      crs2use <- raster::crs(singleComm)
      
      if (st_coordinates(singleComm)[2] != 0) {
        
        # rasterize the displacement buffer around the single community
        singleDens.raster <- rasterizeDisplacement(
          densitybuffer, 
          initialCoords = st_coordinates(singleComm), 
          inputCRS = crs2use
        )
        
        # Trim the weighted raster
        singleDens.raster <- trimProbBuff(singleDens.raster, adminBound = adminBound)
        
        # turn these into a point object
        singleDens.sf <- st_rasterAsPoint(singleDens.raster)
        rm(singleDens.raster)
        
        # Get the most probable nearby health facility
        probableMetrics.result <- weightedMetrics(singleDens.sf, polygons.sf, metrics = metrics)
        probableMetrics.df <- probableMetrics.result$weightedMetrics.df
        probableMetrics.df$DHSID <- rowDHSID
        probableMetrics.df$n_linked <- nrow(probableMetrics.result$hfWeights)
        
        hfWeights.df <- probableMetrics.result$hfWeights
        hfWeights.df$DHSID <- rowDHSID
        
        
        
      } else {
        probableMetrics.df <- NA
        hfWeights.df <- NA
      }
      
      return(list(probableMetrics.df = probableMetrics.df, hfWeights.df = hfWeights.df))
      
    })
  }
  
  
  # tidy things up, put into a single data frame
  #   Extract the objects of interest from the list
  #     probableMetrics
  probableMetrics.list <- lapply(closestHFs, function(x) {
    x$probableMetrics.df
  }) 
  #     hfWeights
  hfWeights.list <- lapply(closestHFs, function(x) {
    x$hfWeights.df
  }) 
  
  rm(closestHFs)
  
  #   bind them into a dataframe
  probableMetrics.df <- iotools::fdrbind(probableMetrics.list)
  hfWeights.df <- iotools::fdrbind(hfWeights.list)
  
  rm(probableMetrics.list)
  rm(hfWeights.list)
  
  
  
  return(list(probableMetrics.df = probableMetrics.df, hfWeights.df = hfWeights.df))
}
