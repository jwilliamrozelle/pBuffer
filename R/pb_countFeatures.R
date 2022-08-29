#' @name pb_countFeatures
#' @rdname pb_countFeatures
#' @title  pb_countFeatures
#'
#' @description For each given displaced community, estimates the raster value based on the probability weighted mean of possible true locations
#' 
#' @returns a vector of the most likely count of point features within a defined radius of the displaced coordinates.
#'
#'
#' @param displaced.sf Displaced coordinates in an sf object type.
#' @param features2count.sf sf point feature that is counted.
#' @param radiusLength a numeric vector to determine the distance from a displaced community within which the number of features2count.sf will be counted.
#' @param densityBuffer A density buffer list object created with the pb_integratedDens function.
#' @param adminBound (Optional) the administrative boundary that circumscribes the displacement.
#' @param adminID The unique ID for the adminBound features, defaults to "ID_2".
#' @param n.cores This allows for parallelization using the futures package. This can be slightly unstable, but normally functions well and dramatically speeds compute time when it does. n.cores specifies the number of cores to use. The default is 1, and does not parallelize.
#'
#' @author J.W. Rozelle
#'
#'
#' @export pb_countFeatures
#' @examples
#'
#' # coming soon!
#' 
#' 
#' 
#' 


# polygons.sf
# features2count.sf <- htSPAge_t2.sf
# n.cores <- 18
# densitybuffer <- pb_integratedDensity(boundaries = 5e3)


pb_countFeatures <- function(displaced.sf, 
                              features2count.sf, 
                              radiusLength = 1e4,
                              displaced.id = "DHSID", 
                              densityBuffer = densitybuffer, 
                              adminBound = NULL, 
                              adminID = NULL,
                              n.cores = 1) {
  
  
  # ERRORS!
  #   displaced.sf is the wrong class
  if (!"sf" %in% class(displaced.sf)) {
    stop("displaced.sf is not of class sf. It must be point features of class sf")
  }
  #   features2count.sf is the wrong class
  if (!"sf" %in% class(features2count.sf)) {
    stop("features2count.sf is not of class sf. It must be point features of class sf")
  }
  #   adminbBound is the wrong class
  if (!is.null(adminBound) & !"sf" %in% class(adminBound)) {
    stop("adminBound is not of class sf. It must be polygon features of class sf")
  }
  
  #   ensure that the adminBound sf object contains the unique ID
  if (!adminID %in% names(adminBound)) {
    stop(paste0(adminID, " is not found in the adminBound sf object."))
  }
  
  if (!is.null(adminBound)) {
    if (sum(duplicated(adminBound[[adminID]])) > 0) {
      stop("adminID must specify a uniquely valid ID for the adminbound object")
    }
  }
  
  #   ensure that the displaced.sf sf object contains the unique ID
  if (!displaced.id %in% names(displaced.sf)) {
    stop(paste0(displaced.id, " is not found in the displaced.sf sf object."))
  }
  
  if (sum(duplicated(displaced.sf[[displaced.id]])) > 0) {
    stop("displaced.id must specify a uniquely valid ID for the displaced.sf object")
  }
  
  
  # Assign unique id name to "DHSID" column
  displaced.sf$DHSID <- displaced.sf[[displaced.id]]
  # assign the adminID character value ot adminID variable name in dataframe
  adminBound$adminID <- st_drop_geometry(adminBound)[[adminID]]
  
  
  if (n.cores > 1) {
    
    tryCatch({
      
      # set up the workers
      cl <- parallel::makeCluster(n.cores)
      plan(cluster, workers = cl)
      
      numFeatures <- future_apply(st_drop_geometry(displaced.sf), 1, function(x) {
        
        rowDHSID <- x[[displaced.id]]
        
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
          
          if (!is.null(adminBound)) {
            # get the single admin boundary for a community
            singleAdminBound <- suppressWarnings(sf::st_intersection(singleComm, adminBound))
            singleAdminBound.poly <- dplyr::filter(adminBound, adminID == singleAdminBound$adminID[1])
            
            rm(singleAdminBound)
            
            # Trim the weighted raster
            singleDens.raster <- trimProbBuff(singleDens.raster, adminBound = singleAdminBound.poly)
            
            rm(singleAdminBound.poly)
            
          }
          
          # turn these into a point object
          singleDens.sf <- st_rasterAsPoint(singleDens.raster)
          rm(singleDens.raster)
          
          singleBuffs.sf <- st_buffer(singleDens.sf, radiusLength)
          
          # intersection.sf <- st_intersection(features2count.sf, singleBuffs.sf)
          
          singleBuffs.sf$pt_count <- lengths(st_intersects(singleBuffs.sf, features2count.sf))
          
          numFeatures.df <- st_drop_geometry(singleBuffs.sf) %>% group_by(pt_count) %>% summarise(weights = sum(layer, na.rm = T)) |> as.data.frame()
          
          
          # Get the maximum likelihood number of features as single row data.frame
          ml_featureCount.df <- filter(numFeatures.df, weights == max(weights, na.rm = T)) |> as.data.frame()
          # add the unique ID to the single row data.frame
          ml_featureCount.df[[displaced.id]] <- rowDHSID
          
          # 
          numFeatures.df[[displaced.id]] <- rowDHSID
          
          
          
          
          
        } else {
          ml_featureCount.df <- NA
          possibleNumFeatures.df <- NA
        }
        
        return(list(ml_featureCount.df = ml_featureCount.df, possibleNumFeatures.df = numFeatures.df))
        
      })
      
    }, 
    finally = {parallel::stopCluster(cl)}
    )
    
    
  } else {
    
    numFeatures <- apply(st_drop_geometry(displaced.sf), 1, function(x) {
      
      rowDHSID <- x[[displaced.id]]
      
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
        
        if (!is.null(adminBound)) {
          # get the single admin boundary for a community
          singleAdminBound <- suppressWarnings(sf::st_intersection(singleComm, adminBound))
          singleAdminBound.poly <- dplyr::filter(adminBound, adminID == singleAdminBound$adminID[1])
          
          rm(singleAdminBound)
          
          # Trim the weighted raster
          singleDens.raster <- trimProbBuff(singleDens.raster, adminBound = singleAdminBound.poly)
          
          rm(singleAdminBound.poly)
          
        }
        
        # turn these into a point object
        singleDens.sf <- st_rasterAsPoint(singleDens.raster)
        rm(singleDens.raster)
        
        singleBuffs.sf <- st_buffer(singleDens.sf, radiusLength)
        
        # intersection.sf <- st_intersection(features2count.sf, singleBuffs.sf)
        
        singleBuffs.sf$pt_count <- lengths(st_intersects(singleBuffs.sf, features2count.sf))
        
        numFeatures.df <- st_drop_geometry(singleBuffs.sf) %>% group_by(pt_count) %>% summarise(weights = sum(layer, na.rm = T)) |> as.data.frame()
        
        
        # Get the maximum likelihood number of features as single row data.frame
        ml_featureCount.df <- filter(numFeatures.df, weights == max(weights, na.rm = T)) |> as.data.frame()
        # add the unique ID to the single row data.frame
        ml_featureCount.df[[displaced.id]] <- rowDHSID
        
        # 
        numFeatures.df[[displaced.id]] <- rowDHSID
        
        
        
        
        
      } else {
        ml_featureCount.df <- NA
        possibleNumFeatures.df <- NA
      }
      
      return(list(ml_featureCount.df = ml_featureCount.df, possibleNumFeatures.df = numFeatures.df))
      
    })
  }
  
  
  # tidy things up, put into a single data frame
  #   Extract the objects of interest from the list
  #     probableMetrics
  ml_featureCount.list <- lapply(numFeatures, function(x) {
    x$ml_featureCount.df
  }) 
  #     hfWeights
  possibleNumFeatures.list <- lapply(numFeatures, function(x) {
    x$possibleNumFeatures.df
  }) 
  
  rm(numFeatures)
  
  #   bind them into a dataframe
  ml_featureCount.df <- iotools::fdrbind(ml_featureCount.list)
  possibleNumFeatures.df <- iotools::fdrbind(possibleNumFeatures.list)
  
  rm(ml_featureCount.list, possibleNumFeatures.list)
  
  return(list(ml_featureCount.df = ml_featureCount.df, possibleNumFeatures.df = possibleNumFeatures.df))
}