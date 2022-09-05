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
#' @param limitDist a numeric vector to determine the distance from a displaced community within which the number of features2count.sf will be counted.
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


pb_radiusLimitJoin <- function(displaced.sf, 
                          features2count.sf, 
                          limitDist = 5e3,
                          displaced.id = "DHSID", 
                          features2count.id = "SPAID",
                          metrics = c("sri_score", 
                                      "sri_basicamenities", 
                                      "sri_basicequip", 
                                      "sri_diagcapacity", 
                                      "sri_infprev", 
                                      "sri_med"),
                          densityBuffer = densitybuffer, 
                          adminBound = NULL, 
                          adminID = NULL,
                          weightedMedian = TRUE,
                          probCuts = c(0.1, 0.9),
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
  
  # ensure that metrics can be found in the 
  if (!F %in% (metrics %in% names(features2count.sf))) {
    stop("One or more column names specified in metrics cannot be found in the features2count.sf")
  }
  
  # Assign unique id name to "DHSID" column
  displaced.sf$DHSID <- displaced.sf[[displaced.id]]
  # assign unique ID name to the "SPAID_use" column
  features2count.sf$SPAID_use <- features2count.sf[[features2count.id]]
  
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
          singleDens.sf$dispid <- row.names(singleDens.sf)
          
          rm(singleDens.raster)
          
          # calculate the distance matrix
          distance.mtx <- st_distance(singleDens.sf, features2count.sf) |> as.matrix()
          
          # name the rows and columns
          row.names(distance.mtx) <- singleDens.sf$dispid
          colnames(distance.mtx) <- features2count.sf$SPAID_use
          
          # get ids back into the distance dataframe
          distance.df <- as.data.frame(distance.mtx)
          distance.df$DHSID <- rowDHSID
          distance.df$dispid <- row.names(distance.mtx)
          
          # clean up after meeself
          rm(distance.mtx)
          
          # merge the weights and the new values
          singleDens.sf <- merge(singleDens.sf, distance.df, by = "dispid")
          
          singleDens.df <- sf::st_drop_geometry(singleDens.sf)
          rm(singleDens.sf)
          
          featureDistances.list <- lapply(features2count.sf[[features2count.id]], function(x) {
            
            # extract only the estimated distance for a pb point and the weight
            singleFeatureDistance.df <- singleDens.df[c(x, "layer")]
            # update name to distance for code readability
            singleFeatureDistance.df$distance <- as.numeric(singleFeatureDistance.df[[x]])
            
            medianDist <- spatstat.geom::weighted.median(singleFeatureDistance.df$distance, singleFeatureDistance.df$layer)
            # credibleInterval
            distCrI <- spatstat.geom::weighted.quantile(singleFeatureDistance.df$distance, 
                                                        singleFeatureDistance.df$layer, 
                                                        # calculate middle credible interval
                                                        probs = probCuts
            )
            # name credibleIntervalvalues
            names(distCrI) <- sapply(names(distCrI), function(x) {paste0("limit_", x)})
            distCrI.df <- as.data.frame(t(distCrI))
            rm(distCrI)
            
            
            probInLimit <- sum(ifelse(singleFeatureDistance.df$distance <= limitDist & 
                                        !is.na(singleFeatureDistance.df$distance) & 
                                        !is.na(singleFeatureDistance.df$layer), 
                                      singleFeatureDistance.df$layer, 
                                      0)
            )
            
            hf_weightedDist.df <- data.frame(likelyDist = medianDist, probInLimit = probInLimit)
            
            hf_weightedDist.df[[features2count.id]] <- x
            
            hf_weightedDist.df <- cbind(hf_weightedDist.df, distCrI.df)
            
            return(hf_weightedDist.df)
          })
          
          # bind this into a single dataframe
          hf_weightedDist.df <- iotools::fdrbind(featureDistances.list)
          # filter to only those that have a possibility of being within the buffer
          hf_weightedDist.df <- filter(hf_weightedDist.df, probInLimit > 0)
          # add in metrics of interest
          
          if (!is.null(metrics)) {
            hf_weightedDist.df <- merge(hf_weightedDist.df, 
                                        sf::st_drop_geometry(features2count.sf)[c(features2count.id, metrics)], 
                                        by = features2count.id, 
                                        all.x = T, 
                                        all.y = F
            )
          }
          
        } else {
          hf_weightedDist.df <- NA
        }
        
        return(hf_weightedDist.df)
        
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
        singleDens.sf$dispid <- row.names(singleDens.sf)
        
        rm(singleDens.raster)
        
        # calculate the distance matrix
        distance.mtx <- st_distance(singleDens.sf, features2count.sf) |> as.matrix()
        
        # name the rows and columns
        row.names(distance.mtx) <- singleDens.sf$dispid
        colnames(distance.mtx) <- features2count.sf$SPAID_use
        
        # get ids back into the distance dataframe
        distance.df <- as.data.frame(distance.mtx)
        distance.df$DHSID <- rowDHSID
        distance.df$dispid <- row.names(distance.mtx)
        
        # clean up after meeself
        rm(distance.mtx)
        
        # merge the weights and the new values
        singleDens.sf <- merge(singleDens.sf, distance.df, by = "dispid")
        
        singleDens.df <- sf::st_drop_geometry(singleDens.sf)
        rm(singleDens.sf)
        
        featureDistances.list <- lapply(features2count.sf[[features2count.id]], function(x) {
          
          # extract only the estimated distance for a pb point and the weight
          singleFeatureDistance.df <- singleDens.df[c(x, "layer")]
          # update name to distance for code readability
          singleFeatureDistance.df$distance <- as.numeric(singleFeatureDistance.df[[x]])
          
          medianDist <- spatstat.geom::weighted.median(singleFeatureDistance.df$distance, singleFeatureDistance.df$layer)
          # credibleInterval
          distCrI <- spatstat.geom::weighted.quantile(singleFeatureDistance.df$distance, 
                                                      singleFeatureDistance.df$layer, 
                                                      # calculate middle credible interval
                                                      probs = probCuts
          )
          # name credibleIntervalvalues
          names(distCrI) <- sapply(names(distCrI), function(x) {paste0("limit_", x)})
          distCrI.df <- as.data.frame(t(distCrI))
          rm(distCrI)
          
          
          probInLimit <- sum(ifelse(singleFeatureDistance.df$distance <= limitDist & 
                                      !is.na(singleFeatureDistance.df$distance) & 
                                      !is.na(singleFeatureDistance.df$layer), 
                                    singleFeatureDistance.df$layer, 
                                    0)
          )
          
          hf_weightedDist.df <- data.frame(likelyDist = medianDist, probInLimit = probInLimit)
          
          hf_weightedDist.df[[features2count.id]] <- x
          
          hf_weightedDist.df <- cbind(hf_weightedDist.df, distCrI.df)
          
          return(hf_weightedDist.df)
        })
        
        # bind this into a single dataframe
        hf_weightedDist.df <- iotools::fdrbind(featureDistances.list)
        # filter to only those that have a possibility of being within the buffer
        hf_weightedDist.df <- filter(hf_weightedDist.df, probInLimit > 0)
        # add in metrics of interest
        
        if (!is.null(metrics)) {
          hf_weightedDist.df <- merge(hf_weightedDist.df, 
                                      sf::st_drop_geometry(features2count.sf)[c(features2count.id, metrics)], 
                                      by = features2count.id, 
                                      all.x = T, 
                                      all.y = F
          )
        }
        
      } else {
        hf_weightedDist.df <- NA
      }
      
      return(hf_weightedDist.df)
      
    })
  }
  
  # put into a nice, tidy table
  hf_weightedDist.df <- iotools::fdrbind(numFeatures)
  
  return(hf_weightedDist.df)
}