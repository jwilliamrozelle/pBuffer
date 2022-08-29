#' @name pb_quickRasterValJoin
#' @rdname pb_quickRasterValJoin
#' @title  pb_quickRasterValJoin
#'
#' @description For each given displaced community, estimates the raster value based on the probability weighted mean of possible true locations
#' 
#' @returns Returns a dataframe of length nrow(displaced.sf), with a DHSID column containing the unique ID of displaced.sf and SPAID containing the unique ID of the most likely polygon that a given point falls in.
#'
#'
#' @param displaced.sf Displaced coordinates in an sf object type.
#' @param polygons.sf Polygons to join the displaced points to. For now, the function assumes there is a column with unique ID for each polygon called "SPAID"
#' @param displaced.id The name of the column that contains a unique ID for the displaced communities.
#' @param densityBuffer A density buffer list object created with the pb_integratedDens or pb_Density functions.
#' @param adminBound (Optional) the administrative boundary that circumscribes the displacement
#' @param adminID The unique ID for the adminBound featurs, defaults to "ID_2".
#' @param n.cores This allows for parallelization using the futures package. This can be slightly unstable, but normally functions well and dramatically speeds compute time when it does. n.cores specifies the number of cores to use. The default is 1, and does not parallelize.
#'
#' @author J.W. Rozelle
#'
#'
#' @export pb_quickRasterValJoin
#' @examples
#'
#' # coming soon!
#' 


pb_quickRasterValJoin <- function(displaced.sf, 
                                   inputRaster, 
                                   bufferlengths = 5000, 
                                   adminBound = NULL, 
                                   adminID = NULL, 
                                   n.cores = 1) {
  
  # ERRORS!
  #   displaced.sf is the wrong class
  if (!"sf" %in% class(displaced.sf)) {
    stop("displaced.sf is not of class sf. It must be point features of class sf")
  }
  #   features2count.sf is the wrong class
  if (!"RasterLayer" %in% class(inputRaster)) {
    stop("inputRaster is not of class 'raster'. It must be point features of class 'raster'")
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
  
  warning("If the input raster has missing values, any probability buffer points that fall on the missing raster values will be excluded from the weighted score.")
  
  # get the resolution of the raster
  raster.res <- raster::res(inputRaster)
  crs2use <- raster::crs(singleComm)
  # pixelArea <- raster.res[1]*raster.res[2]
  
  # assign the adminID character value ot adminID variable name in dataframe
  adminBound$adminID <- st_drop_geometry(adminBound)[[adminID]]
  
  pbBuffer.list <- list()
  
  for (i in 1:length(bufferlengths)) {
    
    # probability density function
    buffer_pb_function <- function(x, y) {
      
      result <- ifelse(sqrt(x^2 + y^2) <= bufferlengths[i], 1/(bufferlengths[i] * 2 * pi * sqrt(x^2+y^2)), 0)
      
      return(result)
      
    }
    
    quadrantRowN <- (ceiling(bufferlengths[i] / raster.res[2]))
    quadrantColN <- (ceiling(bufferlengths[i] / raster.res[1]))
    
    # create a raster of the appropriate size given the bufferlength and raster resolution
    # pbBuffer.matrix <- matrix(data = 0, nrow = quadrantRowN*2, ncol = quadrantColN*2)
    
    
    
    q4_matrix <- matrix(data = 0, nrow = quadrantRowN + 2, ncol = quadrantColN + 2)
    
    xadd <- raster.res[1] / 2
    yadd <- raster.res[2] / 2
    
    
    # do the weight calculation
    for (q4_row in 1:(quadrantRowN+2)) {
      for (q4_col in 1:(quadrantColN+2)) {
        
        # Set limits for integrals
        if (q4_row - 1 > 0 | q4_col - 1 > 0) {
          ylimits <- c(((q4_row-1)*raster.res[2] - yadd) , ((q4_row-1)*raster.res[2] + yadd))
          xlimits <- c(((q4_col-1)*raster.res[1] - xadd) , ((q4_col-1)*raster.res[1] + xadd))
          
          integralValue <- calculus::integral(buffer_pb_function, bounds = list(x = xlimits, y = ylimits), coordinates = "cartesian")
          
          q4_matrix[q4_row, q4_col] <- integralValue$value
          
        } else {
          ylimits <- c(((q4_row-1)*raster.res[2] - yadd) , ((q4_row-1)*raster.res[2] + yadd))
          xlimits <- c(((q4_col-1)*raster.res[1] - xadd-0.01) , ((q4_col-1)*raster.res[1] + xadd))
          
          integralValue <- calculus::integral(buffer_pb_function, bounds = list(x = xlimits, y = ylimits), coordinates = "cartesian")
          
          # 
          q4_matrix[1, 1] <- integralValue$value
        }
        
      }
    }
    
    
    q3_matrix <- q4_matrix[, ncol(q4_matrix):2]
    q2_matrix <- q4_matrix[nrow(q4_matrix):2, ncol(q4_matrix):2]
    q1_matrix <- q4_matrix[nrow(q4_matrix):2, ]
    
    # bind the east and west halves of the north and south hemicircles
    north <- cbind(q2_matrix, q1_matrix)
    south <- cbind(q3_matrix, q4_matrix)
    
    # bind the north and south
    densMatrix <- rbind(north, south)
    
    rm(q1_matrix, q2_matrix, q3_matrix, q4_matrix, north, south)
    
    
    pbBuffer.list[[i]] <- densMatrix
  }
  
  
  
  if (n.cores > 1) {
    
    tryCatch({
      
      # set up the workers
      cl <- parallel::makeCluster(n.cores)
      plan(cluster, workers = cl)
      
      weightedRasterVals <- future_lapply(st_drop_geometry(displaced.sf)$DHSID, function(x) {
        
        rowDHSID <- x
        
        # extract an sf object for each row of the data frame
        singleComm <- filter(displaced.sf, DHSID == rowDHSID)
        
        # throw an error if singleComm has multiple observations
        if (nrow(singleComm) > 1) {
          stop(paste0("There are multiple observations with the same DHSID: ", rowDHSID))
        }
        
        
        if (st_coordinates(singleComm)[2] != 0) {
          
          singleBuffer.list <- list()
          singleBufferRaster.list <- list()
          for (i in 1:length(bufferlengths)){
            
            # create a buffer for the community
            singleBuffer.list[[i]] <- st_buffer(singleComm, bufferlengths[i]+(max(raster.res)*2)) |> suppressWarnings()
            
            # If the administrative boundary is defined
            if (!is.null(adminBound)) {
              # get the single admin boundary for a community
              singleAdminBound <- suppressWarnings(sf::st_intersection(singleComm, adminBound))
              singleAdminBound.poly <- dplyr::filter(adminBound, adminID == singleAdminBound$adminID[1])
              
              # change all values that don't fall within the buffer to 0
              singleBufferRaster.list[[i]] <- raster::mask(inputRaster, singleBuffer.list[[i]]) |> suppressWarnings()
              
              # Trim the weighted raster
              singleBufferRaster.list[[i]] <- trimProbBuff(singleBufferRaster.list[[i]], adminBound = singleAdminBound.poly)
              rm(singleAdminBound, singleAdminBound.poly)
            }
            
            
            # 
            # # mask based on the administrative bounadries
            # singleBufferRaster.list[[i]] <- raster::mask(singleBufferRaster.list[[i]], singleAdminBound.poly)
            
            # trim the excess NA's
            singleBufferRaster.list[[i]] <- raster::trim(singleBufferRaster.list[[i]]) |> suppressWarnings()
            
            # rasterRows 
            rasterRows <- dim(singleBufferRaster.list[[i]])[1]
            rasterCols <- dim(singleBufferRaster.list[[i]])[2]
            
            # RASTER CENTER POINT
            singleCoords <- st_coordinates(singleComm)
            #   cell of center point
            centerCell <- raster::cellFromXY(singleBufferRaster.list[[i]], singleCoords)
            centerRowCol <- raster::rowColFromCell(singleBufferRaster.list[[i]], centerCell)
            
            centerCell_xy <- raster::xyFromCell(singleBufferRaster.list[[i]], centerCell)
            
            pb_singleComm.raster <- raster(pbBuffer.list[[i]])
            
            pb_xmin <- centerCell_xy[1,"x"] - raster.res[1]*(ncol(pbBuffer.list[[i]])/2)
            pb_xmax <- centerCell_xy[1,"x"] + raster.res[1]*(ncol(pbBuffer.list[[i]])/2)
            pb_ymin <- centerCell_xy[1,"y"] - raster.res[2]*(nrow(pbBuffer.list[[i]])/2)
            pb_ymax <- centerCell_xy[1,"y"] + raster.res[2]*(nrow(pbBuffer.list[[i]])/2)
            
            extent(pb_singleComm.raster) <- c(pb_xmin, pb_xmax, pb_ymin, pb_ymax)
            crs(pb_singleComm.raster) <- crs(crs2use)
            raster::values(pb_singleComm.raster)[raster::values(pb_singleComm.raster) == 0] <- NA
            
            # get the rasters to line up
            singleBufferRaster.list[[i]] <- raster::extend(singleBufferRaster.list[[i]], pb_singleComm.raster)
            pb_singleComm.raster <- raster::extend(pb_singleComm.raster, singleBufferRaster.list[[i]])
            
            
            
            rm(pb_xmin, pb_xmax, pb_ymin, pb_ymax)
            
          }
          
          singleCommVals.df <- data.frame(pixelProb = raster::values(pb_singleComm.raster), 
                                          pixelVal = raster::values(singleBufferRaster.list[[i]])
          )
          singleCommVals.df$weights <- singleCommVals.df$pixelProb # / sum(pixelProb, na.rm = T)
          singleCommVals.df$DHSID <- rowDHSID
          
          singleCommVals.df <- filter(singleCommVals.df, !is.na(pixelProb))
          singleCommVals.df <- filter(singleCommVals.df, !is.na(pixelVal))
          
          estimate <- weighted.mean(singleCommVals.df$pixelVal, singleCommVals.df$weights, na.rm = T)
          estimate.df <- data.frame(estimate)
          estimate.df$DHSID <- rowDHSID
          
          
          return(list(estimate.df = estimate.df, commVals.df = singleCommVals.df))
        } else {
          return(NA)
        }
      })
      
    }, 
    
    finally = {parallel::stopCluster(cl)}
    )
    
    
    
    
    
    
    # if not parallelized...  
  } else {
    
    weightedRasterVals <- lapply(st_drop_geometry(displaced.sf)["DHSID"], function(x) {
      
      rowDHSID <- x
      
      # extract an sf object for each row of the data frame
      singleComm <- filter(displaced.sf, DHSID == rowDHSID)
      
      # throw an error if singleComm has multiple observations
      if (nrow(singleComm) > 1) {
        stop(paste0("There are multiple observations with the same DHSID: ", rowDHSID))
      }
      
      
      if (st_coordinates(singleComm)[2] != 0) {
        
        singleBuffer.list <- list()
        singleBufferRaster.list <- list()
        for (i in 1:length(bufferlengths)){
          
          # create a buffer for the community
          singleBuffer.list[[i]] <- st_buffer(singleComm, bufferlengths[i]+(max(raster.res)*2)) |> suppressWarnings()
          
          # If the administrative boundary is defined
          if (!is.null(adminBound)) {
            # get the single admin boundary for a community
            singleAdminBound <- suppressWarnings(sf::st_intersection(singleComm, adminBound))
            singleAdminBound.poly <- dplyr::filter(adminBound, adminID == singleAdminBound$adminID[1])
            
            # change all values that don't fall within the buffer to 0
            singleBufferRaster.list[[i]] <- raster::mask(inputRaster, singleBuffer.list[[i]]) |> suppressWarnings()
            
            # Trim the weighted raster
            singleBufferRaster.list[[i]] <- trimProbBuff(singleBufferRaster.list[[i]], adminBound = singleAdminBound.poly)
            rm(singleAdminBound, singleAdminBound.poly)
          }
          
          
          # 
          # # mask based on the administrative bounadries
          # singleBufferRaster.list[[i]] <- raster::mask(singleBufferRaster.list[[i]], singleAdminBound.poly)
          
          # trim the excess NA's
          singleBufferRaster.list[[i]] <- raster::trim(singleBufferRaster.list[[i]]) |> suppressWarnings()
          
          # rasterRows 
          rasterRows <- dim(singleBufferRaster.list[[i]])[1]
          rasterCols <- dim(singleBufferRaster.list[[i]])[2]
          
          # RASTER CENTER POINT
          singleCoords <- st_coordinates(singleComm)
          #   cell of center point
          centerCell <- raster::cellFromXY(singleBufferRaster.list[[i]], singleCoords)
          centerRowCol <- raster::rowColFromCell(singleBufferRaster.list[[i]], centerCell)
          
          centerCell_xy <- raster::xyFromCell(singleBufferRaster.list[[i]], centerCell)
          
          pb_singleComm.raster <- raster(pbBuffer.list[[i]])
          
          pb_xmin <- centerCell_xy[1,"x"] - raster.res[1]*(ncol(pbBuffer.list[[i]])/2)
          pb_xmax <- centerCell_xy[1,"x"] + raster.res[1]*(ncol(pbBuffer.list[[i]])/2)
          pb_ymin <- centerCell_xy[1,"y"] - raster.res[2]*(nrow(pbBuffer.list[[i]])/2)
          pb_ymax <- centerCell_xy[1,"y"] + raster.res[2]*(nrow(pbBuffer.list[[i]])/2)
          
          extent(pb_singleComm.raster) <- c(pb_xmin, pb_xmax, pb_ymin, pb_ymax)
          crs(pb_singleComm.raster) <- crs(crs2use)
          raster::values(pb_singleComm.raster)[raster::values(pb_singleComm.raster) == 0] <- NA
          
          # get the rasters to line up
          singleBufferRaster.list[[i]] <- raster::extend(singleBufferRaster.list[[i]], pb_singleComm.raster)
          pb_singleComm.raster <- raster::extend(pb_singleComm.raster, singleBufferRaster.list[[i]])
          
          
          
          rm(pb_xmin, pb_xmax, pb_ymin, pb_ymax)
          
        }
        
        singleCommVals.df <- data.frame(pixelProb = raster::values(pb_singleComm.raster), 
                                        pixelVal = raster::values(singleBufferRaster.list[[i]])
        )
        singleCommVals.df$weights <- singleCommVals.df$pixelProb # / sum(pixelProb, na.rm = T)
        singleCommVals.df$DHSID <- rowDHSID
        
        singleCommVals.df <- filter(singleCommVals.df, !is.na(pixelProb))
        singleCommVals.df <- filter(singleCommVals.df, !is.na(pixelVal))
        
        estimate <- weighted.mean(singleCommVals.df$pixelVal, singleCommVals.df$weights, na.rm = T)
        estimate.df <- data.frame(estimate)
        estimate.df$DHSID <- rowDHSID
        
        
        return(list(estimate.df = estimate.df, commVals.df = singleCommVals.df))
      } else {
        return(NA)
      }
      
    })
  }
  
  # tidy things up, put into a single data frame
  #   Extract the objects of interest from the list
  #     probableMetrics
  estimate.list <- lapply(weightedRasterVals, function(x) {
    x$estimate.df
  }) 
  #     hfWeights
  commVals.list <- lapply(weightedRasterVals, function(x) {
    x$commVals.df
  }) 
  
  rm(weightedRasterVals)
  
  #   bind them into a dataframe
  estimate.df <- iotools::fdrbind(estimate.list)
  commVals.df <- iotools::fdrbind(commVals.list)
  
  rm(estimate.list)
  rm(commVals.list)
  
  return(list(estimate.df = estimate.df, commVals.df = commVals.df))
}
