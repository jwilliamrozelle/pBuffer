#' @name pb_rasterExtract
#' @rdname pb_polyJoin
#' @title  pb_polyJoin
#'
#' @description  Estimates th
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
#' @export pb_Density
#' @examples
#'
#' # coming soon!
#' 


pb_rasterExtract <- function(displaced.sf, inputRaster, bufferlengths = 5000, adminBound = haiti_adm2, adminID = "ID_2", n.cores = 1) {
  
  warning("If the input raster has missing values, any probability buffer points that fall on the missing raster values will be excluded from the weighted score.")
  
  # get the resolution of the raster
  raster.res <- raster::res(inputRaster)
  crs2use <- raster::crs(singleComm)
  # pixelArea <- raster.res[1]*raster.res[2]
  
  
  
  if (n.cores > 1) {
    
    tryCatch({
      
      # set up the workers
      cl <- parallel::makeCluster(n.cores)
      plan(cluster, workers = cl)
      
      
      # assign the adminID character value ot adminID variable name in dataframe
      adminBound$adminID <- st_drop_geometry(adminBound)[[adminID]]
      
      weightedRasterVals <- future_lapply(st_drop_geometry(displaced.sf)$DHSID, function(x) {
        
        rowDHSID <- x
        
        # extract an sf object for each row of the data frame
        singleComm <- filter(displaced.sf, DHSID == rowDHSID)
        
        # throw an error if singleComm has multiple observations
        if (nrow(singleComm) > 1) {
          stop("There are multiple observations with the same DHSID.")
        }
        
        if (st_coordinates(singleComm)[2] != 0) {
          
          singleBuffer.list <- list()
          singleBufferRaster.list <- list()
          for (i in 1:length(bufferlengths)){
            
            # define the pb function for the given boundary
            buffer_pb_function <- function(x, y) {
              
              result <- ifelse(sqrt(x^2 + y^2) <= bufferlengths[i], 1/(bufferlengths[i] * 2 * pi * sqrt(x^2+y^2)), 0)
              
              return(result)
              
            }
            
            
            # create a buffer for the community
            singleBuffer.list[[i]] <- st_buffer(singleComm, bufferlengths[i]+max(raster.res)) |> suppressWarnings()
            
            # change all values that don't fall within the buffer to 0
            singleBufferRaster.list[[i]] <- raster::mask(inputRaster, singleBuffer.list[[i]]) |> suppressWarnings()
            
            # trim the excess NA's
            singleBufferRaster.list[[i]] <- raster::trim(singleBufferRaster.list[[i]]) |> suppressWarnings()
            
            # rasterRows 
            rasterRows <- dim(singleBufferRaster.list[[i]])[1]
            rasterCols <- dim(singleBufferRaster.list[[i]])[1]
            
            # RASTER CENTER POINT
            singleCoords <- st_coordinates(singleComm)
            #   cell of center point
            centerCell <- raster::cellFromXY(singleBufferRaster.list[[i]], singleCoords)
            centerRowCol <- raster::rowColFromCell(singleBufferRaster.list[[i]], centerCell)
            
            locInCenterPix <- singleCoords - raster::xyFromCell(singleBufferRaster.list[[i]], centerCell)
            
            # start empty list of pixel probability
            pixelProb <- c()
            pixelVal <- c()
            weightsMx.list <- list()
            
            # for each cell in the raster
            for(cell in 1:(rasterRows*rasterCols)) {
              # get the x and y of the 
              cellRowCol <- rowColFromCell(singleBufferRaster.list[[i]], cell)
              
              
              
              # get the x and y distnace from the center in map units
              cellRelativePos <- singleCoords - raster::xyFromCell(singleBufferRaster.list[[i]], cell)
              
              # get the inner and outer xy's
              outerXY <- cellRelativePos + raster.res/2
              innerXY <- cellRelativePos - raster.res/2
              
              xlimits <- sort(c(outerXY[,"X"], innerXY[,"X"]))
              ylimits <- sort(c(outerXY[,"Y"], innerXY[,"Y"]))
              
              
              if(!setequal(cellRelativePos, c(0,0))) {
                # get the rectangular integral for for the raster cell
                integralValue <- calculus::integral(buffer_pb_function, bounds = list(x = xlimits, y = ylimits), coordinates = "cartesian")
              } else {
                
                xlimits[1] <- xlimits[1]+0.1
                integralValue <- calculus::integral(buffer_pb_function, bounds = list(x = xlimits, y = ylimits), coordinates = "cartesian")
              }
              
              
              
              pixelProb[cell] <- integralValue$value
              pixelVal[cell] <- singleBufferRaster.list[[i]][cellRowCol]
              
              
              
              
              
            }
            
            
          }
          
          singleCommVals.df <- data.frame(pixelProb, pixelVal)
          singleCommVals.df$weights <- pixelProb # / sum(pixelProb, na.rm = T)
          singleCommVals.df$DHSID <- rowDHSID
          
          estimate <- weighted.mean(singleCommVals.df$pixelVal, singleCommVals.df$weights, na.rm = T)
          estimate.df <- data.frame(estimate)
          estimate.df$DHSID <- rowDHSID
          
          
          return(list(estimate.df = estimate.df, commVals.df = singleCommVals.df))
        } else {
          return(NA)
        }
      })
    },
    
    # success or failure, stop the cluster
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
        stop("There are multiple observations with the same DHSID.")
      }
      
      if (st_coordinates(singleComm)[2] != 0) {
        
        singleBuffer.list <- list()
        singleBufferRaster.list <- list()
        for (i in 1:length(bufferlengths)){
          
          # define the pb function for the given boundary
          buffer_pb_function <- function(x, y) {
            
            result <- ifelse(sqrt(x^2 + y^2) <= bufferlengths[i], 1/(bufferlengths[i] * 2 * pi * sqrt(x^2+y^2)), 0)
            
            return(result)
            
          }
          
          
          # create a buffer for the community
          singleBuffer.list[[i]] <- st_buffer(singleComm, bufferlengths[i]) |> suppressWarnings()
          
          # change all values that don't fall within the buffer to 0
          singleBufferRaster.list[[i]] <- raster::mask(inputRaster, singleBuffer.list[[i]]) |> suppressWarnings()
          
          # trim the excess NA's
          singleBufferRaster.list[[i]] <- raster::trim(singleBufferRaster.list[[i]]) |> suppressWarnings()
          
          # rasterRows 
          rasterRows <- dim(singleBufferRaster.list[[i]])[1]
          rasterCols <- dim(singleBufferRaster.list[[i]])[1]
          
          # RASTER CENTER POINT
          singleCoords <- st_coordinates(singleComm)
          #   cell of center point
          centerCell <- raster::cellFromXY(singleBufferRaster.list[[i]], singleCoords)
          centerRowCol <- raster::rowColFromCell(singleBufferRaster.list[[i]], centerCell)
          
          locInCenterPix <- singleCoords - raster::xyFromCell(singleBufferRaster.list[[i]], centerCell)
          
          # start empty list of pixel probability
          pixelProb <- c()
          pixelVal <- c()
          weightsMx.list <- list()
          
          # for each cell in the raster
          for(cell in 1:(rasterRows*rasterCols)) {
            # get the x and y of the 
            cellRowCol <- rowColFromCell(singleBufferRaster.list[[i]], cell)
            
            
            
            # get the x and y distnace from the center in map units
            cellRelativePos <- singleCoords - raster::xyFromCell(singleBufferRaster.list[[i]], cell)
            
            # get the inner and outer xy's
            outerXY <- cellRelativePos + raster.res/2
            innerXY <- cellRelativePos - raster.res/2
            
            xlimits <- sort(c(outerXY[,"X"], innerXY[,"X"]))
            ylimits <- sort(c(outerXY[,"Y"], innerXY[,"Y"]))
            
            
            
            if(!setequal(cellRelativePos, c(0,0))) {
              # get the rectangular integral for for the raster cell
              integralValue <- calculus::integral(buffer_pb_function, bounds = list(x = xlimits, y = ylimits), coordinates = "cartesian")
            } else {
              
              xlimits[1] <- xlimits[1]+0.1
              integralValue <- calculus::integral(buffer_pb_function, bounds = list(x = xlimits, y = ylimits), coordinates = "cartesian")
            }
            
            
            
            pixelProb[cell] <- integralValue$value
            pixelVal[cell] <- singleBufferRaster.list[[i]][cellRowCol]
            
            
            
            
          }
          
          
        }
        
        singleCommVals.df <- data.frame(pixelProb, pixelVal)
        singleCommVals.df$weights <- pixelProb # / sum(pixelProb, na.rm = T)
        singleCommVals.df$DHSID <- rowDHSID
        
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
