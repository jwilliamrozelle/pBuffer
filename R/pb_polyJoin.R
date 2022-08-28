#' @name pb_polyJoin
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
#'





pb_polyJoin <- function(displaced.sf, polygons.sf, displaced.id = "DHSID", densityBuffer = densitybuffer, adminBound = haiti_adm2, adminID = "ID_2", n.cores = 1) {
  
  # Check to make sure that objects which should be in sf format are in sf format
  #   displaced.sf
  if (!"sf" %in% class(displaced.sf)) {
    stop("The specified displaced.sf is not in sf format. It must be in sf format for this function to work.")
  }
  #   polygons.sf
  if (!"sf" %in% class(polygons.sf)) {
    stop("The specified polygons.sf is not in sf format. It must be in sf format for this function to work.")
  }
  #   adminBound
  if (!"sf" %in% class(displaced.sf)) {
    stop("The specified adminBound is not in sf format. It must be in sf format for this function to work.")
  }
  
  
  # ensure that the adminBound sf object contains the unique ID
  if (!adminID %in% names(adminBound)) {
    stop(paste0(adminID, " is not found in the adminBound sf object."))
  }
  
  if (!is.null(adminBound)) {
    if (sum(duplicated(adminBound[[adminID]])) > 0) {
      stop("adminID must specify a uniquely valid ID for the adminbound object")
    }
  }
  
  # ensure that the displaced.sf sf object contains the unique ID
  if (!displaced.id %in% names(displaced.sf)) {
    stop(paste0(displaced.id, " is not found in the displaced.sf sf object."))
  }
  
  if (sum(duplicated(displaced.sf[[displaced.id]])) > 0) {
    stop("displaced.id must specify a uniquely valid ID for the displaced.sf object")
  }
  
  # rename the adminboundary unique ID to adminID
  adminBound$adminID <- adminBound[[adminID]]
  
  
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
        
        # If the administrative boundary is defined
        if (!is.null(adminBound)) {
          # get the single admin boundary for a community
          singleAdminBound <- suppressWarnings(sf::st_intersection(singleComm, adminBound))
          singleAdminBound.poly <- dplyr::filter(adminBound, adminID == singleAdminBound$adminID[1])
        }
        if (st_coordinates(singleComm)[2] != 0) {
          
          # rasterize the displacement buffer around the single community
          singleDens.raster <- rasterizeDisplacement(
            densitybuffer, 
            initialCoords = st_coordinates(singleComm), 
            inputCRS = crs2use
          )
          
          
          # If the administrative boundary is specified, trim the probability buffer by the polygon layer
          if(!is.null(adminBound)) {
            # # Trim the weighted raster !!! 
            # singleDens.raster <- trimProbBuff(singleDens.raster, adminBound = adminBound)
            
            # mask based on the administrative bounadries
            singleDens.raster <- raster::mask(singleDens.raster, singleAdminBound.poly)
          }
          
          
          # turn these into a point object
          singleDens.sf <- st_rasterAsPoint(singleDens.raster)
          rm(singleDens.raster)
          
          # Get the most probable nearby health facility
          probableHF <- ml_polygon(singleDens.sf, polygons.sf)
          
          # put dhsid, spaid, and weight into a merged, single-row data frame
          results.df <- data.frame(DHSID = rowDHSID)
          results.df <- cbind(results.df, probableHF)
        } else {
          results.df <- NA
        }
        
        return(results.df)
        
      })
    },
    
    # always stop the cluster, even if future_apply fails
    finally = {parallel::stopCluster(cl)}
    )
    
  } else {
    
    closestHFs <- apply(st_drop_geometry(displaced.sf), 1, function(x) {
      
      rowDHSID <- x[displaced.id]
      
      # extract an sf object for each row of the data frame
      singleComm <- filter(displaced.sf, DHSID == rowDHSID)
      crs2use <- raster::crs(singleComm)
      
      # If the administrative boundary is defined
      if (!is.null(adminBound)) {
        # get the single admin boundary for a community
        singleAdminBound <- suppressWarnings(sf::st_intersection(singleComm, adminBound))
        singleAdminBound.poly <- dplyr::filter(adminBound, adminID == singleAdminBound$adminID[1])
      }
      if (st_coordinates(singleComm)[2] != 0) {
        
        # rasterize the displacement buffer around the single community
        singleDens.raster <- rasterizeDisplacement(
          densitybuffer, 
          initialCoords = st_coordinates(singleComm), 
          inputCRS = crs2use
        )
        
        
        # If the administrative boundary is specified, trim the probability buffer by the polygon layer
        if(!is.null(adminBound)) {
          # Trim the weighted raster
          singleDens.raster <- trimProbBuff(singleDens.raster, adminBound = adminBound)
          
          # mask based on the administrative bounadries
          singleDens.raster <- raster::mask(singleDens.raster, singleAdminBound.poly)
        }
        
        
        # turn these into a point object
        singleDens.sf <- st_rasterAsPoint(singleDens.raster)
        rm(singleDens.raster)
        
        # Get the most probable nearby health facility
        probableHF <- ml_polygon(singleDens.sf, polygons.sf)
        
        # put dhsid, spaid, and weight into a merged, single-row data frame
        results.df <- data.frame(DHSID = rowDHSID)
        results.df <- cbind(results.df, probableHF)
      } else {
        results.df <- NA
      }
      
      return(results.df)
      
    })
  }
  
  closestHFs.df <- iotools::fdrbind(closestHFs) |> as.data.frame()
  
  return(closestHFs.df)
}
