count_buildings_osm <- function(latitudes, longitudes, radius = 1000) {
  if (length(latitudes) != length(longitudes)) {
    stop("latitudes and longitudes must have the same length")
  }
  
  results <- data.frame(
    latitude = latitudes,
    longitude = longitudes,
    n_buildings = NA_integer_
  )
  
  for (i in seq_along(latitudes)) {
    cat("Processing point", i, "of", length(latitudes), "...\n")
    
    point <- st_sfc(st_point(c(longitudes[i], latitudes[i])), crs = 4326)
    point_proj <- st_transform(point, 3857)
    buffer <- st_buffer(point_proj, radius)
    buffer_wgs84 <- st_transform(buffer, 4326)
    
    bbox <- st_bbox(buffer_wgs84)
    
    query <- opq(bbox = bbox) %>%
      add_osm_feature(key = "building")
    
    buildings <- tryCatch({
      osmdata_sf(query)$osm_polygons
    }, error = function(e) {
      warning(paste("Failed to fetch OSM data for point", i))
      return(NULL)
    })
    
    if (!is.null(buildings) && nrow(buildings) > 0) {
      buildings <- st_transform(buildings, st_crs(buffer_wgs84))
      in_buffer <- buildings[st_intersects(buildings, buffer_wgs84, sparse = FALSE), ]
      results$n_buildings[i] <- nrow(in_buffer)
    } else {
      results$n_buildings[i] <- 0
    }
  }
  
  return(results$n_buildings)
}



# Function to calculate total length of roads near each point
calculate_road_lengths <- function(latitudes, longitudes, radius = 1000) {
  if (length(latitudes) != length(longitudes)) {
    stop("Latitudes and longitudes must have the same length.")
  }
  
  results <- data.frame(Latitude = latitudes,
                        Longitude = longitudes,
                        RoadLength_m = NA_real_)
  
  for (i in seq_along(latitudes)) {
    cat("Processing point", i, "of", length(latitudes), "\n")
    
    # Create a point
    pt <- st_sfc(st_point(c(longitudes[i], latitudes[i])), crs = 4326)
    
    # Transform to metric projection (UTM based on lon)
    pt_utm <- st_transform(pt, crs = 32633)  # use a better UTM zone if needed
    
    # Buffer around the point
    buffer <- st_buffer(pt_utm, dist = radius)
    buffer_wgs <- st_transform(buffer, crs = 4326)  # back to lat/lon for OSM query
    
    # Query OSM roads in bounding box
    bbox <- st_bbox(buffer_wgs)
    
    roads <- tryCatch({
      opq(bbox = bbox) %>%
        add_osm_feature(key = "highway") %>%
        osmdata_sf()
    }, error = function(e) {
      warning(paste("OSM query failed at point", i, ":", e$message))
      return(NULL)
    })
    
    if (!is.null(roads) && !is.null(roads$osm_lines)) {
      road_lines <- st_transform(roads$osm_lines, crs = st_crs(buffer))
      road_lines_in_buffer <- st_intersection(road_lines, buffer)
      
      # Sum lengths
      total_length <- sum(st_length(road_lines_in_buffer))
      total_length <- units::set_units(total_length, "m")
      results$RoadLength_m[i] <- as.numeric(total_length)
    } else {
      results$RoadLength_m[i] <- NA
    }
  }
  
  return(results$RoadLength_m)
}


