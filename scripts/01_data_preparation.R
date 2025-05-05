###---------------------------------------------------###
###        PREPARE SPECIES & ENVIRONMENTAL DATA       ###
###                                                   ###
###       Load, compile, format the input data        ###
###---------------------------------------------------###


library(raster)
library(sf)
library(downloader)
library(ncdf4)
library(httr)
library(jsonlite)
library(xml2)
library(dplyr)
library(digest)
library(virtualspecies)
library(archive)



### A. BASELINE -----------------------------------------------------------------------

## 1. Download temperature, salinity, oxygen data ---------------------------------------

dir.create("./inputs/envdata/baseline", recursive = T)

# List of variables
variables <- data.frame(name = c("temperature", "salinity", "oxygen"),
                        abbrev = c("t", "s", "o"),
                        decade = c("decav", "decav", "all"))  # average of all available years (1955-2017 for t and s ; 1960-2017 for o)
# Seasons
periods <- 13:16 # 13 = South Hemisphere summer, 14 = North Hemisphere spring, 15 = NH summer, 16 = SH spring


for(i in 1:nrow(variables)) {
  
  for(period in periods) {
    
    file_name <- paste0("./inputs/envdata/baseline/woa18_",
                         variables$decade[i],
                         "_",
                         variables$abbrev[i], period,
                         "_01.nc")
    
    # Check if the file is already downloaded
    if (!file.exists(file_name)) {

      download.file(paste0("http://data.nodc.noaa.gov/woa/WOA18/DATA/",
                           variables$name[i],
                           "/netcdf/",
                           variables$decade[i],
                           "/1.00/woa18_",
                           variables$decade[i],
                           "_",
                           variables$abbrev[i], period,
                           "_01.nc"),  # resolution 1°
                    paste0("./inputs/envdata/baseline/woa18_",
                           variables$decade[i],
                           "_",
                           variables$abbrev[i], period,
                           "_01.nc"), mode = "wb")
    }
  }
}



## 2. Download chlorophyll a data ---------------------------------------


#### a. Query list of files ----

# List of variables to look for
variables <- c('chl')

# List of scenarios to look for
scenarios <- c('historical', 'ssp245')

# URL for ESGF queries
base_url <- 'https://esgf-node.llnl.gov/esg-search/search/'

# Common query parameters
common_params <- list(
  project = 'CMIP6',
  variant_label = 'r1i1p1f1',
  nominal_resolution = '1x1 degree',
  frequency = 'mon', # monthly
  type = 'File',
  format = 'application/solr+json',
  limit = 10000
)

# List to store results
results <- list()

for (var in variables) {
  results[[var]] <- list()  # For each variable
  for (scenario in scenarios) {
    cat(sprintf("Processing variable %s for scenario %s/n", var,
                scenario))
    
    # Add scenario parameter
    params <- c(common_params, list(variable_id = var, 
                                    experiment_id = scenario,
                                    activity_id = ifelse(scenario == "historical", 
                                                         "CMIP", 
                                                         "ScenarioMIP")))
    
    # Launch GET query
    response <- GET(base_url, query = params)
    
    # Check query status
    if (status_code(response) == 200) {
      data <- fromJSON(content(response, as = "text"))
      
      # Store the results for eahc variable and scenario
      results[[var]][[scenario]] <- data$response$docs
      
      cat(sprintf("End of processing for variable %s and scenario %s/n",
                  var, scenario))
    } else {
      warning(paste("Error in the query for variable", var,
                    "and scenario", scenario))
    }
  }
}


# Show results
for (var in names(results)) {
  for (scenario in names(results[[var]])) {
    cat(sprintf("/nVariable %s, Scenario %s:/n", var, scenario))
    cat(sprintf("Number of found files: %d/n",
                nrow(results[[var]][[scenario]])))
  }
}


#### b. Check available GCMs ----

# Empty data.frame to store the extracted data
cmip6_files <- data.frame(
  variable = character(),
  scenario = character(),
  gcm = character(),
  title = character(),
  data_node = character(),
  temporal_resolution = character(),
  projection_start = character(),
  projection_end = character(),
  spatial_resolution = character(),
  download_url = character(),
  checksum = character(),
  checksum_type = character(),
  stringsAsFactors = FALSE
)

gr <- function(x) x[1]

# Extract information for each variable and scenario
total_vars <- length(names(results))
for (var_index in seq_along(names(results))) {
  var <- names(results)[var_index]
  total_scenarios <- length(names(results[[var]]))
  
  for (scenario_index in seq_along(names(results[[var]]))) {
    scenario <- names(results[[var]])[scenario_index]
    
    # Show progress
    cat(sprintf("Processing variable %d/%d : %s, scenario %d/%d : %s/n", 
                var_index, total_vars, var, scenario_index, total_scenarios,
                scenario))
    
    # Extract information from the table in one go
    file_info <- results[[var]][[scenario]]
    gcm <- sapply(file_info$source_id, gr)
    temporal_resolution <- sapply(file_info$frequency, gr)
    title <- file_info$title
    spatial_resolution <- sapply(file_info$nominal_resolution, gr)
    download_url <- sapply(file_info$url, gr)  # Utiliser la première URL pour chaque ligne
    data_node <- file_info$data_node
    checksum = sapply(file_info$checksum, gr)
    checksum_type = sapply(file_info$checksum_type, gr)
    
    # Extract start and end dates from the title
    title_parts <- lapply(title, function(t) {
      date_part <- gsub("//.nc$", "", tail(strsplit(t, "_")[[1]], n = 1))
      strsplit(date_part, "-")[[1]]
    })
    
    # Divide dates into projection_start and projection_end
    projection_start <- sapply(title_parts,
                               function(x) x[1])
    projection_end <- sapply(title_parts,
                             function(x) x[2])
    
    # Create a data.frame for this variable and scenario
    df <- data.frame(
      variable = var,
      scenario = scenario,
      gcm = gcm,
      title = title,
      data_node = data_node,
      temporal_resolution = temporal_resolution,
      projection_start = projection_start,
      projection_end = projection_end,
      spatial_resolution = spatial_resolution,
      download_url = download_url,
      checksum = checksum,
      checksum_type = checksum_type,
      stringsAsFactors = FALSE
    )

    # Add the information to the final data.frame
    cmip6_files <- rbind(cmip6_files, df)
  }
}


#### c. Download files ----

library(httr)
library(dplyr)
library(digest)

# Filter the NetCDF files for the desired period, and available in historical and future
files_to_download <- cmip6_files %>%
  mutate(across(c(projection_start, projection_end), as.numeric)) %>%
  filter(projection_end > 195412) %>%  # Download starting 1955
  filter(projection_start < 201612) %>%  # Download until 2017
  mutate(across(c(projection_start, projection_end), as.character)) %>%
  group_by(gcm) %>%    # Group by climate model
  filter(all(c("historical", "ssp245") %in% scenario)) %>%  # Check that there are the 2 scenarios
  ungroup() 
gcm_list <- unique(files_to_download$gcm)

# FIle for download log
download_log_file <- "inputs/envdata/baseline/chl_cmip6_download_log.csv"

# Load or create a download log file
if (file.exists(download_log_file)) {
  download_log <- read.csv(download_log_file, stringsAsFactors = FALSE)
} else {
  download_log <- data.frame(
    title = character(),
    filename = character(),
    status = character(),
    stringsAsFactors = FALSE
  )
}

# File to save the characteristics of downloaded files
file_details_file <- "inputs/envdata/baseline/chl_cmip6_files.csv"

# Load or create a file to save the characteristics of the files
if (file.exists(file_details_file)) {
  file_details_log <- read.csv(file_details_file, stringsAsFactors = FALSE)
} else {
  file_details_log <- data.frame(
    variable = character(),
    scenario = character(),
    gcm = character(),
    temporal_resolution = character(),
    projection_start = character(),
    projection_end = character(),
    spatial_resolution = character(),
    title = character(),
    download_status = character(),
    stringsAsFactors = FALSE
  )
}

# Function to clean the URL and get the NetCDF file URL only
clean_url <- function(url) {
  # Remove everything after the first '|' (characteristic of URL downloads)
  clean <- strsplit(url, "//|")[[1]][1]
  return(clean)
}


# Function to verify the checksum of a file
verify_checksum <- function(destfile, expected_checksum) {
  file_checksum <- digest(destfile, algo = "sha256", file = TRUE)
  return(file_checksum == expected_checksum)
}

# Function to download a file with checksum verification
download_file <- function(url, destfile, expected_checksum, max_retries = 3) {
  clean_url <- clean_url(url)
  
  # Try to download up to max_retries times if the checksum fails
  for (attempt in 1:max_retries) {
    tryCatch({
      # Download file
      response <- GET(clean_url, write_disk(destfile, overwrite = TRUE), progress())
      
      # Check if the download was successful
      if (status_code(response) == 200) {
        # Verify checksum
        if (verify_checksum(destfile, expected_checksum)) {
          message("Successful downloading and checksum verification for : ", destfile)
          return(TRUE)  # Success
        } else {
          message("Incorrect checksum, attempt ", attempt, " for : ", destfile)
        }
      } else {
        message("Failed downloading for URL : ", clean_url)
      }
    }, error = function(e) {
      message("Error in downloading ", destfile, ": ", e)
    })
    
    # If we reach here, it means either the download or the checksum failed
    if (attempt < max_retries) {
      message("New downloading attempt (", attempt + 1, "/", max_retries, ") for : ", destfile)
    }
  }
  
  # If all attempts fail
  message("Impossible to download the file correctly after ", max_retries, " attempts : ", destfile)
  
  return(FALSE)
}

unique_titles <- unique(files_to_download$title)

for (title in unique_titles) {
  # Filter the rows corresponding to this title
  file_info_subset <- files_to_download %>% filter(title == !!title)
  
  # Create the destination path for the file
  filename <- title
  destfile <- file.path("inputs/envdata/baseline", filename)
  
  # Check if the file already exists or is marked as successful in the log
  if (file.exists(destfile) || title %in% download_log$title[download_log$status == "success"]) {
    message("File is already downloaded or marked as successful : ", filename)
    next
  }
  
  # Try to download from each data node
  download_success <- FALSE
  expected_checksum <- file_info_subset$checksum[1]  # Extract the expected checksum
  for (i in 1:nrow(file_info_subset)) {
    url <- file_info_subset$download_url[i]
    message("Downloading attempt from : ", url)
    if (download_file(url, destfile, expected_checksum)) {
      download_success <- TRUE
      break
    }
  }
  
  # Extract file characteristics
  variable <- file_info_subset$variable[1]
  scenario <- file_info_subset$scenario[1]
  gcm <- file_info_subset$gcm[1]
  temporal_resolution <- file_info_subset$temporal_resolution[1]
  projection_start <- file_info_subset$projection_start[1]
  projection_end <- file_info_subset$projection_end[1]
  spatial_resolution <- file_info_subset$spatial_resolution[1]
  
  # Update checksum status
  checksum_status <- ifelse(download_success, "checksum_pass", "checksum_fail")
  
  # Update the file characteristics log
  file_details_log <- rbind(file_details_log, data.frame(
    variable = variable,
    scenario = scenario,
    gcm = gcm,
    temporal_resolution = temporal_resolution,
    projection_start = projection_start,
    projection_end = projection_end,
    spatial_resolution = spatial_resolution,
    title = title,
    download_status = ifelse(download_success, "success", "failed"),
    checksum_status = checksum_status,
    stringsAsFactors = FALSE
  ))
  
  # Update download log
  download_log <- rbind(download_log, data.frame(
    title = title,
    filename = filename,
    status = ifelse(download_success, "success", "failed"),
    checksum_status = checksum_status,
    stringsAsFactors = FALSE
  ))
  
  # Save the log and log files
  write.csv(download_log, download_log_file, row.names = FALSE)
  write.csv(file_details_log, file_details_file, row.names = FALSE)
  
  if (!download_success) {
    message("Impossible to download file : ", filename)
  }
}



## 3. Load bathymetry data ---------------------------------------

bathy <- raster("./inputs/envdata/baseline/topo30_raw.grd")

# Adapt to a 1° x 1° resolution grid
bathy <- aggregate(bathy, fact = 120)
# Extent -180,180 instead of 0,360
bathy <- rotate(bathy)

# Save file
saveRDS(bathy, "./inputs/envdata/baseline/bathy.RDS")
writeRaster(bathy, "./inputs/envdata/baseline/bathy.tif", format = "GTiff", overwrite = TRUE)



## 4. Create the baseline ---------------------------------------

#### a. Temperature, salinity, oxygen -------------

variables <- data.frame(name = c("temperature", "salinity", "oxygen"),
                        abbrev = c("t", "s", "o"),
                        decade = c("decav", "decav", "all"))  # average of all available years (1955-2017 for t and s ; 1960-2017 for o)

depthlevels <- data.frame(level = c(1, 21, 25, 37), # see WOA18 documentation for corresponding levels
                          depth = c(0, 100, 200, 500))

# Assign Spring-Summer for each hemisphere
hemisphere <- data.frame(hemisphere = c("north", "south"),
                         s1 = c(14, 13),  # season 13 = North Hemisphere winter, 14 = NH spring
                         s2 = c(15, 16))  # 15 = NH summer, 16 = NH autumn

baseline <- list(north = stack(),
                 south = stack())

for (hemis in 1:nrow(hemisphere))
{
  for (var in 1:nrow(variables))
  {
    for (depth in 1:nrow(depthlevels))
    {
      cur.s1 <- raster(paste0("./inputs/envdata/baseline/woa18_", 
                              variables$decade[var], "_", 
                              variables$abbrev[var], 
                              hemisphere$s1[hemis], "_01.nc"), 
                       varname = paste0(variables$abbrev[var], "_an"),
                       lvar = 4,
                       band = depthlevels$level[depth])
      cur.s2 <- raster(paste0("./inputs/envdata/baseline/woa18_", 
                              variables$decade[var], "_",  
                              variables$abbrev[var], 
                              hemisphere$s2[hemis], "_01.nc"),
                       varname = paste0(variables$abbrev[var], "_an"),
                       lvar = 4,
                       band = depthlevels$level[depth])
      cur.var <- mean(cur.s1, cur.s2)
      baseline[[hemis]] <- stack(baseline[[hemis]],
                                 cur.var)
      names(baseline[[hemis]])[nlayers(baseline[[hemis]])] <- 
        paste0(variables$name[var],
               depthlevels$depth[depth])
    }
  }
}


#### b. Bathymetry -------------

bathy <- readRDS("./inputs/envdata/baseline/bathy.RDS")

baseline[["north"]] <- stack(baseline[["north"]],
                             bathy)
baseline[["south"]] <- stack(baseline[["south"]],
                             bathy)

names(baseline[["north"]])[nlayers(baseline[["north"]])] <- 
  names(baseline[["south"]])[nlayers(baseline[["south"]])] <- "bathymetry"


#### c. Chlorophyll a -------------

gcm_list <- c("GFDL-CM4", "GFDL-ESM4")

chla_gcm_south_list <- list()
chla_gcm_north_list <- list()

depthlevels <- data.frame(level = c(1, 7, 10, 14),  # Level 1 for 2.5m, level 10 for 200m
                          depth = c(2.5, 100, 200, 500))

for (gcm in gcm_list) {
  
  cat(sprintf("Processing model %s/n", gcm))
  
  raster_list <- list.files("inputs/envdata/baseline", full.names = TRUE, pattern = paste0("chl_Omon_", gcm))
  
  chla_gcm_south_list_depth <- list()  
  chla_gcm_north_list_depth <- list()
  
  # Load and filter by depth and period for each depth
  for (depth in 1:nrow(depthlevels)) {
    
    cat(sprintf("Processing depth %f m/n", depthlevels$depth[depth]))
    
    # Create a RasterStack for the GCM at the given depth
    chla_gcm <- stack(lapply(raster_list, function(raster) {
      brick(raster, level = depthlevels$level[depth])  # Filter on given depth
    }))
    
    
    # Cut to period 1955-2017
    layer_names <- names(chla_gcm) # Extract layer names
    dates <- as.Date(gsub("X", "", layer_names), format="%Y.%m.%d") # Convert to dates (format YYYY.MM.DD)
    start_date <- as.Date("1955-01-01")
    end_date <- as.Date("2017-12-31")
    selected_layers <- which(dates >= start_date & dates <= end_date) # Find the indices of the layers corresponding to the period
    chla_gcm <- stack(chla_gcm[[selected_layers]]) # Create new RasterStack with those layers
    
    
    # Spring-Summer season
    dates_1955_2017 <- dates[selected_layers]  # Available dates
    months <- as.numeric(format(dates_1955_2017, "%m"))  # Convert to numerical months
    selected_months_south <- which(months %in% c(10, 11, 12, 1, 2, 3)) # Filter the indices corresponding to October-March for the Southern Hemisphere
    selected_months_north <- which(months %in% c(4, 5, 6, 7, 8, 9)) # Filter the indices corresponding to April-September for the Northern Hemisphere
    chla_gcm_south <- stack(chla_gcm[[selected_months_south]])  # Create new RasterStack with those layers
    chla_gcm_north <- stack(chla_gcm[[selected_months_north]])  # Create new RasterStack with those layers
    
    
    # Average value
    chla_gcm_final_south <- mean(chla_gcm_south)
    chla_gcm_final_north <- mean(chla_gcm_north)
    
    # Add to the lists for both hemispheres
    chla_gcm_south_list_depth[[as.character(depthlevels$depth[depth])]] <- chla_gcm_final_south
    chla_gcm_north_list_depth[[as.character(depthlevels$depth[depth])]] <- chla_gcm_final_north
  }
  
  # Add to the lists for both hemispheres
  chla_gcm_south_list[[gcm]] <- chla_gcm_south_list_depth
  chla_gcm_north_list[[gcm]] <- chla_gcm_north_list_depth
}

# Average over GCMs for each depth level
chla_south_2.5 <- rotate(mean(stack(lapply(chla_gcm_south_list, function(x) x$'2.5'))))
chla_north_2.5 <- rotate(mean(stack(lapply(chla_gcm_north_list, function(x) x$'2.5'))))

chla_south_100 <- rotate(mean(stack(lapply(chla_gcm_south_list, function(x) x$'100'))))
chla_north_100 <- rotate(mean(stack(lapply(chla_gcm_north_list, function(x) x$'100'))))

chla_south_200 <- rotate(mean(stack(lapply(chla_gcm_south_list, function(x) x$'200'))))
chla_north_200 <- rotate(mean(stack(lapply(chla_gcm_north_list, function(x) x$'200'))))

chla_south_500 <- rotate(mean(stack(lapply(chla_gcm_south_list, function(x) x$'500'))))
chla_north_500 <- rotate(mean(stack(lapply(chla_gcm_north_list, function(x) x$'500'))))


# Rename layers
names(chla_south_2.5) <- names(chla_north_2.5) <- "chla2.5"
names(chla_south_100) <- names(chla_north_100) <- "chla100"
names(chla_south_200) <- names(chla_north_200) <- "chla200"
names(chla_south_500) <- names(chla_north_500) <- "chla500"


# Add to baseline (for the 2 depths)
baseline[["south"]] <- addLayer(baseline[["south"]],
                                extend(chla_south_2.5, baseline[["south"]])) 
baseline[["north"]] <- addLayer(baseline[["north"]],
                                extend(chla_north_2.5, baseline[["north"]])) 
baseline[["south"]] <- addLayer(baseline[["south"]],
                                extend(chla_south_100, baseline[["south"]])) 
baseline[["north"]] <- addLayer(baseline[["north"]],
                                extend(chla_north_100, baseline[["north"]])) 
baseline[["south"]] <- addLayer(baseline[["south"]],
                                extend(chla_south_200, baseline[["south"]])) 
baseline[["north"]] <- addLayer(baseline[["north"]],
                                extend(chla_north_200, baseline[["north"]])) 
baseline[["south"]] <- addLayer(baseline[["south"]],
                                extend(chla_south_500, baseline[["south"]])) 
baseline[["north"]] <- addLayer(baseline[["north"]],
                                extend(chla_north_500, baseline[["north"]])) 


## 4. Differentiate the hemispheres ---------------------------------------

## Crop data for each hemisphere
baseline[["north"]] <- crop(baseline[["north"]],
                            extent(-180, 180, 0, 90)) # Northern Hemisphere
baseline[["south"]] <- crop(baseline[["south"]],
                            extent(-180, 180, -90, 0)) # Southern Hemisphere
baseline[["south30"]] <- crop(baseline[["south"]],
                            extent(-180, 180, -90, -30)) # Southern Hemisphere below 30°S

## Set the same initial projection for the whole baseline
crs(baseline[["south30"]]) <- crs(baseline[["south"]]) <- crs(baseline[["north"]]) <-  crs(bathy)


## Synchronize NAs among layers
# Create a mask that marks the cells containing NAs in all layers
na_mask_north <- calc(baseline[["north"]], fun = function(x) {
  return(any(is.na(x)))     # If a layer contains NA, then we mark the whole cell as NA
})
na_mask_south <- calc(baseline[["south"]], fun = function(x) {
  return(any(is.na(x))) 
})
na_mask_south30 <- calc(baseline[["south30"]], fun = function(x) {
  return(any(is.na(x)))
})

# Apply the mask to all layers of the RasterBrick
for (i in 1:nlayers(test)) {
  values(baseline[["north"]][[i]])[values(na_mask_north) == 1] <- NA
}
for (i in 1:nlayers(test)) {
  values(baseline[["south"]][[i]])[values(na_mask_south) == 1] <- NA
}
for (i in 1:nlayers(test)) {
  values(baseline[["south30"]][[i]])[values(na_mask_south30) == 1] <- NA
}



## Merge the two hemispheres into full baseline
fullbaseline <- do.call(merge, list(baseline$north, baseline$south))
names(fullbaseline) <- names(baseline$north)

## Save the whole world baseline
saveRDS(fullbaseline, "./inputs/envdata/baseline/processed/fullbaseline.RDS")

## Save -30° baseline (lon-lat version)
saveRDS(baseline[["south30"]], "./inputs/envdata/baseline/processed/baseline_30_lonlat.RDS")


## Pixel size in long-lat projection
plot(area(baseline[["south"]]))
area(baseline[["south"]])



## 5. Polar projection ---------------------------------------

## Lambert azimuthal equal-area in order to have same size pixels
baseline[["north"]] <- projectRaster(baseline[["north"]],
                                     crs = "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
baseline[["south"]] <- projectRaster(baseline[["south"]],
                                     crs = "+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
baseline[["south30"]] <- projectRaster(baseline[["south30"]],
                                     crs = "+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")


## Pixel size in polar LAEA projection
plot(area(baseline[["south"]]))
area(baseline[["south"]])

## Save the whole world baseline in polar projection
saveRDS(baseline, "./inputs/envdata/baseline/processed/baseline.RDS")



## 6. Study region below -30° latitude ---------------------------------------

# Save as tif using terra package
baseline_30 <- rast(baseline[["south30"]])
terra::writeRaster(baseline_30, filename = "./inputs/envdata/baseline/processed/env_baseline_1deg.tif", overwrite = TRUE)



### B. SPECIES -----------------------------------------------------------------------

## 1. Load data below -30° latitude ---------------------------------------

fulldb <- readRDS("./inputs/occdata/processed/fulldb_after_checking.RDS")
baseline <- readRDS("./inputs/envdata/baseline/processed/baseline_30.RDS")

class(fulldb) <- "data.frame"   # Raster functions don't work with class "spatialvalid" objects

## Study area = below -30° latitude 
sp_occ <- fulldb[fulldb$decimalLatitude <= -30, ] 


## 2. Transform coordinates for polar LAEA projection ---------------------------------------

## Convert dataframe to sf object
sp_occ_sf <- st_as_sf(x = sp_occ,
                      coords = c("decimalLongitude", "decimalLatitude"),
                      crs = "+proj=longlat +ellps=WGS84")

## Convert coordinates
sp_occ_sf <- st_transform(x = sp_occ_sf,
                          crs = CRS("+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))


## Replace spaces by points in species name
sp_occ_sf$scientificName <- gsub(" ", ".",
                                 sp_occ_sf$scientificName, fixed = TRUE) # replace spaces by points

## Aggregate records on a 1° grid 
sp_stack <- raster(baseline)

for (sp in unique(sp_occ_sf$scientificName)) {
  sp_stack <- addLayer(sp_stack,
                       rasterize(x = sp_occ_sf[sp_occ_sf$scientificName %in% sp, ],
                                 y = baseline,
                                 field = "scientificName",
                                 fun = function (x, ...) 1))  # Assign "1" to a cell everytime one or more occurrences fall in this cell
}

names(sp_stack) <- unique(sp_occ_sf$scientificName)


# Save as tif using terra package
sp_stack <- rast(sp_stack)
terra::writeRaster(sp_stack, filename = "./inputs/occdata/processed/mycto_occ_1deg.tif", overwrite = TRUE)



## 2. Save as dataframe with lat/lon coordinates ---------------------------------------

# Extraire les coordonnées des points
coords <- st_coordinates(sp_occ_sf)

# Obtenir les indices des cellules du raster pour chaque point
cell_ids <- cellFromXY(baseline, coords)

# Ajouter les indices au dataframe original
sp_occ_sf$cell_id <- cell_ids

# Obtenir les coordonnées du centre de chaque cellule
cell_coords <- xyFromCell(baseline, cell_ids)

# Ajouter les coordonnées au dataframe original
sp_occ_sf$cell_coord <- cell_coords

# Projection
sp_occ_latlon <- st_transform(sp_occ_sf, crs = 4326)

coords <- st_coordinates(sp_occ_latlon)
sp_occ_latlon$decimalLongitude <- coords[, 1]
sp_occ_latlon$decimalLatitude <- coords[, 2]

sp_occ_latlon <- sp_occ_latlon %>% 
  as.data.frame() %>%
  dplyr::select(scientificName, decimalLongitude, decimalLatitude, everything(),
                -c(geometry, cell_coord, cell_id, id)) %>%
  mutate(source = ifelse(source == "MyctoDB", "Atlas", source))

# Save as CSV file
write.csv2(southdb_latlon, "./inputs/occdata/processed/mycto_occ_1deg.csv", row.names = FALSE)



### C. LAND BASE MAP -----------------------------------------------------------------------

land <- read_sf("./inputs/sigdata/ne_50m_land.shp")

## Polar LAEA projection
landpolar <- st_transform(x = land,
                          crs = CRS("+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))

saveRDS(landpolar, "./inputs/sigdata/landpolar.RDS")






### D. CMIP6 FUTURE PROJECTIONS -----------------------------------------------------------------------

GCMs <- c("AWI-CM-1-1-MR",
          "FGOALS-g3",
          "CanESM5",
          "CMCC-ESM2",
          "ACCESS-CM2",
          "FIO-ESM-2-0",
          "IPSL-CM6A-LR",
          "MPI-ESM1-2-LR",
          "MRI-ESM2-0",
          "CIESM")
SSPs <- c("ssp126", "ssp245", "ssp585")
horizons <- c("2060", "2100")


hist_stack <- stack()

for(GCM in GCMs) {
  print(GCM)
  
  ## 1. Historical data ----
  
  hist_file <- paste0("./inputs/envdata/future/raw/", GCM, "_hist_final.nc")
  
  if (!file.exists(hist_file)) {
    # Unzip archived file
    archive_extract(archive = paste0("./inputs/envdata/future/raw/", GCM, "_hist_final.7z"), 
                    dir = paste0("./inputs/envdata/future/raw/"))
  }
  
  # Depth = 200 m (level 10)
  hist <- brick(hist_file,
                varname = "thetao", level = 10)
  # Seasons = October-March
  season <- c("10", "11", "12", "01", "02", "03") 
  hist_season <- hist[[which(format(as.Date(getZ(hist), format = "X%Y.%m.%d"), format = "%m") %in% season)]]
  # Time period = 1950-2014
  hist_year <- hist_season[[which(format(as.Date(names(hist_season), format = "X%Y.%m.%d")) >= as.Date("1950-01-01") & 
                                    format(as.Date(names(hist_season), format = "X%Y.%m.%d")) <= as.Date("2014-12-31"))]]
  # Values = average on the time period
  hist_mean <- mean(hist_year, na.rm = TRUE)
  
  
  ## 2. Future data ----
  
  for (ssp in SSPs) {
    print(ssp)
    
    future_file <- paste0("./inputs/envdata/future/raw/", GCM, "_", ssp, "_final.nc")
    
    if (!file.exists(future_file)) {
      # Unzip archived file
      archive_extract(archive = paste0("./inputs/envdata/future/raw/", GCM, "_", ssp, "_final.7z"),
                      dir = paste0("./inputs/envdata/future/raw/"))
    }
    
    # Depth = 200 m (level 10)
    absfuture <- brick(future_file, 
                       varname = "thetao", level = 10)
    # Seasons = October-March
    season <- c("10", "11", "12", "01", "02", "03") 
    absfuture_season <- absfuture[[which(format(as.Date(getZ(absfuture), format = "X%Y.%m.%d"), format = "%m") %in% season)]]
    
    for (horizon in horizons) {
      print(horizon)
      
      if (!(horizon == "2100" & ssp == "ssp245" & GCM == "FGOALS-g3") &       # Missing values for these time periods
          !(horizon == "2100" & ssp == "ssp126" & GCM == "IPSL-CM6A-LR")) {
        
        # Select time period
        if (horizon == "2060"){
          print("2060")
          absfuture_year <- absfuture_season[[which(format(as.Date(names(absfuture_season), format = "X%Y.%m.%d")) >= as.Date("2041-01-01") & 
                                                      format(as.Date(names(absfuture_season), format = "X%Y.%m.%d")) <= as.Date("2060-12-31"))]]
        }
        if (horizon == "2100") {  
          print("2100") 
          absfuture_year <- absfuture_season[[which(format(as.Date(names(absfuture_season), format = "X%Y.%m.%d")) >= as.Date("2081-01-01") & 
                                                      format(as.Date(names(absfuture_season), format = "X%Y.%m.%d")) <= as.Date("2100-12-31"))]]
        }
        
        
        # Values = average on the time period
        absfuture_mean <- mean(absfuture_year, na.rm = TRUE)
        # Anomaly
        anomaly <- absfuture_mean - hist_mean
        # Set the same extent, origin, projection as the baseline
        anomaly <- projectRaster(rotate(anomaly), baseline_LL, method = "bilinear")
        anomaly <- projectRaster(anomaly, baseline_30, method = "bilinear")
        
        
        # Sum the anomaly and the baseline to get future values
        future <- subset(baseline_30, "temperature200") + anomaly
        names(future) <- "temperature200"
        
        # Save as tif using terra package
        future <- rast(future)
        terra::writeRaster(future, filename = paste0("./inputs/envdata/future/processed/temp200_", ssp, "_", horizon,"_", GCM, "_1deg.tif"), overwrite = TRUE)
        
      }
    }
  }
}

