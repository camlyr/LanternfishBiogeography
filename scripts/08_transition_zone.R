###-------------------------------------###
###           TRANSITION ZONE           ###
###                                     ###
###  Based on Transition index formula  ###
###-------------------------------------###

library(raster)
library(tmap)
library(sf)


baseline_30 <- brick("./inputs/envdata/baseline/processed/env_baseline_1deg.tif")
land <- read_sf("./inputs/sigdata/ne_50m_land.shp")


fronts <- read_sf("./inputs/sigdata/ATLAS_Template/Fronts.shp")
pf <- read_sf("./inputs/sigdata/ATLAS_Template/PF.shp")
fronts <- fronts[!fronts$NAME %in% c("PF", "PF2", "SACCB"), ]  # Wrong paths for polar front in this shapefile
fronts <- fronts[!fronts$NAME %in% c("ICE_FEB", "ICE_OCT"), ]  # Remove ice cover limits
names(fronts) <- sub("NAME", "Fronts", names(fronts))  

palette_fronts <- c("orchid", "cyan3", "darkolivegreen3")

graticules <- st_graticule(ndiscr = 10000, 
                          lat = seq(-90, -30, 20),
                          lon = seq(0, 360, 30)) %>%
  st_transform("+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") %>% 
  st_geometry()


## 1. Present -----------------------------------------

# Load ensemble models for clusters 1 & 2
clust1 <- raster("./outputs/modelling/models/cluster1/cluster1/proj_baseline/proj_baseline_cluster1_em.grd")/1000
clust2 <- raster("./outputs/modelling/models/cluster2/cluster2/proj_baseline/proj_baseline_cluster2_em.grd")/1000

# Get coordinates where there is cluster information
coorXY <- xyFromCell(clust1, 1:ncell(clust1))
clust1_val  <- getValues(clust1)
clust2_val  <- getValues(clust2)
coorXY <- coorXY[-which(is.na(clust1_val)), ]

clust1_val  <- na.omit(clust1_val)
clust2_val  <- na.omit(clust2_val)

# Assign probabilities of cluster occurrence to each pixel
clust_proba <- data.frame(coorXY, clust1_val, clust2_val)

# Transition index formula
clust_proba$clust1_prop <- clust_proba$clust1_val/(clust_proba$clust1_val + clust_proba$clust2_val)
clust_proba$clust2_prop <- clust_proba$clust2_val/(clust_proba$clust1_val + clust_proba$clust2_val)
clust_proba$TI <- -(clust_proba$clust1_prop*log(clust_proba$clust1_prop) + clust_proba$clust2_prop*log(clust_proba$clust2_prop))              

clust_TI <- rasterize(x = clust_proba[, c("x", "y")], 
                           y = baseline_30,
                           field = clust_proba[, "TI"]) 
# Save
terra::writeRaster(clust_TI, "./outputs/transition_zone/clust_TI.tif")


# Maps
proj_transition <- tm_shape(clust_TI) +
  tm_raster(title = "Transition index") +
  tm_shape(land) +
  tm_borders(col = "grey") +
  tm_shape(graticules) +
  tm_lines(col = "grey", alpha = 0.5) +
  tm_layout(asp = 1, legend.outside = T)

tmap_save(tm = proj_transition,
          filename = paste0("./outputs/plots/transition_zone/transition_zone.png"),
          outer.margins = c(0.01, 0.01, 0.01, 0))

proj_transition_fronts <- tm_shape(clust_TI) +
  tm_raster(title = "Transition index") +
  tm_shape(land) +
  tm_borders(col = "grey") +
  tm_shape(fronts) +
  tm_lines(col = "Fronts", palette = palette_fronts, scale = 2, legend.col.show = F) +
  tm_shape(pf) +
  tm_lines(col = "firebrick1", scale = 2) +
  tm_shape(graticules) +
  tm_lines(col = "grey", alpha = 0.5) +
  tm_add_legend(type = c("line"), col = c("darkolivegreen3", "cyan3", "firebrick1", "orchid"), 
                labels = c("STF", "SAF", "PF", "SACCF"), title = "Fronts") +
  tm_layout(asp = 1, legend.outside = T)

tmap_save(tm = proj_transition_fronts,
          filename = paste0("./outputs/plots/transition_zone/transition_zone_fronts.png"),
          outer.margins = c(0.01, 0.01, 0.01, 0))


proj_transition01 <- tm_shape(clust_TI/maxValue(clust_TI)) +
  tm_raster(title = "Transition index \n(TI/TImax)") +
  tm_shape(land) +
  tm_borders(col = "grey") +
  tm_shape(graticules) +
  tm_lines(col = "grey", alpha = 0.5) +
  tm_layout(asp = 1, legend.outside = T)

tmap_save(tm = proj_transition01,
          filename = paste0("./outputs/plots/transition_zone/transition_zone_norm.png"),
          width = 2000, height = 1200,
          outer.margins = c(0.01, 0.01, 0.01, 0))


proj_transition01_fronts <- tm_shape(clust_TI/maxValue(clust_TI)) +
  tm_raster(title = "Transition index \n(TI/TImax)") +
  tm_shape(land) +
  tm_borders(col = "grey") +
  tm_shape(fronts) +
  tm_lines(col = "Fronts", palette = palette_fronts, scale = 2, legend.col.show = F) +
  tm_shape(pf) +
  tm_lines(col = "firebrick1", scale = 2) +
  tm_add_legend(type = c("line"), col = c("darkolivegreen3", "cyan3", "firebrick1", "orchid"), 
                labels = c("STF", "SAF", "PF", "SACCF"), title = "Fronts") +
  tm_shape(graticules) +
  tm_lines(col = "grey", alpha = 0.5) +
  tm_layout(asp = 1, legend.outside = T, main.title = "Transition zone", main.title.size = 1)

tmap_save(tm = proj_transition01_fronts,
          filename = paste0("./outputs/plots/transition_zone/transition_zone_norm_fronts.png"),
          width = 2000, height = 1200,
          outer.margins = c(0, 0.01, 0, 0))


## 2. Future CMIP6 -----------------------------------------

SSPs <- c("ssp126", "ssp245", "ssp585")
horizons <- c("2060", "2100")

for (ssp in SSPs) {
  
  for(horizon in horizons) {
    
    proj <- paste0(ssp, "_", horizon)
    print(proj)

    clust1 <- raster(paste0("./outputs/modelling/models/cluster1/cluster1/proj_", proj, "_cluster1_em.grd"))/1000
    clust2 <- raster(paste0("./outputs/modelling/models/cluster2/cluster2/proj_", proj, "_cluster2_em.grd"))/1000
    
    # Get coordinates where there is cluster information
    coorXY <- xyFromCell(clust1, 1:ncell(clust1))
    clust1_val  <- getValues(clust1)
    clust2_val  <- getValues(clust2)
    coorXY <- coorXY[-which(is.na(clust1_val)), ]
    
    clust1_val <-  na.omit(clust1_val)
    clust2_val <-  na.omit(clust2_val)
    
    # Assign probabilities of cluster occurrence to each pixel
    clust_proba <- data.frame(coorXY, clust1_val, clust2_val)
    
    # Transition index formula
    clust_proba$clust1_prop <- clust_proba$clust1_val/(clust_proba$clust1_val + clust_proba$clust2_val)
    clust_proba$clust2_prop <- clust_proba$clust2_val/(clust_proba$clust1_val + clust_proba$clust2_val)
    clust_proba$TI <- -(clust_proba$clust1_prop*log(clust_proba$clust1_prop) + clust_proba$clust2_prop*log(clust_proba$clust2_prop))              
    
    clust_TI <- rasterize(x = clust_proba[, c("x", "y")], 
                               y = baseline_30,
                               field = clust_proba[, "TI"]) 
    
    # Map TI/TImax
    proj_transition <- tm_shape(clust_TI/maxValue(clust_TI)) +
      tm_raster(title = "Transition index \n(TI/TImax)") +
      tm_shape(land) +
      tm_borders(col = "grey") +
      tm_shape(graticules) +
      tm_lines(col = "grey", alpha = 0.5) +
      tm_layout(main.title = paste0("Transition zone ", horizon, " (", ssp, ")"), asp = 1, legend.outside = TRUE)
    
    tmap_save(tm = proj_transition,
              filename = paste0("./outputs/plots/transition_zone/transition_zone_", proj,".png"),
              outer.margins = c(0, 0.01, 0, 0))
    
    # Save
    terra::writeRaster(clust_TI, paste0("./outputs/transition_zone/clust_TI_", proj, ".tif"))
    
  }
}



## 3. Present-future evolution CMIP6-----------------------------------------     

SSPs <- c("ssp126", "ssp245", "ssp585")
horizons <- c("2060", "2100")

palettediff <- colorRampPalette(colors = c("deepskyblue3", "lemonchiffon", "firebrick"))(12)

TI_present <- readRDS("./outputs/transition_zone/clust_TI.RDS")

for (ssp in SSPs) {
  
  for(horizon in horizons) {
    
    proj <- paste0(ssp, "_", horizon)
    print(proj)
    
    TI_future <- readRDS(paste0("./outputs/transition_zone/clust_TI_", proj, ".RDS"))
    
    TI_diff <- tm_shape((TI_future/maxValue(TI_future)) - (TI_present/maxValue(TI_present))) +
      tm_raster(title = "\u394Probabilities", palette = palettediff, 
                breaks = c(seq(-1, 0, length.out = 6), seq(0, 1, length.out = 6))) +
      tm_shape(land) +
      tm_borders(col = "grey") +
      tm_shape(graticules) +
      tm_lines(col = "grey", alpha = 0.5) +
      tm_layout(main.title = paste0("Difference (", horizon, " - present) ", ssp), asp = 1, legend.outside = TRUE)
    
    tmap_save(tm = TI_diff,
              filename = paste0("./outputs/plots/transition_zone/diff_", proj,"_transitionzone.png"),
              outer.margins = c(0, 0.01, 0, 0))
    
  }
}



## 4. Compare area of transition zone : present-future ----------------

library(raster)

# Option 1 - without Eastern Pacific

## All study area
present_TI <- raster("./outputs/transition_zone/clust_TI.tif")
future_TI <- raster("./outputs/transition_zone/clust_TI_ssp585_2100.tif")

present_TI <- present_TI/maxValue(present_TI)
future_TI <- future_TI/maxValue(future_TI)

lon_min_eastpacific <- -130
lon_max_eastpacific <- -65

coord <- xyFromCell(present_TI, 1:ncell(present_TI)) # Extract longitude and latitude values from raster's coordinates
pts_laea <- SpatialPoints(coord, proj4string = CRS("+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
pts_longlat <- spTransform(pts_laea, CRS("+proj=longlat +datum=WGS84")) # convert to decimal degrees
coord_longlat<- data.frame(coordinates(pts_longlat))
indices_to_remove <- which(coord_longlat$coords.x1 < lon_max_eastpacific & 
                             coord_longlat$coords.x1 > lon_min_eastpacific)

present_TI[indices_to_remove] <- NA # Set the corresponding pixels to NA
future_TI[indices_to_remove] <- NA 

## Compare area > threshold between present and future

threshold <- 0.8

# Present
filter_present <- present_TI < threshold
present <- mask(present_TI, filter_present, maskvalue = 1)
plot(present, main = "Present")
# plot(present_TI)
nbcell_present <- cellStats(present, function(i, ...) sum(!is.na(i)))

# Future
filter_future <- future_TI < threshold
future <- mask(future_TI, filter_future, maskvalue = 1)
plot(future, main = "Future")
# plot(future_TI)
nbcell_future <- cellStats(future, function(i, ...) sum(!is.na(i)))

diff_area_TI <- nbcell_future * 6242.4 / 10^6 - nbcell_present * 6242.4 / 10^6

cat("Transition zone (without Eastern Pacific) area difference:", diff_area_TI, "km2\n")



# Option 2 - only Eastern Pacific

## All study area
present_TI <- raster("./outputs/transition_zone/clust_TI.tif")
future_TI <- raster("./outputs/transition_zone/clust_TI_ssp585_2100.tif")

present_TI <- present_TI/maxValue(present_TI)
future_TI <- future_TI/maxValue(future_TI)

lon_min_eastpacific <- -130
lon_max_eastpacific <- -65

coord <- xyFromCell(present_TI, 1:ncell(present_TI)) # Extract longitude and latitude values from raster's coordinates
pts_laea <- SpatialPoints(coord, proj4string = CRS("+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
pts_longlat <- spTransform(pts_laea, CRS("+proj=longlat +datum=WGS84")) # convert to decimal degrees
coord_longlat<- data.frame(coordinates(pts_longlat))
indices_to_remove <- which(coord_longlat$coords.x1 > lon_max_eastpacific | 
                             coord_longlat$coords.x1 < lon_min_eastpacific)

present_TI[indices_to_remove] <- NA # Set the corresponding pixels to NA
future_TI[indices_to_remove] <- NA 


## Compare area > threshold between present and future

threshold <- 0.8

# Present
filter_present <- present_TI < threshold
present <- mask(present_TI, filter_present, maskvalue = 1)
plot(present, main = "Present")
# plot(present_TI)
nbcell_present <- cellStats(present, function(i, ...) sum(!is.na(i)))

# Future
filter_future <- future_TI < threshold
future <- mask(future_TI, filter_future, maskvalue = 1)
plot(future, main = "Future")
# plot(future_TI)
nbcell_future <- cellStats(future, function(i, ...) sum(!is.na(i)))

diff_area_TI <- nbcell_future * 6242.4 / 10^6 - nbcell_present * 6242.4 / 10^6

cat("Transition zone (only Eastern Pacific) area difference:", diff_area_TI, "km2\n")
