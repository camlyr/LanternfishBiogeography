###-------------------------------------------------###
###                   MODELLING                     ###
###                                                 ###
###       Species cluster distribution models       ###
###-------------------------------------------------###

library(biomod2)
library(raster)
library(tmap)
library(ggplot2)
library(viridis)
library(sf)

## Data downloading
baseline_30 <- brick("./inputs/envdata/baseline/processed/env_baseline_1deg.tif")

fronts <- read_sf("./inputs/sigdata/ATLAS_Template/Fronts.shp")
pf <- read_sf("./inputs/sigdata/ATLAS_Template/PF.shp")
fronts <- fronts[!fronts$NAME %in% c("PF", "PF2", "SACCB"), ]  # Wrong paths for polar front in this shapefile
fronts <- fronts[!fronts$NAME %in% c("ICE_FEB", "ICE_OCT"), ]  # Remove ice cover limits
names(fronts) <- sub("NAME", "Fronts", names(fronts))  

land <- read_sf("./inputs/sigdata/ne_50m_land.shp")

clust_colors <- readRDS("./outputs/bioregionalisation/clust_colors.RDS")

graticules <- st_graticule(ndiscr = 10000, 
                          lat = seq(-90, -30, 20),
                          lon = seq(0, 360, 30)) %>%
  st_transform("+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") %>% 
  st_geometry()



# --- # --- # --- # --- # --- # --- # --- # --- #
#
# /!\ Comment/uncomment the options below (1, 2a, 2b) to run the script with the desired parameters /!\
#
# --- # --- # --- # --- # --- # --- # --- # --- #


### Opt.1 Before variable selection ----

# clust_var <- readRDS("./outputs/modelling/var_to_select.RDS")
# file <- "var_selection"



### Opt.2 After variable selection ----

load("./outputs/modelling/selected_variables.RDS")
clust_var <- sel_var
file <- "models"



models <- c("GLM", "GAM", "ANN", "MARS", "FDA")

## Clusters to modelize
list_clust <- c("cluster1", "cluster2", "cluster3", "cluster4")
# list_clust <- c("cluster1", "cluster2")



initial.wd <-  getwd()


##  1. Model calibration -----------------------------------------

try(  # skip error if there is any
  for (clst in list_clust) {

    # Download occurrence data
    clust_occ <-  readRDS("./outputs/bioregionalisation/clust_df.RDS")

    # Selected variables
    cur_var <- clust_var[[clst]]
    clust_env_stack <- raster::stack(baseline_30[[cur_var]])

    # Occurrence formatting
    coorxy <- clust_occ[, c("x", "y")]
    obs <- clust_occ[, clst]

    # Create a specific file for the cluster
    if(!exists(paste0("./outputs/modelling/", file, "/", clst)))
    {
      dir.create(paste0("./outputs/modelling/", file, "/", clst), recursive = T)
    }
    setwd(paste0("./outputs/modelling/", file, "/", clst))

    # Biomod initialisation
    run_data <- BIOMOD_FormatingData(resp.var = obs,
                                     expl.var = clust_env_stack,
                                     resp.name = clst,
                                     resp.xy = coorxy,
                                     PA.nb.rep = 0, # nb of pseudoabsences
                                     PA.nb.absences = 0) # nb of pseudoabsence repetitions
    save(run_data, file = "run_data")

    # Overview of the results
    run_data
    str(run_data, max.level = 3)


    # Download blocks for cross-validation
    cur.folds <- readRDS(paste0(initial.wd, "/outputs/block_crossvalidation/", clst, "_blocks.RDS"))

    # Biomod calibration
    model_runs <- BIOMOD_Modeling(data = run_data,
                                  models =  models,
                                  DataSplitTable = cur.folds,  # folds for block cross-validation
                                  VarImport = 10,  # nb of permutation to estimate variable importance
                                  models.eval.meth = c("TSS", "ROC"),  # model evaluation metrics (for presence-absence data)
                                  SaveObj = F,
                                  Prevalence = NULL, # assign weight to presences
                                  Yweights = NULL,
                                  rescal.all.models = FALSE,
                                  do.full.models = FALSE)
    save(model_runs, file = "model_runs")



    # Overview of the results
    model_runs
    str(model_runs)
    get_variables_importance(model_runs)

    
    # Projection baseline
    projection_runs <- BIOMOD_Projection(modeling.output = model_runs, 
                                         new.env = clust_env_stack, # Environmental data used for the projection
                                         proj.name = "baseline",
                                         selected.models = 'all', # Which models to project
                                         binary.meth = NULL, # no binary transformation, keep continuous probability
                                         filtered.meth = NULL, # no threshold to assign 0 to a probability
                                         build.clamping.mask = TRUE, # Create a clamping mask to show regions where predictions are out of the data used for calibration
                                         # do.stack = F) # Save the projections (all calibration models and CV runs) in separate files
                                         do.stack = T) # Save the projections in a single stack
    save(projection_runs, file = paste0("./", clst, "/proj_baseline/projection_runs"))
    
    setwd(initial.wd)
  }
)


## 2. Model projection -----------------------------------------

for (clst in list_clust) {
  
  # Download projection data 
  baseline_cont <- stack(paste0("./outputs/modelling/", file, "/", clst, "/", clst, "/proj_baseline/proj_baseline_", clst, ".grd"))  # continuous projection on baseline
  
  # Use the cluster's color to create a palette
  palette <- viridis(6)
  
  # Title
  bioregion <- ifelse(clst == "cluster1", "Subtropical region",
                      ifelse(clst == "cluster2", "Southern region",
                             clst))
  
  # Map
  proj <- tm_shape(baseline_cont/1000) +
    tm_raster(title = "Probability of \nbioregion occurrence", 
              palette = palette) +
    tm_facets(ncol = length(names(baseline_cont))/4, scale.factor = 1.5) +  # one column per model (4 runs)
    tm_shape(land) +
    tm_borders(col = "grey") +
    tm_layout(main.title = paste("Individual projections of the", bioregion), main.title.size = 1, asp = 2, legend.outside.size = 0.15)
  
  
  tmap_save(tm = proj,
            filename = paste0("./outputs/plots/modelling/model_projections/proj_baseline_", clst, ".png"),
            width = 2000, height = 1000,
            outer.margins = c(0, 0.01, 0, 0))
  
  # Projection of clamping mask (where predictions are out of the data used for calibration)
  clamping_mask <- stack(paste0("./outputs/modelling/", file, "/", clst, "/", clst, "/proj_baseline/proj_baseline_ClampingMask.grd"))
  
  png(file = paste0("./outputs/plots/modelling/clamping_mask/proj_baseline_ClampingMask_", clst, ".png"), width = 2000, height = 1200)
  plot(clamping_mask, main = paste("Clamping Mask", clst), cex.main = 3, cex.axis = 2)
  plot(land, add = TRUE, border = "grey")
  dev.off()

}



## 3. Ensemble model -----------------------------------------


## Create ensemble model
for (clst in list_clust) {
  
  setwd(paste0("./outputs/modelling/", file, "/", clst, "/", clst, "/proj_baseline"))
  
  # Create ensemble model as the average of individual models
  model_stack <- stack(paste0("proj_baseline_", clst, ".grd"))
  cur_em <- mean(model_stack)

  # Uncertainty (standard deviation)
  cur_sd <- calc(model_stack, sd) 
  
  # Penalised probability of presence
  cur_em_penalised <- cur_em - cur_sd 
  cur_em_penalised[cur_em_penalised < 0] <- 0
  
  # Save
  writeRaster(cur_em, paste0("proj_baseline_", clst, "_em.grd"), overwrite = T)
  writeRaster(cur_sd, paste0("proj_baseline_", clst, "_sd.grd"), overwrite = T)
  writeRaster(cur_em_penalised, paste0("proj_baseline_", clst, "_em_pen.grd"), overwrite = T)
  
  setwd(initial.wd)
}


## Project ensemble model
for (clst in list_clust) {
  
  # Download projection data 
  baseline_em <- stack(paste0("./outputs/modelling/", file, "/", clst, "/", clst, "/proj_baseline/proj_baseline_", clst, "_em.grd"))
  baseline_em_pen <- stack(paste0("./outputs/modelling/", file, "/", clst, "/", clst, "/proj_baseline/proj_baseline_", clst, "_em_pen.grd"))
  baseline_sd <- stack(paste0("./outputs/modelling/", file, "/", clst, "/", clst, "/proj_baseline/proj_baseline_", clst, "_sd.grd"))
  
  # Use the cluster's color to create a palette
  palette <- viridis(6)
  palette_fronts <- c("orchid", "firebrick1", "darkorange")
  
  # Title
  bioregion <- ifelse(clst == "cluster1", "Subtropical region",
                      ifelse(clst == "cluster2", "Southern region",
                             clst))
  
  # Maps
  proj_em <- tm_shape(baseline_em/1000) +
    tm_raster(title = "Probability of \nbioregion occurrence",
              palette = palette, breaks = c(0, seq(0, 1, length.out = 6))) +
    # Continents
    tm_shape(land) +
    tm_borders(col = "grey") +
    # Fronts
    tm_shape(fronts) +
    tm_lines(col = "Fronts", palette = palette_fronts, scale = 2, legend.col.show = F) +
    tm_shape(pf) +
    tm_lines(col = "cyan1", scale = 2) +
    tm_add_legend(type = c("line"), col = c("darkorange", "firebrick1", "cyan1", "orchid"), 
                  labels = c("STF", "SAF", "PF", "SACCF"), title = "Fronts") +
    tm_shape(graticules) +
    tm_lines(col = "grey", alpha = 0.5) +
    tm_layout(main.title = paste(bioregion, "- Ensemble model"), main.title.size = 1, asp = 1, legend.outside = T)
  
  proj_em_pen <- tm_shape(baseline_em_pen/1000) +
    tm_raster(title = "Probability of \nbioregion occurrence", 
              palette = palette, breaks = c(0, seq(0, 1, length.out = 6))) +
    tm_shape(land) +
    tm_borders(col = "grey") +
    tm_shape(fronts) +
    tm_lines(col = "Fronts", palette = palette_fronts, scale = 2, legend.col.show = F) +
    tm_shape(pf) +
    tm_lines(col = "cyan1", scale = 2) +
    tm_add_legend(type = c("line"), col = c("darkorange", "firebrick1", "cyan1", "orchid"), 
                  labels = c("STF", "SAF", "PF", "SACCF"), title = "Fronts") +
    tm_layout(main.title = paste(bioregion, "- Penalised ensemble model"), main.title.size = 1, asp = 1, legend.outside = T)
  
  proj_sd <- tm_shape(baseline_sd/1000) +
    tm_raster(title = "Probability of \nbioregion occurrence", 
              palette = palette) +
    tm_shape(land) +
    tm_borders(col = "grey") +
    tm_shape(fronts) +
    tm_lines(col = "Fronts", palette = palette_fronts, scale = 2, legend.col.show = F) +
    tm_shape(pf) +
    tm_lines(col = "cyan1", scale = 2) +
    tm_add_legend(type = c("line"), col = c("darkorange", "firebrick1", "cyan1", "orchid"), 
                  labels = c("STF", "SAF", "PF", "SACCF"), title = "Fronts") +
    tm_layout(main.title = paste(bioregion, "- Standard deviation"), main.title.size = 1, asp = 1, legend.outside = T)
  
  
  # Save maps
  tmap_save(tm = proj_em,
            filename = paste0("./outputs/plots/modelling/model_projections/proj_baseline_", clst, "_em.png"),
            width = 2000, height = 1200,
            outer.margins = c(0, 0.01, 0, 0))

  tmap_save(tm = proj_em_pen,
            filename = paste0("./outputs/plots/modelling/model_projections/proj_baseline_", clst, "_em_pen.png"),
            width = 2000, height = 1200,
            outer.margins = c(0, 0.01, 0, 0))

  tmap_save(tm = proj_sd,
            filename = paste0("./outputs/plots/modelling/model_projections/proj_baseline_", clst, "_sd.png"),
            width = 2000, height = 1200,
            outer.margins = c(0, 0.01, 0, 0))
}
 
