###-------------------------------------------------###
###             BLOCK CROSS-VALIDATION              ###
###                                                 ###
###        Generate blocks for each cluster         ###
###-------------------------------------------------###

library(blockCV)


## Data downloading
baseline_30 <- brick("./inputs/envdata/baseline/processed/env_baseline_1deg.tif")
clust_sf <- readRDS("./outputs/bioregionalisation/clust_sf.RDS")



## 1. Search the optimal block size -----------------------------------------

# Non correlated variables
clust_var <- readRDS("./outputs/modelling/var_to_select.RDS")

blocksize <-  spatialAutoRange(rasterLayer = baseline_30[[clust_var[[1]]]],
                               sampleNumber = 5000, # number of cells to be used
                               doParallel = TRUE,
                               showPlots = TRUE)

plot(blocksize)
blocksize$range # optimal range 


## 2. Create folds -----------------------------------------
for (clust in unique(clust_sf$cluster)){ 
  print(clust)
  cur.blocks <- spatialBlock(speciesData = clust_sf, # Presence-absence data
                             species = paste0("cluster", clust), # Column with presence-absence information
                             rasterLayer = baseline_30[["temperature200"]], # Base map
                             theRange = blocksize$range, # Block size
                             k = 4, # Number of folds
                             selection = "random", # Method to assign blocks to folds 
                             iteration = 100, # Nb of runs to get optimal distribution into folds
                             biomod2Format = TRUE) # Create table in biomod2 format
  
  ## Save whole object
  # saveRDS(cur.blocks, paste0("./outputs/block_crossvalidation/cluster", clust, "_blocks_whole.RDS"))
  ## Save only fold table for biomod
  saveRDS(cur.blocks$biomodTable, paste0("./outputs/block_crossvalidation/cluster", clust, "_blocks.RDS"))
  
  ## Save map
  png(paste0("./outputs/plots/block_crossvalidation/blocks_cluster", clust, ".png"))
  plot(cur.blocks)
  dev.off()

}

