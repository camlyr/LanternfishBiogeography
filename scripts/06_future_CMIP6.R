###-------------------------------------###
###             MODELLING               ###
###                                     ###
###         Future projections          ###
###-------------------------------------###


library(ncdf4)
library(raster)
library(biomod2)
library(tmap)
library(ggplot2)
library(sf)
library(downloader)
library(virtualspecies)
library(dplyr)
library(viridis)
library(rstatix)

initial.wd <- getwd()

baseline_30 <- brick("./inputs/envdata/baseline/processed/env_baseline_1deg.tif")
baseline_LL <- readRDS("./inputs/envdata/baseline/processed/baseline_30_lonlat.RDS")

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

land <- read_sf("./inputs/sigdata/ne_50m_land.shp")


## 1. Model projection -----------------------------------------     

clust_colors <- readRDS("./outputs/bioregionalisation/clust_colors.RDS")

graticules = st_graticule(ndiscr = 10000, 
                          lat = seq(-90, -30, 20),
                          lon = seq(0, 360, 30)) %>%
  st_transform("+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") %>% 
  st_geometry()

models <- c('GLM', "GAM", "ANN", "MARS", "GBM", "FDA")

list_clust <- c("cluster1", "cluster2")

file <- "models"

listStacks <- expand.grid(SSPs, horizons, GCMs)
listStacks <- apply(listStacks, 1,
                    function(x) paste0(x[1], "_", x[2], "_", x[3] ))
listStacks <- listStacks[! listStacks %in% c("ssp245_2100_FGOALS-g3", "ssp126_2100_IPSL-CM6A-LR")]  # Missing climatic data

initial.wd <- getwd()

for (clst in list_clust) {
  
  load(paste0("./outputs/modelling/models/", clst, "/model_runs"))
  
  for(proj in listStacks) {
    
    cat(paste("---- ", Sys.time(), "Projection:", proj, "----"))
    clim_data <- raster("./inputs/envdata/future/processed/temp200_", proj, "_1deg.tif")
    
    setwd(paste0("./outputs/modelling/models/", clst))
    
    projection_runs <- BIOMOD_Projection(modeling.output = model_runs, 
                                         new.env = stack(clim_data), # Environmental data used for the projection
                                         proj.name = proj,
                                         selected.models = 'all', # Which models to project
                                         binary.meth = NULL, # no binary transformation, keep continuous probability
                                         filtered.meth = NULL, # no threshold to assign 0 to a probability
                                         build.clamping.mask = TRUE, # Create a clamping mask to show regions where predictions are out of the data used for calibration
                                         do.stack = F) # Save the projections (all calibration models and CV runs) in separate files
                                         # do.stack = T) # Save the projections in a single file
    
    # Create a specific folder
    if(!exists(paste0(clst, "/proj_", proj))) 
    {
      dir.create(paste0(clst, "/proj_", proj), recursive = T) 
    }
    
    save(projection_runs, file = paste0(clst, "/proj_", proj, "/projection_runs"))
    
    cat(paste("---- ", Sys.time(), "Projection:", proj, "finished ----\n"))
    
    # Clean memory
    rm(list = c("projection_runs"))    
    unlink(rasterOptions()$tmpdir, force = T)
    
    setwd(initial.wd)
  }
} 



for (clst in list_clust) {
  
  print(clst)
  bioregion <- ifelse(clst == "cluster1", "Subtropical", "Southern")
  setwd(paste0("./outputs/modelling/models/", clst, "/", clst))
  
  for(ssp in SSPs) {
    
    print(ssp)
    scenario <- ifelse(ssp == "ssp126", "SSP1-2.6", ifelse(ssp == "ssp245", "SSP2-4.5", "SSP5-8.5"))
    
    for(horizon in horizons) {
      
      print(horizon)
      model_stack <- stack()
      
      for(GCM in GCMs) {
        
        print(GCM)
        proj <- paste0(ssp, "_", horizon, "_", GCM)
        
        if (!proj %in% c("ssp245_2100_FGOALS-g3", "ssp126_2100_IPSL-CM6A-LR")) {  # Missing climatic data
          # Download projection data 
          cur_proj <- stack(paste0("proj_", proj, "/proj_", proj, "_", clst, ".grd"))  # continuous projection on baseline
          model_stack <- stack(model_stack, cur_proj)
        }
      }
      
      # Ensemble model as average of all calibration models and climatic models
      cur_em <- mean(model_stack)
      
      # Uncertainty (standard deviation)
      cur_sd <- calc(model_stack, sd) 
      
      # Penalised probability of presence
      cur_em_penalised <- cur_em - cur_sd 
      cur_em_penalised[cur_em_penalised < 0] <- 0
      
      # Save
      # writeRaster(cur_em, paste0("proj_", ssp, "_", horizon, "_", clst, "_em.grd"), overwrite = T)
      writeRaster(cur_sd, paste0("proj_", ssp, "_", horizon, "_", clst, "_sd.grd"), overwrite = T)
      writeRaster(cur_em_penalised, paste0("proj_", ssp, "_", horizon, "_", clst, "_em_pen.grd"), overwrite = T)
      
      
      # Maps
      palette <- viridis(6)
      
      proj_em <- tm_shape(cur_em/1000) +
        tm_raster(title = "Probability of \nbioregion occurrence", palette = palette) +
        tm_shape(land) +
        tm_borders(col = "grey") +
        tm_shape(graticules) +
        tm_lines(col = "grey", alpha = 0.5) +
        tm_layout(main.title = paste0("Projection ", bioregion, " bioregion ", horizon, " (", scenario, ")"), 
                  main.title.size = 1, asp = 1, legend.outside = TRUE)
      
      proj_em_pen <- tm_shape(cur_em_penalised/1000) +
        tm_raster(title = "Probability of \nbioregion occurrence", palette = palette) +
        tm_shape(land) +
        tm_borders(col = "grey") +
        tm_layout(main.title = paste0("Penalised projection ", bioregion, " bioregion ", horizon, " (", scenario, ")"), 
                  main.title.size = 1, asp = 2, legend.outside.size = 0.15)
       
      proj_sd <- tm_shape(cur_sd/1000) +
        tm_raster(title = "Probability of \nbioregion occurrence", palette = palette) +
        tm_shape(land) +
        tm_borders(col = "grey") +
        tm_layout(main.title = paste0("Standard deviation for the ", bioregion, " bioregion ", horizon, " (", scenario, ")"), 
                  main.title.size = 1, asp = 2, legend.outside.size = 0.15)
      
      # Save maps
      tmap_save(tm = proj_em,
                filename = paste0(initial.wd, "./outputs/plots/modelling/model_projections/proj_", ssp, "_", horizon, "_", clst, "_em.png"),
                width = 2000, height = 1200,
                outer.margins = c(0, 0.01, 0, 0))
      
      tmap_save(tm = proj_em_pen,
                filename = paste0(initial.wd, "./outputs/plots/modelling/model_projections/proj_", ssp, "_", horizon, "_", clst, "_em_pen.png"),                width = 2000, height = 1200,
                outer.margins = c(0, 0.01, 0, 0))

      tmap_save(tm = proj_sd,
                filename = paste0(initial.wd, "./outputs/plots/modelling/model_projections/proj_", ssp, "_", horizon, "_", clst, "_sd.png"),                width = 2000, height = 1200,
                outer.margins = c(0, 0.01, 0, 0))
    }
  }
  setwd(initial.wd)
}


## 2. Present-future difference -----------------------------------------     

list_clust <- c("cluster1", "cluster2")

for (clst in list_clust) {
  
  proj_present <- stack(paste0("./outputs/modelling/models/", clst, "/", clst, "/proj_baseline/proj_baseline_", clst, "_em.grd"))
  proj_future <- stack(paste0("./outputs/modelling/models/", clst, "/", clst, "/proj_ssp585_2100_", clst, "_em.grd"))
  
  palettediff <- colorRampPalette(colors = c("deepskyblue3", "lemonchiffon", "firebrick"))(12)
  
  proj_diff <- tm_shape((proj_future - proj_present)/1000) +
    tm_raster(title = "\u394Probabilities", palette = palettediff, 
              breaks = c(seq(-1, 0, length.out = 6), seq(0, 1, length.out = 6))) +
    tm_shape(land) +
    tm_borders(col = "grey") +
    tm_shape(graticules) +
    tm_lines(col = "grey", alpha = 0.5) +
    tm_layout(main.title =  "Difference (2100 - present) SSP5-8.5", main.title.size = 1.3, asp = 1, legend.outside = TRUE)
  
  tmap_save(tm = proj_diff,
            filename = paste0(initial.wd, "/outputs/plots/modelling/diff_ssp585_2100_", clst, "_em.png"),
            width = 2000, height = 1200,
            outer.margins = c(0, 0.01, 0, 0))
}


## 3. Present-future evolution ENVIRONMENT -----------------------------------------     

palette <- colorRampPalette(colors = c("dodgerblue3", "olivedrab3", "gold", "firebrick"))(6)

for (ssp in SSPs) {
  
  for(horizon in horizons) {
    T_stack <- stack()
  
    for(GCM in GCMs) {
      proj <- paste0(ssp, "_", horizon, "_", GCM)
      
      if (!proj %in% c("ssp245_2100_FGOALS-g3", "ssp126_2100_IPSL-CM6A-LR")) {  # Missing climatic data
        # Download projection data 
        proj_GCM <- readRDS(paste0("./outputs/", proj, ".RDS"))
        T_stack <- stack(T_stack, proj_GCM)
      }
    }
    
    # Average of all GCM values    
    T_mean <- mean(T_stack)
    
    map <- tm_shape(T_mean) +
      tm_raster(title = "Temperature at 200 m (°C)", style = "cont", midpoint = NA, 
                palette = palette, breaks = seq(-2, 20, 2), legend.is.portrait = FALSE) +
      tm_shape(land) +
      tm_borders(col = "grey") +
      tm_shape(graticules) +
      tm_lines(col = "grey", alpha = 0.5) +
      tm_layout(main.title = paste(ssp, "-", horizon), asp = 1, legend.outside = TRUE)
    
    tmap_save(tm = map,
              filename = paste0("./outputs/plots/environmental_variables/Temp200_", ssp, "_", horizon, "_review.png"),
              outer.margins = c(0, 0.01, 0, 0))
  }
}



## 4. Quantify probability DIFFERENCE ----------------------------------------- 

## Load data

models <- c("GLM", "GAM", "ANN", "MARS", "FDA")
cv_runs <- c("RUN1", "RUN2", "RUN3", "RUN4")
list_clust <- c("cluster1", "cluster2")

list_models <- expand.grid(cv_runs, models)
list_models <- apply(list_models, 1,
                     function(x) paste0(x[1], "_", x[2]))


# Significance threshold (for statistical test)
signif_thresh <- 0.05


## Plot

for (clst in list_clust) {
  
  print(paste("------- Currently working on", clst, "-------"))
  bioreg <- ifelse(clst == "cluster1", "Subtropical", "Southern")
  setwd(paste0("./outputs/modelling/models/", clst, "/", clst))
  
  
  ## A. PRESENT ---
  
  ## STAT TEST - Create empty dataframe to store all the present counts
  df_counts_present <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(df_counts_present) <- c("breaks", "time", "counts")
  
  ## A1. Mean  ---
  ## DIFF - Load present data (mean values = ensemble model) 
  # to calculate the difference future-present
  proj_present_mean <- raster(paste0("./proj_baseline/proj_baseline_", clst, "_em.grd"))/1000
  hist_present_mean <- hist(proj_present_mean, main = "Present", 
                            breaks = 5, xlim = range(0,1), 
                            xlab = "Probability of region occurrence", 
                            plot = F) # Create histogram for the present and get the data
  df_present_mean <- data.frame(hist_present_mean$breaks[-1], hist_present_mean$counts)
  colnames(df_present_mean) <- c("breaks", "counts_present")
  
  ## A2. All values ---
  
  ## SD & STAT TEST
  # Load present data for all calibration models and CV runs projections 
  # to calculate the standard deviation and statistical test
  for(model in list_models) {  
    
    ## A2a. Load data ---
    
    # Download projection data
    proj_present <- stack(paste0("./proj_baseline/individual_projections/proj_baseline_", clst, "_AllData_", model,".grd"))/1000 # continuous projection on baseline
    
    # Create histogram for the present and get the data
    hist_present <- hist(proj_present, main = paste("Present -", model), 
                         breaks = 5, xlim = range(0,1), 
                         xlab = "Probability of region occurrence",
                         plot = F)
    df_present_model <- data.frame(hist_present$breaks[-1], "present", hist_present$counts)
    colnames(df_present_model) <- c("breaks", "time", "counts")
    df_present_model$model_cvrun <- model
    
    df_counts_present <- rbind(df_counts_present, df_present_model)
  }
  
  ## A2b. Standard deviation ---
  
  # SD - Calculate standard deviation fo all present values
  df_present_ed <- aggregate(counts ~ breaks, data = df_counts_present, FUN = sd)
  colnames(df_present_ed) <- c("breaks", "sd_present")
  
  
  ## B. FUTURE ---
  
  for(horizon in horizons) {
    
    print(paste("-------------- Future horizon :", horizon))
    
    # Create empty dataframe to store all the future counts
    ## DIFF & SD
    df_sd_future <- data.frame()
    df_mean_evol <- data.frame()
    # Create empty dataframe to store the results of the statistical tests
    test_results <- data.frame(matrix(ncol = 4, nrow = 0))
    colnames(test_results) <- c("breaks", "ssp", "test", "p_value")
    
    
    for(ssp in SSPs) {
      
      print(paste("------- Scenario :", ssp))
      ## STAT TEST
      df_counts_future <- data.frame(matrix(ncol = 5, nrow = 0))
      colnames(df_counts_future) <- c("breaks", "time", "counts", 
                                      "model_cvrun", "ssp")
      
      ## B1. Mean ---
      
      ## DIFF - Load future data (mean values = ensemble model) 
      # to calculate the difference future-present
      proj_future_mean <- raster(paste0(initial.wd, "/outputs/modelling/models/", clst, "/", clst, "/proj_", ssp, "_", horizon, "_", clst, "_em.grd"))/1000
      hist_future_mean <- hist(proj_future_mean, main = paste0(horizon, " (", ssp, ")"), 
                               breaks = 5, xlim = range(0,1), 
                               xlab = "Probability of region occurrence",
                               plot = F) # Create histogram for the future and get the data
      df_future_mean <- data.frame(hist_future_mean$breaks[-1], hist_future_mean$counts)
      colnames(df_future_mean) <- c("breaks", "counts_future")
      
      
      ## B2. Difference F-P ---
      
      # Calculate difference future-present
      df_future_mean <- inner_join(df_present_mean, df_future_mean, by = "breaks")
      df_future_mean$diff <- df_future_mean$counts_future - df_future_mean$counts_present
      df_future_mean$ssp <- ssp
      
      df_mean_evol <- rbind(df_mean_evol, dplyr::select(df_future_mean, breaks, diff, ssp))
      
      
      ## B3. All values ---
      
      for(GCM in GCMs) {
        
        print(GCM)
        proj <- paste0(ssp, "_", horizon, "_", GCM)
        
        if (!proj %in% c("ssp245_2100_FGOALS-g3", "ssp126_2100_IPSL-CM6A-LR")) {  # Missing climatic data
          
          # SD & STAT TEST - Load future data for all calibration models and CV run projections
          # to calculate the standard deviation and statistical test
          for(model in list_models) {  
            
            ## B3a. Load data ---
            
            # Download projection data 
            proj_future <- raster(paste0("proj_", proj, "/individual_projections/proj_", proj, "_", clst, "_AllData_", model, ".grd"))/1000  # continuous projection on baseline
            
            # Create histogram for the future and get the data
            hist_future <- hist(proj_future, main = paste(horizon, "(", ssp, ") -", model), 
                                breaks = 5, xlim = range(0,1), 
                                xlab = "Probability of region occurrence",
                                plot = F)
            df_future_model <- data.frame(hist_future$breaks[-1], "future", hist_future$counts)
            colnames(df_future_model) <- c("breaks", "time", "counts")
            df_future_model$model_cvrun <- model
            df_future_model$ssp <- ssp
            
            ## SD
            df_sd_future <- rbind(df_sd_future, df_future_model)
            ## STAT TEST
            df_counts_future <- rbind(df_counts_future, 
                                      df_future_model)
          }
        }
      }
      
      
      ## B3b. Statistical test ---
      
      ## STAT TEST
      # Mean of all GCMs
      df_counts_future <- aggregate(df_counts_future$counts,
                                    by = list(df_counts_future$breaks,
                                              df_counts_future$time,
                                              df_counts_future$model_cvrun,
                                              df_counts_future$ssp), 
                                    FUN = mean)
      colnames(df_counts_future) <- c("breaks", "time", "model_cvrun", "ssp", "counts")
      
      # Present and future counts in the same dataframe
      df_counts_evol <- rbind(dplyr::select(df_counts_present, -model_cvrun), 
                              dplyr::select(df_counts_future, -model_cvrun, -ssp))
      
      # Statistical test 
      for (proba in unique(df_counts_evol$breaks)) {
        df_counts_evol_proba <- df_counts_evol %>%
          filter(breaks == proba)
        
        # Verify statistical assumptions - normality and outliers
        df_counts_diff <- left_join(filter(df_counts_present, breaks == proba),
                                    filter(dplyr::select(df_counts_future, -ssp), breaks == proba),
                                    by = c("breaks", "model_cvrun"))
        df_counts_diff$diff <- df_counts_diff$counts.y - df_counts_diff$counts.x
        
        norm_test <- shapiro_test(df_counts_diff, diff)$p
        outlier_test <- identify_outliers(df_counts_diff, diff)$is.outlier
        
        # Paired t-test
        if (norm_test > 0.05 & !any(outlier_test)) {
          
          result <- t_test(counts ~ time, 
                           data = df_counts_evol_proba, 
                           paired = T) %>%
            add_significance()
          test_results <- rbind(test_results, 
                                data.frame(breaks = proba, 
                                           ssp = ssp,
                                           test = "t-test",
                                           p_value = result$p,
                                           signif = ifelse(result$p < signif_thresh,
                                                           yes = "*",
                                                           no = "")))
          # Sign test
        } else {
          
          result <- df_counts_evol_proba %>% 
            sign_test(counts ~ time) %>%
            add_significance()
          
          test_results <- rbind(test_results, 
                                data.frame(breaks = proba, 
                                           ssp = ssp,
                                           test = "sign test",
                                           p_value = result$p,
                                           signif = ifelse(result$p < signif_thresh,
                                                           yes = "*",
                                                           no = "")))
        }
        
      }
      
      
      ## DIFF & STAT TEST
      # Add p-value to dataframe with difference values
      df_evol <- inner_join(df_mean_evol, test_results)
    }
    
    
    ## B3c. Standard deviation ---
    
    ## SD
    # Calculate standard deviation for each scenario
    df_ssp126 <- subset(df_sd_future, ssp == "ssp126")
    df_sd_ssp126 <- aggregate(df_ssp126$counts, by = list(df_ssp126$breaks), FUN = sd) # Standard deviation
    colnames(df_sd_ssp126) <- c("breaks", "sd_future")
    df_sd_ssp126$ssp <- "ssp126"
    
    df_ssp245 <- subset(df_sd_future, ssp == "ssp245")
    df_sd_ssp245 <- aggregate(df_ssp245$counts, by = list(df_ssp245$breaks), FUN = sd)
    colnames(df_sd_ssp245) <- c("breaks", "sd_future")
    df_sd_ssp245$ssp <- "ssp245"
    
    df_ssp585 <- subset(df_sd_future, ssp == "ssp585")
    df_sd_ssp585 <- aggregate(df_ssp585$counts, by = list(df_ssp585$breaks), FUN = sd)
    colnames(df_sd_ssp585) <- c("breaks", "sd_future")
    df_sd_ssp585$ssp <- "ssp585"
    
    # Group the sd in one data frame
    df_future_sd <- rbind(df_sd_ssp126, df_sd_ssp245, df_sd_ssp585)
    
    
    ## C. Gather the values ---
    
    # Group the sd values (uncertainty between future model predictions) and the mean difference values
    df_evol <- inner_join(df_evol, df_future_sd, by = c("breaks", "ssp"))
    # Add the sd values for the present
    df_evol <- inner_join(df_evol, df_present_ed, by = "breaks")
    
    df_evol$center <- 0   # A trick to center the present variability bars around 0
    df_evol$breaks <- as.character(df_evol$breaks)
    
    # Calculate total area (in million km2)
    df_evol$diff <- df_evol$diff * 6242.4 / 10^6
    df_evol$sd_future <- df_evol$sd_future * 6242.4 / 10^6
    df_evol$sd_present  <- df_evol$sd_present  * 6242.4 / 10^6
    
    # Rename probability classes
    df_evol$breaks <- ifelse(df_evol$breaks == "0.2", "very low", ifelse(df_evol$breaks == "0.4", "low", ifelse(df_evol$breaks == "0.6", "medium", ifelse(df_evol$breaks == "0.8", "high", "very high"))))
    df_evol$breaks <- factor(df_evol$breaks, levels = c("very low", "low", "medium", "high", "very high"))
    
    
    ## D. Plot ---
    
    color_bioregion <- ifelse(bioreg == "Southern", yes = "#56B4E9", no = "#D55E00") # color of the bars according to the bioregion
    panel_tag <- ifelse(i == 1, "(a)", 
                        ifelse(i == 2, "(b)",
                               ifelse(i == 3, "(c)", "(d)")))


    plot_evol <- ggplot(df_evol) +
      geom_pointrange(aes(breaks, center,
                          ymin = center - sd_present/2, ymax = center + sd_present/2),
                      fill = color_bioregion, color = color_bioregion,
                      shape = 32, alpha = 0.2,  size = 18, linewidth = 15) +
      geom_pointrange(aes(breaks, diff, col = ssp, 
                          ymin = diff - sd_future/2, ymax = diff + sd_future/2),
                      position = position_dodge(width = 0.3), size = 0.7) +
      # Add * if the difference is significant (t-test)
      geom_text(aes(x = breaks, y = diff + 2.5, label = signif, group = ssp),
                position = position_dodge(width = 0.3), size = 6) +
      geom_hline(yintercept = 0, linetype = "dashed", size = 0.3) +
      ylim(-18.3, 18.3) +
      xlab("Probability of region occurrence") +
      ylab(paste0("\u394 Area (million km\u00B2)")) +
      scale_color_manual(name = "Scenario", 
                         labels = c("SSP1-2.6", "SSP2-4.5", "SSP5-8.5"),
                         values = c("#66a61e", "#7570b3", "#e7298a")) +
      ggtitle(paste0(bioreg, " region by ", horizon)) +
      theme(
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        plot.title = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 15), 
        legend.text = element_text(size = 14),
        legend.key = element_rect(fill = "transparent", colour = "transparent"), # Remove frame on legend
        plot.tag = element_text(size = 17, face = "bold"),
        panel.background = element_rect(fill = "white", colour = "grey"),
        panel.grid.major.y = element_line(colour = "grey", size = 0.2),
        panel.grid.major.x = element_blank())
    
    
    ggsave(paste0(initial.wd, "/outputs/plots/modelling/diffPlot_", horizon, "_", clst, ".png"), 
           plot_evol, height = 1200, width = 2000, units = "px")

  }
  setwd(initial.wd)
}





## 5a. Statistical test present - future difference (t-test) -----------------------------------------
## V1. Non paired test

models <- c("GLM", "GAM", "ANN", "MARS", "FDA")
cv_runs <- c("RUN1", "RUN2", "RUN3", "RUN4")
list_clust <- c("cluster1", "cluster2")

list_models <- expand.grid(cv_runs, models)
list_models <- apply(list_models, 1,
                     function(x) paste0(x[1], "_", x[2]))

# Create empty dataframe to store the results of th statistical tests
test_results <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(test_results) <- c("bioregion", "temp_horizon", "scenario", "class_proba", "p_value")

for (clst in list_clust) {
  
  print(paste("------- Currently working on", clst, "-------"))
  bioreg <- ifelse(clst == "cluster1", "Subtropical", "Southern")
  setwd(paste0("./outputs/modelling/models/", clst, "/", clst))
  
  # Create empty dataframe to store all the counts
  df_counts_present <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(df_counts_present) <- c("breaks", "time", "counts")
  
  
  for(model in list_models) {  
    
    # Download projection data
    proj_present <- stack(paste0("./proj_baseline/individual_projections/proj_baseline_", clst, "_AllData_", model,".grd"))/1000 # continuous projection on baseline
    
    # Create histogram for the present and get the data
    hist_present <- hist(proj_present, main = paste("Present -", model), 
                         breaks = 5, xlim = range(0,1), 
                         xlab = "Probability of bioregion occurrence",
                         plot = F)
    df_present_model <- data.frame(hist_present$breaks[-1], "present", hist_present$counts)
    colnames(df_present_model) <- c("breaks", "time", "counts")
    
    df_counts_present <- rbind(df_counts_present, df_present_model)
  }
  
  for(horizon in horizons) {
    
    print(paste("-------------- Future horizon :", horizon))
    df_sd_future <- data.frame()
    df_evol <- data.frame()
    
    for(ssp in SSPs) {
      
      print(paste("------- Scenario :", ssp))
      df_counts_evol <- df_counts_present
      
      for(GCM in GCMs) {
        
        print(GCM)
        proj <- paste0(ssp, "_", horizon, "_", GCM)
        
        if (!proj %in% c("ssp245_2100_FGOALS-g3", "ssp126_2100_IPSL-CM6A-LR")) {  # Missing climatic data
          
          # Future data for all calibration models and CV runs projections (to calculate the standard deviation)
          for(model in list_models) {  
            
            # Download projection data 
            proj_future <- raster(paste0("proj_", proj, "/individual_projections/proj_", proj, "_", clst, "_AllData_", model, ".grd"))/1000  # continuous projection on baseline
            
            # Create histogram for the future and get the data
            hist_future <- hist(proj_future, main = paste(horizon, "(", ssp, ") -", model), 
                                breaks = 5, xlim = range(0,1), 
                                xlab = "Probability of bioregion occurrence",
                                plot = F)
            df_future_model <- data.frame(hist_future$breaks[-1], "future", hist_future$counts)
            colnames(df_future_model) <- c("breaks", "time", "counts")
            
            df_counts_evol <- rbind(df_counts_evol, df_future_model)
          }
        }
      }
      
      # Statistical test 
      for (proba in unique(df_counts_evol$breaks)) {
        df_counts_evol_proba <- df_counts_evol %>%
          filter(breaks == proba)
        result <- t_test(counts ~ time, 
                         data = df_counts_evol_proba, 
                         paired = F) %>%
          add_significance()
        test_results <- rbind(test_results, 
                              data.frame(bioregion = bioreg, 
                                         temp_horizon = horizon, 
                                         scenario = ssp, 
                                         class_proba = proba, 
                                         p_value = result$p))
      }
    }
  }
  setwd(initial.wd)
}

# Significance threshold 
signif_thresh <- 0.5
test_results$signif <- ifelse(test_results$p_value < signif_thresh,
                              yes = "Significant",
                              no = "Not significant")


## 5b. Statistical test present - future difference (t-test) -----------------------------------------
## V2. Paired test

models <- c("GLM", "GAM", "ANN", "MARS", "FDA")
cv_runs <- c("RUN1", "RUN2", "RUN3", "RUN4")
list_clust <- c("cluster1", "cluster2")

list_models <- expand.grid(cv_runs, models)
list_models <- apply(list_models, 1,
                     function(x) paste0(x[1], "_", x[2]))

# Create empty dataframe to store the results of th statistical tests
test_results <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(test_results) <- c("bioregion", "temp_horizon", "scenario", "class_proba", "p_value")

for (clst in list_clust) {
  
  print(paste("------- Currently working on", clst, "-------"))
  bioreg <- ifelse(clst == "cluster1", "Subtropical", "Southern")
  setwd(paste0("./outputs/modelling/models/", clst, "/", clst))
  
  # Create empty dataframe to store all the present counts
  df_counts_present <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(df_counts_present) <- c("breaks", "time", "counts")
  
  
  for(model in list_models) {  
    
    # Download projection data
    proj_present <- stack(paste0("./proj_baseline/individual_projections/proj_baseline_", clst, "_AllData_", model,".grd"))/1000 # continuous projection on baseline
    
    # Create histogram for the present and get the data
    hist_present <- hist(proj_present, main = paste("Present -", model), 
                         breaks = 5, xlim = range(0,1), 
                         xlab = "Probability of bioregion occurrence",
                         plot = F)
    df_present_model <- data.frame(hist_present$breaks[-1], "present", hist_present$counts)
    colnames(df_present_model) <- c("breaks", "time", "counts")
    # df_present_model$model_cvrun <- model
    
    df_counts_present <- rbind(df_counts_present, df_present_model)
  }
  
  for(horizon in horizons) {
    
    print(paste("-------------- Future horizon :", horizon))
    # df_sd_future <- data.frame()
    # df_evol <- data.frame()
    
    for(ssp in SSPs) {
      
      print(paste("------- Scenario :", ssp))
      # Create empty dataframe to store all the future counts
      df_counts_future <- data.frame(matrix(ncol = 5, nrow = 0))
      colnames(df_counts_future) <- c("breaks", "time", "counts", 
                                      "model_cvrun", "scenario")
        
      for(GCM in GCMs) {
        
        print(GCM)
        proj <- paste0(ssp, "_", horizon, "_", GCM)
        
        if (!proj %in% c("ssp245_2100_FGOALS-g3", "ssp126_2100_IPSL-CM6A-LR")) {  # Missing climatic data
          
          # Future data for all calibration models and CV runs projections (to calculate the standard deviation)
          for(model in list_models) {  
            
            # Download projection data 
            proj_future <- raster(paste0("proj_", proj, "/individual_projections/proj_", proj, "_", clst, "_AllData_", model, ".grd"))/1000  # continuous projection on baseline
            
            # Create histogram for the future and get the data
            hist_future <- hist(proj_future, main = paste(horizon, "(", ssp, ") -", model), 
                                breaks = 5, xlim = range(0,1), 
                                xlab = "Probability of bioregion occurrence",
                                plot = F)
            df_future_model <- data.frame(hist_future$breaks[-1], "future", hist_future$counts)
            colnames(df_future_model) <- c("breaks", "time", "counts")
            df_future_model$model_cvrun <- model
            df_future_model$scenario <- ssp
            
            df_counts_future <- rbind(df_counts_future, df_future_model)
          }
        }
      }
      
      # Mean of all GCMs
      df_counts_future <- aggregate(df_counts_future$counts,
                                    by = list(df_counts_future$breaks,
                                              df_counts_future$time,
                                              df_counts_future$model_cvrun,
                                              df_counts_future$scenario), 
                                    FUN = mean)
      colnames(df_counts_future) <- c("breaks", "time", "model_cvrun", "scenario", "counts")
      
      # Present and future counts in the same dataframe
      df_counts_evol <- rbind(df_counts_present, 
                              select(df_counts_future, -model_cvrun, -scenario))
      
      # Statistical test 
      for (proba in unique(df_counts_evol$breaks)) {
        df_counts_evol_proba <- df_counts_evol %>%
          filter(breaks == proba)
        result <- t_test(counts ~ time, 
                         data = df_counts_evol_proba, 
                         paired = T) %>%
          add_significance()
        test_results <- rbind(test_results, 
                              data.frame(bioregion = bioreg, 
                                         temp_horizon = horizon, 
                                         scenario = ssp, 
                                         class_proba = proba, 
                                         p_value = result$p))
      }
    }
  }
  setwd(initial.wd)
}

# Significance threshold 
signif_thresh <- 0.05
test_results$signif <- ifelse(test_results$p_value < signif_thresh,
                              yes = "*",
                              no = "")



## 6. Compare area of bioregions : present-future ----------------

## All study area
present_subtrop <- stack("./outputs/modelling/models/cluster1/cluster1/proj_baseline/proj_baseline_cluster1_em.grd")/1000
present_southern <- stack("./outputs/modelling/models/cluster2/cluster2/proj_baseline/proj_baseline_cluster2_em.grd")/1000

future_subtrop <- stack("./outputs/modelling/models/cluster1/cluster1/proj_ssp585_2100_cluster1_em.grd")/1000
future_southern <- stack("./outputs/modelling/models/cluster2/cluster2/proj_ssp585_2100_cluster2_em.grd")/1000


## Compare area > threshold between present and future

threshold <- 0.8

# Present - subtropical
filter_present_subtrop <- present_subtrop < threshold
present_subtrop <- mask(present_subtrop, filter_present_subtrop, maskvalue = 1) # remove pixels with probability < threshold
plot(present_subtrop, main = paste0("Present subtropical bioregion (p > ", threshold, ")"))
nbcell_present_subtrop <- cellStats(present_subtrop, function(i, ...) sum(!is.na(i)))

# Future - subtropical
filter_future_subtrop <- future_subtrop < threshold
future_subtrop <- mask(future_subtrop, filter_future_subtrop, maskvalue = 1)
plot(future_subtrop, main = paste0("Future subtropical bioregion (p > ", threshold, ")"))
nbcell_future_subtrop <- cellStats(future_subtrop, function(i, ...) sum(!is.na(i)))

diff_area_subtrop <- nbcell_future_subtrop * 6242.4 / 10^6 - nbcell_present_subtrop * 6242.4 / 10^6


# Present - southern
filter_present_southern <- present_southern < threshold
present_southern <- mask(present_southern, filter_present_southern, maskvalue = 1) # remove pixels with probability < threshold
plot(present_southern, main = "Present southern bioregion (p > 0.5)")
nbcell_present_southern <- cellStats(present_southern, function(i, ...) sum(!is.na(i)))

# Future - southern
filter_future_southern<- future_southern < threshold
future_southern <- mask(future_southern, filter_future_southern, maskvalue = 1)
plot(future_southern, main = "Future southern bioregion (p > 0.5)")
nbcell_future_southern <- cellStats(future_southern, function(i, ...) sum(!is.na(i)))

diff_area_southern <- nbcell_future_southern * 6242.4 / 10^6 - nbcell_present_southern * 6242.4 / 10^6


cat("Subtropical region area difference:", diff_area_subtrop, "km2\n")
cat("Southern region area difference:", diff_area_southern, "km2\n")



## 7. Taylor diagram - compare historical T200 datasets (observed baseline/CMIP6 projections) ----------------

library(raster)
library(archive)
library(plotrix)


## Baseline dataset (1955-2017)
baseline_30 <- brick("./inputs/envdata/baseline/processed/env_baseline_1deg.tif")
baseline_val <- getValues(baseline_30$temperature200)

## CMIP6 historical dataset (1950-2014)
list_gcm <- c("AWI-CM-1-1-MR", "FGOALS-g3", "CanESM5", "CMCC-ESM2", "ACCESS-CM2",
              "FIO-ESM-2-0",  "IPSL-CM6A-LR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "CIESM")
cmip_val <- list()

for(GCM in list_gcm) {
  
  print(paste("Processing", GCM))
  file <- paste0("./inputs/envdata/future/raw/", GCM, "_hist_final")
  
  if (!file.exists(paste0(file, ".nc"))) {
    if (file.exists(paste0(file, ".7z"))) {
      # Unzip file
      archive_extract(paste0(file, ".7z"), dir = "./inputs/envdata/future/raw/")
    } else {
      # Deal with issue of missing files - recover hist data from final data
      print(paste("There is no such file for", GCM))
      print("Recovering historical data from final projection files...")
      
      final <- readRDS(paste0("./inputs/envdata/future/processed/ssp126_2060_", GCM, ".RDS"))
      absfuture <- brick(paste0("./inputs/envdata/future/raw/", GCM, "_ssp126_final.nc"), varname = "thetao", level = 10)
      season <- c("10", "11", "12", "01", "02", "03") 
      absfuture_season <- absfuture[[which(format(as.Date(getZ(absfuture), format = "X%Y.%m.%d"), format = "%m") %in% season)]]
      absfuture_year <- absfuture_season[[which(format(as.Date(names(absfuture_season), format = "X%Y.%m.%d")) >= as.Date("2041-01-01") & format(as.Date(names(absfuture_season), format = "X%Y.%m.%d")) <= as.Date("2060-12-31"))]]
      absfuture_mean <- mean(absfuture_year, na.rm = TRUE)
      future_proj <- projectRaster(rotate(absfuture_mean), baseline_30, method = "bilinear")
      
      hist_cmip6 <- baseline_30$temperature200 + future_proj - final
      cmip_val[[GCM]] <- getValues(hist_cmip6)
    }
  } else {
    
    # Temperature at 200 m (level 10)
    hist <- brick(paste0(file, ".nc"), varname = "thetao", level = 10)
    # Seasons = October-March
    season <- c("10", "11", "12", "01", "02", "03") 
    hist_season <- hist[[which(format(as.Date(getZ(hist), format = "X%Y.%m.%d"), format = "%m") %in% season)]]
    # Time period = 1950-2014
    hist_year <- hist_season[[which(format(as.Date(names(hist_season), format = "X%Y.%m.%d")) >= as.Date("1950-01-01") & 
                                      format(as.Date(names(hist_season), format = "X%Y.%m.%d")) <= as.Date("2014-12-31"))]]
    # Values = average on the time period
    hist_mean <- mean(hist_year, na.rm = TRUE)
    # Set the same extent, origin, projection as the baseline
    hist_cmip6 <- projectRaster(rotate(hist_mean), baseline_30, method = "bilinear")
    
    cmip_val[[GCM]] <- getValues(hist_cmip6)
    
  }
}  

# Save Taylor diagram
png("./outputs/plots/environmental_variables/taylordiagram_histdata.png")
taylor.diagram(baseline_val, cmip_val[["AWI-CM-1-1-MR"]], col = "#D55E00", ref.sd = T, normalize = T)
taylor.diagram(baseline_val, cmip_val[["FGOALS-g3"]], col = "#CC79A7", normalize = T, add = TRUE)
taylor.diagram(baseline_val, cmip_val[["CanESM5"]], col = "#0072B2", normalize = T, add = TRUE)
taylor.diagram(baseline_val, cmip_val[["CMCC-ESM2"]], col = "#F0E442", normalize = T, add = TRUE)
taylor.diagram(baseline_val, cmip_val[["ACCESS-CM2"]], col = "#009E73", normalize = T, add = TRUE)
taylor.diagram(baseline_val, cmip_val[["FIO-ESM-2-0"]], col = "#D55E00", pch = 17, normalize = T, add = TRUE)
taylor.diagram(baseline_val, cmip_val[["IPSL-CM6A-LR"]], col = "#CC79A7", pch = 17, normalize = T, add = TRUE)
taylor.diagram(baseline_val, cmip_val[["MPI-ESM1-2-LR"]], col = "#0072B2", pch = 17, normalize = T, add = TRUE)
taylor.diagram(baseline_val, cmip_val[["MRI-ESM2-0"]], col = "#F0E442", pch = 17, normalize = T, add = TRUE)
taylor.diagram(baseline_val, cmip_val[["CIESM"]], col = "#009E73", pch = 17, normalize = T, add = TRUE)
legend(1.55, 2,legend = list_gcm[1:5], pch = 16, col = c("#D55E00", "#CC79A7", "#0072B2", "#F0E442", "#009E73"))
legend(1.55, 1.55,legend = list_gcm[6:10], pch = 17, col = c("#D55E00", "#CC79A7", "#0072B2", "#F0E442", "#009E73"))
dev.off()
