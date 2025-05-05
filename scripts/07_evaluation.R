###------------------------------------------###
###            MODEL EVALUATION              ###
###                                          ###
###  Jaccard index, AUC & variable response  ###
###------------------------------------------###


library(biomod2)
library(ggplot2)
library(dplyr)
library(cowplot)

clust_colors <- readRDS("./outputs/bioregionalisation/clust_colors.RDS")

list_clust <- c("cluster1", "cluster2", "cluster3", "cluster4")

initial.wd <- getwd()


##  1a. Variable response -----------------------------------------

for (clst in list_clust){
  
  setwd(paste0("./outputs/modelling/models/", clst))
  
  # Load formated data and calibrated models
  load("run_data")
  load("model_runs")
  
  # Variables used for calibration
  cur_vars <- model_runs@expl.var.names
  
  # Names of calibrated models
  models_to_plot <- BIOMOD_LoadModels(model_runs)
  
  resp <- response.plot2(models = models_to_plot,
                         Data = get_formal_data(model_runs, "expl.var"),   # environmental variables for each p/A point
                         fixed.var.metric = "median",  # method for fixing the other variables
                         show.variables = cur_vars, 
                         data_species = get_formal_data(model_runs, 'resp.var'))   # P/A data
  colnames(resp) <- c("Index", "Variable", "Var.value", "Model", "Response")
  
  p1 <- ggplot(resp, aes(x = Var.value, y = Response))+ 
    geom_line(alpha = 0.2, aes(group = Model)) +     # high transparency for individual models
    stat_smooth() +    # average curve
    facet_wrap(~Variable, scales = "free_x") + 
    theme_bw() + 
    ylim(0, 1.1) +   # To avoid flattened curves
    xlab("Variable value") +
    ggtitle(paste0("Average model response to environmental variables - ", clst))
  
  png(paste0(initial.wd, "/outputs/plots/model_evaluation/response_plot_average_", clst, ".png"), width = 2000, height = 1800, res = 300)
  print(p1)
  dev.off()

  for (var in cur_vars) {
    
    p2 <- ggplot(resp, aes(x = Var.value, y = Response))+ 
      geom_line() +    
      facet_wrap(~Model, scales = "free_x", ncol = length(models_to_plot)/4) + 
      theme_bw() + 
      ylim(0, 1.1) + 
      xlab("Variable value") +
      ggtitle(paste0("Model response to ", var, " - ", clst))
    
    png(paste0(initial.wd, "/outputs/plots/model_evaluation/response_plot_",var, "_", clst, ".png"), width = 2400, height = 1600, res = 200)
    print(p2)
    dev.off()
  }

  setwd(initial.wd)
  
}


## Same plot

resp_list <- list()

for (i in 1:2){
  
  setwd(paste0("./outputs/modelling/models/cluster", i))
  
  # Load formated data and calibrated models
  load("run_data")
  load("model_runs")
  
  # Variables used for calibration
  cur_vars <- model_runs@expl.var.names
  
  # Names of calibrated models
  models_to_plot <- BIOMOD_LoadModels(model_runs)
  
  resp <- response.plot2(models = models_to_plot,
                         Data = get_formal_data(model_runs, "expl.var"),   # environmental variables for each p/A point
                         fixed.var.metric = "median",  # method for fixing the other variables
                         show.variables = cur_vars, 
                         data_species = get_formal_data(model_runs, 'resp.var'))   # P/A data
  colnames(resp) <- c("Index", "Variable", "Var.value", "Model", "Response")
  
  resp_list[[i]] <- resp
  
  setwd(initial.wd)
}


# Get the fitted curves
ggp1 <- ggplot(data = resp_list[[1]], aes(x = Var.value, y = Response)) + 
  geom_line() +
  stat_smooth()
ggp_data1 <- ggplot_build(ggp1)$data[[2]]
ggp_data1$Bioregion <- "Subtropical"

ggp2 <- ggplot(data = resp_list[[2]], aes(x = Var.value, y = Response)) + 
  geom_line() +
  stat_smooth()
ggp_data2 <- ggplot_build(ggp2)$data[[2]]
ggp_data2$Bioregion <- "Southern"


# 2 clusters
ggp_data <- rbind(ggp_data1, ggp_data2)

ggp_data <- rename(ggp_data, "temperature200" = "x", "Response" = y)



## Colors
mycto_clusters <- readRDS("./outputs/bioregionalisation/mycto_df.RDS")
# Get matching colors & clusters
clust_colors <- dplyr::select(mycto_clusters, lvl1, color)
clust_colors <- filter(clust_colors, clust_colors$lvl1 %in% c("1", "2"))
clust_colors$Bioregion <- ifelse(clust_colors$lvl1 == "1", "Subtropical", "Southern")

my_colors <- clust_colors$color
names(my_colors) <- clust_colors$Bioregion


png("./outputs/plots/model_evaluation/response_plot_temperature200_2clusters.png", width = 2000, height = 1600, res = 300)
ggplot(ggp_data, aes(temperature200, Response, col = Bioregion)) +   
  geom_line() +
  scale_color_manual(values = my_colors) +
  theme(panel.background = element_rect(fill = "white", colour = "grey"),
        panel.grid.major = element_line(colour = "grey", size = 0.2)) +
  ggtitle("Model response of clusters to temperature200")
dev.off()


##  1b. Variable realisation -----------------------------------------

## Sites
clustnet <- readRDS("./outputs/bioregionalisation/clust_df.RDS")
# Only plot clusters 1, 2
clustnet <- filter(clustnet, clustnet$cluster %in% c("1", "2"))
clustnet$Bioregion <- ifelse(clustnet$cluster == "1", "Subtropical", "Southern")


## Species

# Match species and environmental data
fulldb <- readRDS("./inputs/occdata/processed/fulldb_after_checking.RDS")
class(fulldb) <- "data.frame" 
fullbaseline <- readRDS("./inputs/envdata/baseline/processed/fullbaseline.RDS")
env_data <- raster::extract(fullbaseline, fulldb[, c("decimalLongitude", "decimalLatitude")])
sp_env <- cbind(fulldb[, c("scientificName", "decimalLongitude", "decimalLatitude")], env_data)
sp_env <- na.omit(sp_env)  # Remove occurrences with NA (generally on land)
sp_env <- sp_env[sp_env$decimalLatitude < -30, ]  # Only in the study zone below 30°S



# Add cluster information
mycto_clusters <- readRDS("./outputs/bioregionalisation/mycto_clusters.RDS")
sp_clusters <- dplyr::select(mycto_clusters, Name, lvl1, color)
colnames(sp_clusters) <- c("scientificName", "cluster", "color")
sp_clusters$scientificName <- gsub(".", " ", sp_clusters$scientificName, fixed = TRUE) # replace points by spaces
sp_env_cluster <- left_join(sp_env, sp_clusters)

# Only plot clusters 1, 2
sp_env_cluster <- filter(sp_env_cluster, sp_env_cluster$cluster %in% c("1", "2"))
sp_env_cluster$Bioregion <- ifelse(sp_env_cluster$cluster == "1", "Subtropical", "Southern")

# Plots
p1 <- ggplot(ggp_data, aes(temperature200, Response, col = Bioregion)) +   
  geom_line(size = 2) +
  scale_color_manual(values = my_colors) +
  xlab("") +
  ylab("Response") +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 22),
        legend.title = element_text(size = 18), 
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white", colour = "grey"),
        panel.grid.major = element_line(colour = "grey", size = 0.2)) +
  ggtitle("Predicted functional response of bioregion models")

p2 <- ggplot() +
  geom_violin(data = clustnet, 
               aes(temperature200, Bioregion, fill = Bioregion)) +
  scale_fill_manual(values = my_colors) +
  xlab("") +
  ylab("") +  
  theme(legend.position = "none",
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 22),
        panel.background = element_rect(fill = "white", colour = "grey"),
        panel.grid.major = element_line(colour = "grey", size = 0.2)) +
  ggtitle("Distribution of sites")

p3 <- ggplot() +
  geom_violin(data = sp_env_cluster, 
              aes(temperature200, Bioregion, fill = Bioregion)) +
  scale_fill_manual(values = my_colors) +
  xlab("Temperature at 200 m") +
  ylab("") +  
  theme(legend.position = "none",
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 22),
        panel.background = element_rect(fill = "white", colour = "grey"),
        panel.grid.major = element_line(colour = "grey", size = 0.2)) +
  ggtitle("Distribution of species")


png("./outputs/plots/model_evaluation/response_distribution_temp200.png", width = 1600, height = 1400, res = 150)
plot_grid(p1 + xlim(-3, 18), 
          p2 + xlim(-3, 18),
          p3 + xlim(-3, 18),
          align = 'v', rel_heights = c(3, 2, 2), nrow = 3, axis = "rl")
dev.off()



##  2. Jaccard index -----------------------------------------

load("./outputs/modelling/models/cluster1/model_runs") # to get the number of runs of cross-validation (values stored on disk following the link)

# Create an array to store Jaccard evaluations of all models
allclust_jaccard <- array(dim = c(length(list_clust),   # nb species
                               length(dimnames(model_runs@models.evaluation@val)[[3]]),  # nb models = nb names of dim 3 of the array
                               length(dimnames(model_runs@models.evaluation@val)[[4]])),  # nb cross-validation runs
                       dimnames = list(cluster = list_clust,
                                       model = dimnames(model_runs@models.evaluation@val)[[3]],
                                       cv.run = dimnames(model_runs@models.evaluation@val)[[4]]))

# Create an array to store optimal cutoffs of conversion to presence-absence
# to use them in the ensemble model
allclust_cutoffs <- allclust_jaccard

jaccard_list <- list()


for (i in 1:length(list_clust)) {
  
  clst <- paste0("cluster", i)
  cat(paste("----", Sys.time(), "cluster", i, " evaluation initialised ----\n", sep = " "))
  
  setwd(paste0("./outputs/modelling/models/cluster", i))
  
  # Get presence, calibration and evaluation coordinates
  load("model_runs")
  load(model_runs@formated.input.data@link)  # presence
  load(model_runs@calib.lines@link)  # calibration
  load(model_runs@models.prediction@link)  # evaluation
  
  for (cv in dimnames(model_runs@models.evaluation@val)[[4]]) {
    
    # Observed data (P/A) that were not used for calibration
    obs_data <- data@data.species[which(!calib.lines[, paste0("_", cv), "_AllData"])]

    # Get the lines that were used for calibration
    cur_eval <- calib.lines[, paste0("_", cv), "_AllData"] 
    cur_eval <- cur_eval[which(!is.na(cur_eval))]  # Pas de NAs
    
    
    for (mod in dimnames(model_runs@models.evaluation@val)[[3]]) {
      
      # Evaluation of only the lines that were not used for calibration
      cur_preds <- models.prediction[, mod, cv, "AllData"][!cur_eval] 
    
      # Sometimes the run fails and we have NAs in the predictions, check if this is the case
      if(!any(is.na(cur_preds))) { 
      
        jaccard_test <- NULL
        # Calculate the Jaccard index for all thresholds between 0 and 1000 with a step of 1
        for(cutoff in seq(0, 1000, by = 1)) {
          
          pred_pa <- cur_preds  # predictions as presence/absence
          pred_pa[pred_pa < cutoff] <- 0 # absence
          pred_pa[pred_pa >= cutoff] <- 1 # presence
          
          TP <- length(which(obs_data == 1 & pred_pa == 1))  # nb true positive
          FN <- length(which(obs_data == 1 & pred_pa == 0))  # nb false negative
          FP <- length(which(obs_data == 0 & pred_pa == 1))  # nb false positive
          
          # Jaccard index formula
          jaccard <- TP / (TP + FP + FN)
          jaccard_test <- rbind.data.frame(jaccard_test,
                                           data.frame(cutoff = cutoff,
                                                      TP = TP,
                                                      FN = FN,
                                                      FP = FP,
                                                      jaccard = jaccard))
        }
        
        # Store in the list
        jaccard_list[[clst]][[paste0("AllData_", cv, "_", mod)]] <- jaccard_test
        
        # Keep the cutoff with the maximum Jaccard index
        # If there are several cutoffs with the same Jaccard index, take the mean
        allclust_cutoffs[clst, mod, cv] <- mean(jaccard_test$cutoff[which(jaccard_test$jaccard == max(jaccard_test$jaccard))])
        
        # Keep this index value for model evaluation
        allclust_jaccard[clst, mod, cv] <- jaccard_test$jaccard[which(jaccard_test$cutoff == round(allclust_cutoffs[clst, mod, cv]))]
      } 
      else {   # If there are NAs in cur_preds
        
        jaccard_list[[clst]][[paste0("AllData_", cv, "_", mod)]] <- NA
        
        allclust_cutoffs[clst, mod, cv] <- NA 
        
        allclust_jaccard[clst, mod, cv] <- NA
      }
    }
  }
  
  setwd(initial.wd)
}



# Save cutoffs
saveRDS(allclust_cutoffs, file = "./outputs/evaluation/jaccard_cutoffs.RDS")
# Save Jaccard evaluations
saveRDS(allclust_jaccard, file = "./outputs/evaluation/jaccard_evals.RDS")
saveRDS(jaccard_list, file = "./outputs/evaluation/jaccard_tests.RDS")



## Plot ---

allclust_jaccard <- readRDS("./outputs/evaluation/jaccard_evals.RDS")

ggjaccard <- reshape2::melt(allclust_jaccard) 

png("./outputs/plots/model_evaluation/jaccard_evals.png", width = 2200, height = 1600, res = 300)
ggplot(ggjaccard, aes(x = model, y = value, col = cv.run)) +
  geom_point() + 
  facet_wrap (~ cluster) +
  labs(col = "Cross validation runs") +
  xlab("Models") +
  ylab("Jaccard index") +
  theme(legend.title = element_text(size = 10)) 
dev.off()




##  3. AUC -----------------------------------------

# Empty list
auc_by_model <- list()

for (clst in list_clust) {
  
  load(paste0("./outputs/modelling/models/", clst, "/model_runs"))
  
  biomod_eval <- get_evaluations(model_runs, as.data.frame = TRUE)

  # Keep AUC (ROC) values
  auc_data <- biomod_eval[biomod_eval$Eval.metric == "ROC", ]
  
  # SDM names
  model_names <- unique(gsub("_RUN\\d+_AllData", "", auc_data$Model.name))
  
  # Add to list
  for (model in model_names) {
    auc_by_model[[clst]][[model]] <- auc_data$Testing.data[grep(model, auc_data$Model.name)]
  }
  
}

saveRDS(auc_by_model, file = "./outputs/evaluation/auc_evals.RDS")
