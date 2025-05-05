###-------------------------------------------------###
###              SELECTION OF VARIABLES             ###
###                                                 ###
###            Environmental variables              ###
###-------------------------------------------------###


library(virtualspecies)
library(tidyr)
library(Rarity)
library(ggplot2)
library(biomod2)
library(terra)


## Data downloading
baseline_30 <- brick("./inputs/envdata/baseline/processed/env_baseline_1deg.tif")
mycto_clusters <- readRDS("./outputs/bioregionalisation/mycto_clusters.RDS")
clust_colors <- readRDS("./outputs/bioregionalisation/clust_colors.RDS")

## Variables distribution
hist(baseline_30)  # if gaussian, then we will use the Pearson coefficient


## 1. Correlation between variables -----------------------------------------
baseline_30 <- terra::rast(baseline_30)

png("./outputs/plots/var_selection/collinearity_groups.png", width = 1800, height = 1400, res = 300)
groups <- removeCollinearity(baseline_30, multicollinearity.cutoff = 0.7,
                             method = "spearman",
                             plot = T)
dev.off()


## Store the selected variables for all clusters
clust_names <- paste0("cluster", unique(mycto_clusters$lvl1))

clust_var <- lapply(clust_names,
                    function (x) c("temperature200", "salinity200", "oxygen200", "chla2.5", "bathymetry",
                                   "chla200", "oxygen500", "salinity500", "chla100"))
names(clust_var) <- clust_names


## Save
saveRDS(clust_var, "./outputs/modelling/var_to_select.RDS")



## 2. Variable importance -----------------------------------------

## /!\ We must have ran the models in script 05 with the "Before variable selection" parameters ##


clust_var <- readRDS("./outputs/modelling/var_to_select.RDS")
list_clust <- c("cluster1", "cluster2", "cluster3", "cluster4")

## Select important variables
sel_var <- list()

# Create list of plots
plot_list <- list()

for (clst in list_clust) {
  
  # Download model data
  load(paste0("./outputs/modelling/var_selection/", clst, "/model_runs"))
  
  # get variable importance values
  varimp <- reshape2::melt(get_variables_importance(model_runs))
  colnames(varimp) <- c("Variable", "Model", "CV.Run", "PA.Run", "Variable.importance")
  varimp$Variable <- reorder(varimp$Variable,    # Order by importance median
                             varimp$Variable.importance,
                             FUN = median,
                             na.rm = TRUE)
  
  # Plot
  p <- ggplot(varimp, aes(x = Variable, y = Variable.importance)) +
    geom_boxplot(fill = clust_colors$color[clust_colors$cluster == clst]) + 
    theme_bw() +
    ggtitle(paste("Region", substr(clst, 8, 8))) +
    theme(axis.text.x = element_text(size = 20, angle = 60, hjust = 1),
          axis.text.y = element_text(size = 20),
          axis.title = element_text(size = 23),
          plot.title = element_text(size = 25))
  
  png(file = paste0("./outputs/plots/var_selection/variable_importance_", clst, ".png"), width = 2000, height = 1400, res = 200)
  print(p)
  dev.off()
  
  # Add plot to list
  plot_list[[clst]] <- p
  
  # Get median values for all models
  median_imp <- aggregate(Variable.importance ~ Variable, data = varimp, FUN = median)
  moy_imp <- aggregate(Variable.importance ~ Variable, data = varimp, FUN = mean) # Added for specific case of equal median values
  
  # Select variables with importance >= 0.1 for at least 50% of the models
  imp <- as.character(median_imp$Variable[which(median_imp$Variable.importance >= 0.1)]) 
  sel_var[[clst]] <- imp
  
  # If several variables, check if there are negative interactions
  var_interactions <- spread(varimp, Variable, Variable.importance)
  
  png(file = paste0("./outputs/plots/var_selection/variable_interactions_", clst, ".png"), width = 700, height = 500)
  corPlot(var_interactions[, -c(1:3)], method = "pearson")
  dev.off()

  
  if (length(imp) > 1) { 
    corr <- as.data.frame(cor(var_interactions[, imp], use = "na.or.complete"))
    # Get index of negatively correlated variables 
    ind <- which(corr <= -0.3, arr.ind = TRUE)
    
    # Between two variables negatively correlated, keep the most important variable
    while (length(which(corr <= -0.3, arr.ind = TRUE)) != 0) {
      
      # Remove the negatively correlated variables from the list of selected variables
      sel_var[[clst]] <- sel_var[[clst]][-which(sel_var[[clst]] %in% rownames(ind))]
       
      # Get column and row names of negatively correlated variables 
      neg_corr <- as.data.frame(cbind(rownames(corr)[ind[,"row"]], colnames(corr)[ind[,"col"]])) 
      
      for (i in 1:nrow(neg_corr)) {
        # For each couple, select the most important variable
        most_imp <- median_imp$Variable[median_imp$Variable.importance == 
                                                max(median_imp$Variable.importance[median_imp$Variable == neg_corr[i, 1]],
                                                    median_imp$Variable.importance[median_imp$Variable == neg_corr[i, 2]])]
        # Specific case where equal median values, select the highest mean value
        if (length(most_imp) > 1) {
          most_imp <- moy_imp$Variable[moy_imp$Variable.importance == 
                                          max(moy_imp$Variable.importance[moy_imp$Variable == neg_corr[i, 1]],
                                              moy_imp$Variable.importance[moy_imp$Variable == neg_corr[i, 2]])]
        }
        
        neg_corr[i, 3] <- most_imp
        
      }
      
      sel_var[[clst]] <- append(sel_var[[clst]], as.character(unique(neg_corr[, 3])))
      
      # Repeat the operation if the selected variables are negatively correlated
      corr <- cor(dplyr::select(var_interactions, unique(neg_corr[, 3])), use = "na.or.complete")
      ind <- which(corr <= -0.3, arr.ind = TRUE)
    }
  }
}

## Save
save(sel_var, file = paste0("./outputs/modelling/selected_variables.RDS"))


## Create a plot with all variable importance plots
png(file = "./outputs/plots/var_selection/variable_importance_all.png", width = 2000, height = 2000, res = 200)
gridExtra::grid.arrange(grobs = plot_list, ncol = 2)
dev.off()



