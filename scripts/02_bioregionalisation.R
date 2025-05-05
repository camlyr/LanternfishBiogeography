###-------------------------------------------------###
###               BIOREGIONALISATION                ###
###                                                 ###
###            Clusters, network, analyses          ###
###                                                 ###
###-------------------------------------------------###


library(raster)
library(dplyr)
library(reshape2)
library(biogeonetworks)
library(RColorBrewer)
library(stringr)
library(tmap)
library(ggplot2)
library(sf)


## 1. Data preparation -----------------------------------------

## Load data
sp_stack <- stack("./inputs/occdata/processed/mycto_occ_1deg.tif")
sp_occ <- read.csv2("./inputs/occdata/processed/mycto_occ_1deg.csv")
baseline_30 <- brick("./inputs/envdata/baseline/processed/env_baseline_1deg.tif")
land <- st_read("./inputs/sigdata/ne_50m_land.shp")


# Group data into sites (raster cells)
mycto_df <- getValues(sp_stack)

mycto_df <- mycto_df %>% 
  melt() %>%
  na.omit() %>%
  dplyr::select(Var1, Var2) %>%
  rename("cell" = "Var1",
         "scientificName" = "Var2")


## 2. Write the network in Pajek format -----------------------------------------

writePajek(mycto_df, 
           site.field = "cell", 
           species.field = "scientificName", 
           abundance.field = NULL,
           filename = "./outputs/bioregionalisation/mycto_clusters.net")



## 3. Clusters research by Map Equation algorithm -----------------------------------------

# /!\ Install Infomap following https://www.mapequation.org/infomap/

# Arguments : Infomap, network, output, tree, nb of runs
system("./Infomap ./outputs/bioregionalisation/mycto_clusters.net ./outputs/bioregionalisation/ --tree -N 100") 



## 4. Read Map Equation clusters -----------------------------------------

mycto_clusters <- readInfomapTree("./outputs/bioregionalisation/mycto_clusters.tree",
                                  network.summary = TRUE,
                                  replace.leaf.names = FALSE) # Changes node numbers to actual names for terminal nodes (i.e. site & species names)



## 5. Analyse the results -----------------------------------------

# Number of elements per cluster
count(mycto_clusters, lvl1)

# Number of sites/species per cluster
sites <- getSiteTable(db = mycto_df,
                      site.field = 'cell',
                      network = mycto_clusters)
species <- getSpeciesTable(db = mycto_df,
                           species.field = 'scientificName',
                           network = mycto_clusters)

count(sites, lvl1)
count(species, lvl1)

# Color the small clusters
# mycto_clusters <- mycto_clusters[mycto_clusters$lvl1 %in% c(6, 7, 8, 9, 10, 11), ]


# Color the main clusters
mycto_clusters <- attributeColors(network = mycto_clusters,
                                  nb.max.colors = 5,   # nb of significant clusters
                                  lvl = "lvl1",
                                  # cluster.order = "sites+species", # order to select the most significant clusters
                                  sh.grey = TRUE,  # other clusters will be colored in shades of grey
                                  db = mycto_df,
                                  site.field = "cell",
                                  species.field = "scientificName")

# Change the colors of the main clusters
mycto_clusters$color <- ifelse(mycto_clusters$lvl1 == 1, "#D55E00",
                               ifelse(mycto_clusters$lvl1 == 2, "#56B4E9",
                                      ifelse(mycto_clusters$lvl1 == 3, "#F0E442",
                                             ifelse(mycto_clusters$lvl1 == 4, "#009E73", 
                                                    ifelse(mycto_clusters$lvl1 == 5, "#CC79A7", 
                                                           mycto_clusters$color)))))

## 6. Write the network in GDF format for Gephi -----------------------------------------

writeGDF(db = mycto_df, 
         site.field = "cell", 
         species.field = "scientificName",
         color.field = "color",
         network = mycto_clusters, 
         filename = "./outputs/bioregionalisation/mycto_clusters.gdf") 


## 7. Maps of clusters -----------------------------------------

# Get the colors of the cells
mycto_cells <- getSiteTable(db = mycto_df, 
                            site.field = "cell", 
                            network = mycto_clusters)

# Get the colors of the clusters
cols <- unique(mycto_cells[, c("lvl1", "color")])
cols$ID <- as.numeric(cols$lvl1)

# Make an empty raster based on the species rasters
rasterclust <- raster(sp_stack, layer = 0)
# Assign the cluster to which each cell belongs
mycto_cells$cellnumbers <- as.numeric(as.character(mycto_cells$Name))
rasterclust[mycto_cells$cellnumbers] <- as.numeric(mycto_cells$lvl1)
# Transform the raster into a categorical raster
rasterclust <- ratify(rasterclust)
levels(rasterclust)[[1]] <-  merge(levels(rasterclust)[[1]], cols)
names(rasterclust) <- "Bioregions"

# Latitude and longitude lines
graticules = st_graticule(ndiscr = 10000, 
                          lat = seq(-90, -30, 20),
                          lon = seq(0, 360, 30)) %>%
  st_transform("+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") %>% 
  st_geometry()


# Create the map
clusters_map <- tm_shape(rasterclust) +
  tm_raster(col = "Bioregions",
            palette = levels(rasterclust)[[1]]$color,
            legend.show = F) +
  tm_shape(land) +
  tm_borders(col = "grey", lwd = 8) +
  tm_shape(graticules) +
  tm_lines(col = "grey", alpha = 0.5) +
  # Add legend manually
  tm_add_legend("fill", title = "Regions",
                labels = seq(1, 11, 1), 
                col = levels(rasterclust)[[1]]$color[order(as.numeric(levels(levels(rasterclust)[[1]]$lvl1)))]) +
  tm_legend(legend.title.size = 4, legend.text.size	= 3.3, legend.position = c("left","bottom"), scale = 0.3)


# Save map
tmap_save(tm = clusters_map,
          filename = "./outputs/plots/bioregionalisation/clusters_map.png")



## 8. Participation coefficient -----------------------------------------

mycto_clusters <- participationCoefficient(network = mycto_clusters,
                                           db = mycto_df, 
                                           site.field = "cell", 
                                           species.field = "scientificName",
                                           lvl = "lvl1") 

# Re-extract the cell table to get the participation coefficient column
mycto_cells <- getSiteTable(db = mycto_df, 
                            site.field = "cell",
                            network = mycto_clusters)

# Make an empty raster based on the species rasters
rasterPC <- raster(sp_stack, layer = 0)
# Assign the participation coefficient to each cell
mycto_cells$cellnumbers <- as.numeric(as.character(mycto_cells$Name))
rasterPC[mycto_cells$cellnumbers] <- mycto_cells$participation.coef
names(rasterPC) <- "Participation coefficient"

# Create the map
part_coeff <- tm_shape(rasterPC) +
  tm_raster() +
  tm_shape(land) +
  tm_borders(col = "grey") +
  tm_legend(legend.title.size = 2, legend.text.size	= 2)

# Save map
png(file = "./outputs/plots/bioregionalisation/partcoeff_map.png", width = 2000, height = 1200)
part_coeff
dev.off()



## 9. Cluster metrics -----------------------------------------

mycto_metrics <- clusterMetrics(db = mycto_df, 
                                network = mycto_clusters, 
                                site.field = "cell", 
                                species.field = "scientificName",
                                level = "lvl1")

# Endemism proportion per cluster
endemism <- dplyr::select(mycto_metrics$region.stats, cluster, char.richness, end.richness)
endemism$end.prop <- endemism$end.richness/endemism$char.richness

fidelity <- do.call(data.frame, aggregate(mycto_metrics$species.stats$Occ.Fi,
                                          by = list(mycto_metrics$species.stats$cluster),
                                          FUN = function(x) c(mean = mean(x), median = median(x), min = min(x), max = max(x))))


# Fidelity and frequency of species in the five clusters
count_rec <- as.data.frame(table(sp_occ$scientificName))
colnames(count_rec) <- c("species", "nb.records")

cols <- cols %>%
  select(lvl1, color) %>%
  rename("cluster" = "lvl1")

count_rec <- count_rec %>% 
  merge(dplyr::select(mycto_metrics$species.stats, species, cluster, Occ.Fi)) %>%
  merge(cols) %>%
  filter(cluster %in% c(1, 2, 3, 4, 5))

# Fidelity plot
fidelity_freq <- ggplot(data = count_rec) +
  geom_point(aes(x = nb.records, y = Occ.Fi), size = 6, col = count_rec$color) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 25)) +
  ylab("Species fidelity") + 
  xlab("Number of records") +
  ggtitle("Fidelity and frequency of species in the five clusters")

# Save plot
png(file = "./outputs/plots/bioregionalisation/fidelity_freq.png", width = 1600, height = 1000, , res = 100)
fidelity_freq
dev.off()


# Fidelity boxplot
my_colors <- mycto_clusters$color
names(my_colors) <- mycto_clusters$lvl1

clust_fid <- dplyr::select(mycto_metrics$species.stats, cluster, Occ.Fi)
clust_fid <- clust_fid[order(clust_fid$cluster), ]

fidelity_clust <- ggplot(data = clust_fid) +
  geom_boxplot(aes(x = Occ.Fi, y = reorder(cluster, Occ.Fi, median), fill = cluster)) +
  scale_fill_manual(values = my_colors) +
  ylab("Clusters") +
  xlab("Species fidelity") +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 25),
        legend.position = "none") +
  ggtitle("Species fidelity to each cluster")

# Save plot
png(file = "./outputs/plots/bioregionalisation/fidelity_clust.png", width = 2000, height = 1400, res = 220)
fidelity_clust
dev.off()




## Save data -----------------------------------------

# Save dataframe
saveRDS(mycto_clusters, "./outputs/bioregionalisation/mycto_clusters.RDS")

# Save raster
saveRDS(rasterclust, "./outputs/bioregionalisation/clust_raster.RDS")


## Assign colors to clusters -----------------------------------------

# Create a reference dataframe for cluster colors
clust_colors <- unique(mycto_clusters[, c("lvl1", "color")])
colnames(clust_colors) <- c("cluster", "color")
clust_colors$cluster <- paste0("cluster", clust_colors$cluster)

saveRDS(clust_colors, "./outputs/bioregionalisation/clust_colors.RDS")



## 10. Prepare cluster data for modelling -----------------------------------------

## Extract coordinates and id
clust_df <- as.data.frame(rasterToPoints(rasterclust)) 
names(clust_df) <- c("x", "y", "id")

## Get the matching cluster numbers
clust_df$cluster <- rasterclust@data@attributes[[1]]$lvl1[match(clust_df$id, rasterclust@data@attributes[[1]]$ID)]

## Add matching environmental data
clust_df <- data.frame(clust_df, 
                       raster::extract(baseline_30, clust_df[, c("x", "y")]))

## Remove sites with no environmental data
clust_df <- clust_df[-unique(which(is.na(clust_df), arr.ind = TRUE)[, 1]), ]


## Presence-absence information for each cluster
for (clust in sort(unique(clust_df$cluster))) {
  clust_df[paste0("cluster", clust)] <- ifelse(clust_df$cluster == clust, 1, 0)
}

## Convert to sf object
clust_sf <- st_as_sf(clust_df,
                     coords = c("x", "y"),
                     crs = as.character(crs(baseline_30)))


## Save
saveRDS(clust_df, "./outputs/bioregionalisation/clust_df.RDS")
saveRDS(clust_sf, "./outputs/bioregionalisation/clust_sf.RDS")




