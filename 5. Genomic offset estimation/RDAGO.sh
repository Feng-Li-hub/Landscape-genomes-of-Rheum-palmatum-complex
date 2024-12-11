########RDA——GO
##pan_all
library(pegas)
library(ggplot2)
library(raster)
library(rgdal)
library(LEA)
library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)
library(ggpubr)
library(vegan)
library(qvalue)
library(robust)
library(WMDB)
library(ggVennDiagram)
library(cowplot)
library(corrplot)
library(rgeos)
library(data.table)
AllFreq <- read.table("pan_loci2_AF.csv", sep = ",", header = T, row.names=1)

## Environmental variable processing
# Historical environment
library(raster)

# Get the list of files
file_list <- list.files("C:\\Users\\feng\\Desktop\\景观\\RDA_GO\\RDA_GO_171_3144\\历史\\适生区asc", pattern = ".asc$", full.names = TRUE)

# Load each file and check
raster_list <- lapply(file_list, raster)

# Check the attributes of each raster
for (i in seq_along(raster_list)) {
  cat("File", file_list[i], "has", nrow(raster_list[[i]]), "rows and", ncol(raster_list[[i]]), "columns.\n")
}

# Try to stack the files
tryCatch({
  ras <- stack(raster_list)
}, error = function(e) {
  message("Error:", e)
})

# Function to remove NAs from stack
remove.NAs.stack <- function(rast.stack) {
  nom <- names(rast.stack)
  test1 <- calc(rast.stack, fun=sum)
  test1[!is.na(test1)] <- 1
  test2 <- rast.stack * test1
  test2 <- stack(test2)
  names(test2) <- nom
  return(test2)
}

ras <- remove.NAs.stack(ras) # Apply the remove.NAs.stack function to the raster stack created earlier

### 585
### 585_50
# Get the list of files
file_list <- list.files("C:\\Users\\feng\\Desktop\\景观\\RDA_GO\\RDA_GO_171_3144\\ssp585_50\\58550_适生区_asc", pattern = ".asc$", full.names = TRUE)

# Load each file and check
raster_list <- lapply(file_list, raster)

# Check the attributes of each raster
for (i in seq_along(raster_list)) {
  cat("File", file_list[i], "has", nrow(raster_list[[i]]), "rows and", ncol(raster_list[[i]]), "columns.\n")
}

# Try to stack the files
tryCatch({
  ras_58550 <- stack(raster_list)
}, error = function(e) {
  message("Error:", e)
})

names(ras_58550)
ras_58550 <- remove.NAs.stack(ras_58550)

#### 58590
# Get the list of files
file_list <- list.files("C:\\Users\\feng\\Desktop\\景观\\RDA_GO\\RDA_GO_171_3144\\ssp585_90\\58590_适生区_asc", pattern = ".asc$", full.names = TRUE)

# Load each file and check
raster_list <- lapply(file_list, raster)

# Check the attributes of each raster
for (i in seq_along(raster_list)) {
  cat("File", file_list[i], "has", nrow(raster_list[[i]]), "rows and", ncol(raster_list[[i]]), "columns.\n")
}

# Try to stack the files
tryCatch({
  ras_58590 <- stack(raster_list)
}, error = function(e) {
  message("Error:", e)
})

names(ras_58590)
ras_58590 <- remove.NAs.stack(ras_58590)

Coordinates <- read.table("43pop_longlat.csv", sep = ",", header = TRUE, row.names = 1)

# Match coordinate data based on AllFreq's row names
Coordinates <- Coordinates[match(row.names(AllFreq), Coordinates$id2, nomatch = 0),]

# Rename columns of the coordinate data
colnames(Coordinates) <- c("Population", "Longitude", "Latitude", "Elevation")

# Extract environmental variables from raster 'ras' corresponding to the coordinates and create a dataframe
Env <- data.frame(extract(ras, Coordinates[, 2:3]))

## Standardize the variables
Env <- scale(Env, center = TRUE, scale = TRUE)  # center=TRUE and scale=TRUE are default settings for scale() function

## Recover the scaling coefficients
scale_env <- attr(Env, 'scaled:scale')  # Get the scaling factor
center_env <- attr(Env, 'scaled:center')  # Get the center values

## Convert the standardized environmental data back into a dataframe
Env <- as.data.frame(Env)
row.names(Env) <- c(Coordinates$Population)  # Set the row names as population names

# Source other R scripts
source("C:\\Users\\feng\\Desktop\\景观\\RDA_GO\\RDA_GO_171_3144\\输入脚本\\genomic_offset.R")
source("C:\\Users\\feng\\Desktop\\景观\\RDA_GO\\RDA_GO_171_3144\\输入脚本\\provgar_offset.R")

# Perform Principal Component Analysis (PCA)
pca <- rda(AllFreq, scale = TRUE)

# Plot the PCA eigenvalue plot
screeplot(pca, type = "barplot", npcs = 10, main = "PCA Eigenvalues")

### Extract the first few principal components (PCs)
PCs <- scores(pca, choices = c(1:3), display = "sites", scaling = 0)

# Create a structure dataframe containing populations and principal components
PopStruct <- data.frame(POP = row.names(AllFreq), PCs)
colnames(PopStruct) <- c("POP", "PC1", "PC2", "PC3") # Rename columns
write.csv(PopStruct, file = "PopStruct.csv", quote = FALSE)

# Combine coordinates, population structure, and environmental variables into one dataframe
Variables <- data.frame(Coordinates, PopStruct[,-1], Env)

# Perform Redundancy Analysis (RDA)
RDA_outliers <- rda(AllFreq ~ BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15, Variables)

## Run the genomic_offset function for the specified populations (58550 and 58590)
res_RDA_proj58590 <- genomic_offset(RDA_outliers, K = 3, env_pres = ras, env_fut = ras_58590, range = NULL, method = "loadings", scale_env = scale_env, center_env = center_env)
res_RDA_proj58550 <- genomic_offset(RDA_outliers, K = 3, env_pres = ras, env_fut = ras_58550, range = NULL, method = "loadings", scale_env = scale_env, center_env = center_env)

## Table for global genetic offset predicted for 58550 and 58590
RDA_proj_offset <- data.frame(rbind(rasterToPoints(res_RDA_proj58550$Proj_offset_global), rasterToPoints(res_RDA_proj58590$Proj_offset_global)), Date = c(rep("58550", nrow(rasterToPoints(res_RDA_proj58550$Proj_offset_global))), rep("58590", nrow(rasterToPoints(res_RDA_proj58590$Proj_offset_global)))))
write.csv(RDA_proj_offset, file = "RDA_proj_offset.csv", quote = F)

# Read the RDA projected offset results for population 58550
RDA_proj_offset_58550 <- read.csv("RDA_proj_offset-58550.csv")
raster_data <- rasterFromXYZ(RDA_proj_offset_58550[, c("x", "y", "Global_offset")])
asc_file <- "RDA_proj_offset_58550.asc"
writeRaster(raster_data, filename = asc_file, format = "ascii")

# Read the RDA projected offset results for population 58590
RDA_proj_offset_58590 <- read.csv("RDA_proj_offset-58590.csv")
raster_data <- rasterFromXYZ(RDA_proj_offset_58590[, c("x", "y", "Global_offset")])
asc_file <- "RDA_proj_offset_58590.asc"
writeRaster(raster_data, filename = asc_file, format = "ascii")



####126
## Climate data for 12690
# Get the list of files
file_list <- list.files("C:\\Users\\feng\\Desktop\\景观\\RDA_GO\\RDA_GO_171_3144\\ssp126_90\\12690_适生区", pattern = ".tif$", full.names = TRUE)
# Load each file and check
raster_list <- lapply(file_list, raster)
# Check the attributes of each raster
for (i in seq_along(raster_list)) {
  cat("File", file_list[i], "has", nrow(raster_list[[i]]), "rows and", ncol(raster_list[[i]]), "columns.\n")
}
# Try stacking the files
tryCatch({
  ras_12690 <- stack(raster_list)
}, error = function(e) {
  message("Error:", e)
})
names(ras_12690)
ras_12690 <- remove.NAs.stack(ras_12690)

####12650
## Climate data for 12650
# Get the list of files
file_list <- list.files("C:\\Users\\feng\\Desktop\\景观\\RDA_GO\\RDA_GO_171_3144\\ssp126_50\\12650_适生区", pattern = ".tif$", full.names = TRUE)
# Load each file and check
raster_list <- lapply(file_list, raster)
# Check the attributes of each raster
for (i in seq_along(raster_list)) {
  cat("File", file_list[i], "has", nrow(raster_list[[i]]), "rows and", ncol(raster_list[[i]]), "columns.\n")
}
# Try stacking the files
tryCatch({
  ras_12650 <- stack(raster_list)
}, error = function(e) {
  message("Error:", e)
})
names(ras_12650)
ras_12650 <- remove.NAs.stack(ras_12650)

# Perform genomic offset analysis for 12690 and 12650
res_RDA_proj12690 <- genomic_offset(RDA_outliers, K = 3, env_pres = ras, env_fut = ras_12690, range = NULL, method = "loadings", scale_env = scale_env, center_env = center_env)
res_RDA_proj12650 <- genomic_offset(RDA_outliers, K = 3, env_pres = ras, env_fut = ras_12650, range = NULL, method = "loadings", scale_env = scale_env, center_env = center_env)

# Create a table for the global genetic offset for 12650 and 12690
RDA_proj_offset <- data.frame(rbind(rasterToPoints(res_RDA_proj12650$Proj_offset_global), rasterToPoints(res_RDA_proj12690$Proj_offset_global)), Date = c(rep("12650", nrow(rasterToPoints(res_RDA_proj12650$Proj_offset_global))), rep("12690", nrow(rasterToPoints(res_RDA_proj12690$Proj_offset_global)))))
write.csv(RDA_proj_offset, "RDA_proj_offset_12650_90.csv", row.names = FALSE)

# Read the offset results for 12650 and convert to raster
RDA_proj_offset_12650 <- read.csv("RDA_proj_offset_12650.csv")
raster_data <- rasterFromXYZ(RDA_proj_offset_12650[, c("x", "y", "Global_offset")])
asc_file <- "RDA_proj_offset_12650.asc"
writeRaster(raster_data, filename = asc_file, format = "ascii")

# Read the offset results for 12690 and convert to raster
RDA_proj_offset_12690 <- read.csv("RDA_proj_offset_12690.csv")
raster_data <- rasterFromXYZ(RDA_proj_offset_12690[, c("x", "y", "Global_offset")])
asc_file <- "RDA_proj_offset_12690.asc"
writeRaster(raster_data, filename = asc_file, format = "ascii")

########## 43pop_future_climate_GO
library(dplyr)
library(fuzzyjoin)

# Read the latitude and longitude data file
df_pop <- read.csv("43pop.csv", header = TRUE)

# Read the GO values data file
df_go <- read.csv("RDA_proj_offset_12650.csv", header = TRUE)

# Perform a fuzzy join on the two datasets by latitude and longitude
merged_data <- df_pop %>%
  fuzzy_inner_join(df_go, 
                   by = c("longitude" = "x", "latitude" = "y"),
                   match_fun = list(`>=`, `<=`))  # Use "greater than or equal" and "less than or equal" for fuzzy matching

# Calculate the distance between each population and raster point
merged_data <- merged_data %>%
  mutate(distance = sqrt((longitude - x)^2 + (latitude - y)^2))

# Save the merged data to a CSV file
write.csv(merged_data, "43_12650.csv", row.names = FALSE)

# Find the row with the minimum distance for each population and select the corresponding GO value
min_distance_rows <- merged_data %>%
  group_by(pop) %>%
  slice(which.min(distance)) %>%
  select(pop, Global_offset)

# Save the minimum distance rows to a CSV file
write.csv(min_distance_rows , "43_12650_min.csv", row.names = FALSE)


# Read the GO value data file
df_go <- read.csv("RDA_proj_offset_12690.csv", header = TRUE)
merged_data <- df_pop %>%
     fuzzy_inner_join(df_go, 
                      by = c("longitude" = "x", "latitude" = "y"),
                      match_fun = list(`>=`, `<=`))  # Use "greater than or equal" and "less than or equal" for fuzzy matching					  
merged_data <- merged_data %>%
     mutate(distance = sqrt((longitude - x)^2 + (latitude - y)^2))
write.csv(merged_data, "43_12690.csv", row.names = FALSE)
# Find the row with the minimum distance for each population and select the corresponding GO value
min_distance_rows <- merged_data %>%
  group_by(pop) %>%
  slice(which.min(distance)) %>%
  select(pop, Global_offset)
write.csv(min_distance_rows , "43_12690_min.csv", row.names = FALSE)

# Read the GO value data file
df_go <- read.csv("RDA_proj_offset-58550.csv", header = TRUE)
merged_data <- df_pop %>%
     fuzzy_inner_join(df_go, 
                      by = c("longitude" = "x", "latitude" = "y"),
                      match_fun = list(`>=`, `<=`))  # Use "greater than or equal" and "less than or equal" for fuzzy matching					  
merged_data <- merged_data %>%
     mutate(distance = sqrt((longitude - x)^2 + (latitude - y)^2))
write.csv(merged_data, "43_58550.csv", row.names = FALSE)
# Find the row with the minimum distance for each population and select the corresponding GO value
min_distance_rows <- merged_data %>%
  group_by(pop) %>%
  slice(which.min(distance)) %>%
  select(pop, Global_offset)
write.csv(min_distance_rows , "43_58550_min.csv", row.names = FALSE)

# Read the GO value data file
df_go <- read.csv("RDA_proj_offset-58590.csv", header = TRUE)
merged_data <- df_pop %>%
     fuzzy_inner_join(df_go, 
                      by = c("longitude" = "x", "latitude" = "y"),
                      match_fun = list(`>=`, `<=`))  # Use "greater than or equal" and "less than or equal" for fuzzy matching					  
merged_data <- merged_data %>%
     mutate(distance = sqrt((longitude - x)^2 + (latitude - y)^2))
write.csv(merged_data, "43_58590.csv", row.names = FALSE)
# Find the row with the minimum distance for each population and select the corresponding GO value
min_distance_rows <- merged_data %>%
  group_by(pop) %>%
  slice(which.min(distance)) %>%
  select(pop, Global_offset)
write.csv(min_distance_rows , "43_58590_min.csv", row.names = FALSE)


library(sp)
library(raster)
DH5026_raster <- raster("RDA_proj_offset_12650.asc")
DH9026_raster <- raster("RDA_proj_offset_12690.asc")
DH5085_raster <- raster("RDA_proj_offset_58550.asc")
DH9085_raster <- raster("RDA_proj_offset_58590.asc")

pdf(file = "panrdaDH_all_GO.pdf")
par(mfrow = c(2, 2), mar = c(2, 2, 2, 1), oma = c(1, 1, 1, 2))
# First subplot
plot(DH5026_raster, zlim = c(0.01, 2), legend = FALSE)  # Plot DH5026_raster, set z-axis limits, hide legend
points_data <- read.csv("pan_rdago.csv")  # Read GO value data from CSV file
library(sp)  # Load spatial data handling library
# Create spatial points data frame
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object
# Set coordinate reference system (assume WGS84 for longitude and latitude)
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
# Calculate maximum and minimum GO values
max_GO <- max(points_data$GO_5026)  # Get maximum GO value
min_GO <- min(points_data$GO_5026)  # Get minimum GO value
# Calculate symbol radius
points_data$radius <- 0.05 + 0.3 * (points_data$GO_5026 - min_GO) / (max_GO - min_GO)
# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(a)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label
# Second subplot
plot(DH9026_raster, zlim = c(0.01, 2), legend = FALSE)  # Plot DH9026_raster
points_data <- read.csv("pan_rdago.csv")  # Read GO value data from CSV file
# Create spatial points data frame
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")  # Set coordinate reference system
# Calculate maximum and minimum GO values
max_GO <- max(points_data$GO_9026)
min_GO <- min(points_data$GO_9026)
points_data$radius <- 0.05 + 0.3 * (points_data$GO_9026 - min_GO) / (max_GO - min_GO)
# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(b)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label
# Third subplot
plot(DH5085_raster, zlim = c(0.01, 2), legend = FALSE)  # Plot DH5085_raster
points_data <- read.csv("pan_rdago.csv")  # Read GO value data from CSV file
# Create spatial points data frame
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")  # Set coordinate reference system
# Calculate maximum and minimum GO values
max_GO <- max(points_data$GO_5085)
min_GO <- min(points_data$GO_5085)
points_data$radius <- 0.05 + 0.3 * (points_data$GO_5085 - min_GO) / (max_GO - min_GO)
# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(c)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label
# Fourth subplot
plot(DH9085_raster, zlim = c(0.01, 2))  # Plot DH9085_raster, display legend
points_data <- read.csv("pan_rdago.csv")  # Read GO value data from CSV file
# Create spatial points data frame
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")  # Set coordinate reference system
# Calculate maximum and minimum GO values
max_GO <- max(points_data$GO_9085)
min_GO <- min(points_data$GO_9085)
points_data$radius <- 0.05 + 0.3 * (points_data$GO_9085 - min_GO) / (max_GO - min_GO)
# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(d)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label
dev.off()

#######################################################################################################################
###Contemporary Adaptive Landscapes

source("C:\\Users\\feng\\Desktop\\景观\\RDA_GO\\RDA_GO_171_3144\\输入脚本\\adaptive_index.R")
library(pegas)
library(ggplot2)
library(raster)
library(rgdal)
library(LEA)
library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)
library(ggpubr)
library(vegan)
library(qvalue)
library(robust)
library(WMDB)
library(ggVennDiagram)
library(cowplot)
library(corrplot)
library(rgeos)
summary(RDA_outliers)
eigens <- eigenvals(RDA_outliers)
eigens <- RDA_outliers$CA$eig
total_variance <- sum(eigens)
variance_explained <- eigens / total_variance * 100
# View the proportion of variance explained by each principal component
variance_explained
total_variance <- sum(eigens)
# Calculate the proportion of variance explained by each principal component
variance_explained <- eigens / total_variance * 100
# Print the results
print(variance_explained)

# Extract species and environmental variable scores from RDA analysis results
TAB_loci <- as.data.frame(scores(RDA_outliers, choices=c(1:2), display="species", scaling="none"))
TAB_var <- as.data.frame(scores(RDA_outliers, choices=c(1:2), display="bp"))

ggplot() +
  # Add horizontal and vertical reference lines
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  # Plot species points
  geom_point(data = TAB_loci, aes(x=RDA1*3, y=RDA2*3), colour = "#EB8055FF", size = 2, alpha = 0.8) +
  # Plot environmental variable arrows
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, 
               arrow=arrow(length = unit(0.02, "npc"))) + 
  # Add labels for environmental variables
  geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  xlab("RDA 1 (Corresponding value)") + ylab("RDA 2 (Corresponding value)") +
  facet_wrap(~"Adaptively enriched RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), 
        strip.text = element_text(size=11))

# Calculate the adaptive index and extract RDA results
res_RDA_proj_current <- adaptive_index(RDA = RDA_outliers, K = 3, env_pres = ras, range = NULL, method = "loadings", scale_env = scale_env, center_env = center_env)
RDA_proj <- list(res_RDA_proj_current$RDA1, res_RDA_proj_current$RDA2)
# Convert RDA results to point data
RDA_proj <- lapply(RDA_proj, function(x) rasterToPoints(x))
# Normalize the values
for(i in 1:length(RDA_proj)){
  RDA_proj[[i]][,3] <- (RDA_proj[[i]][,3]-min(RDA_proj[[i]][,3]))/(max(RDA_proj[[i]][,3])-min(RDA_proj[[i]][,3]))
}
# Combine RDA results into a data frame
TAB_RDA <- as.data.frame(do.call(rbind, RDA_proj[1:2]))
colnames(TAB_RDA)[3] <- "value"
TAB_RDA$variable <- factor(c(rep("RDA1", nrow(RDA_proj[[1]])), rep("RDA2", nrow(RDA_proj[[2]]))), levels = c("RDA1","RDA2"))

ggplot(data = TAB_RDA) + 
  geom_raster(aes(x = x, y = y, fill = cut(value, breaks=seq(0, 1, length.out=10), include.lowest = TRUE))) + 
  scale_fill_viridis_d(alpha = 0.8, direction = -1, option = "A", 
                       labels = c("Negative scores", "", "", "", "Intermediate scores", "", "", "", "Positive scores")) +
  coord_sf(xlim = c(85, 115), ylim = c(20, 45), expand = FALSE) +  
  xlab("Longitude") + 
  ylab("Latitude") +
  guides(fill = guide_legend(title = "Adaptive index")) +
  facet_grid(~ variable) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.grid = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        strip.text = element_text(size = 11))















############################################################################################################
####core
AllFreq <- read.table("core_loci2_AF.csv", sep = ",", header = T, row.names=1)

# Perform Principal Component Analysis (PCA)
pca <- rda(AllFreq, scale = TRUE)
# Plot the eigenvalue graph of PCA
screeplot(pca, type = "barplot", npcs = 10, main = "PCA Eigenvalues")
### Extract the first few principal components (PCs)
PCs <- scores(pca, choices = c(1:3), display = "sites", scaling = 0)
# Create a data frame containing populations and principal components
PopStruct <- data.frame(POP = row.names(AllFreq), PCs)
colnames(PopStruct) <- c("POP", "PC1", "PC2", "PC3") # Rename columns
write.csv(PopStruct, file = "PopStruct.csv", quote = FALSE)
# Combine coordinates, population structure, and environmental variables into one data frame
Variables <- data.frame(Coordinates, PopStruct[,-1], Env)
# Perform Redundancy Analysis (RDA)
RDA_outliers <- rda(AllFreq ~ BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15, Variables)
## Run genomic_offset function for specified populations (58550 and 58590)
res_RDA_proj58590 <- genomic_offset(RDA_outliers, K = 3, env_pres = ras, env_fut = ras_58590, range = NULL, method = "loadings", scale_env = scale_env, center_env = center_env)
res_RDA_proj58550 <- genomic_offset(RDA_outliers, K = 3, env_pres = ras, env_fut = ras_58550, range = NULL, method = "loadings", scale_env = scale_env, center_env = center_env)
## Table global genetic offset predicted for 58550 and 58590
RDA_proj_offset <- data.frame(rbind(rasterToPoints(res_RDA_proj58550$Proj_offset_global), rasterToPoints(res_RDA_proj58590$Proj_offset_global)), Date = c(rep("58550", nrow(rasterToPoints(res_RDA_proj58550$Proj_offset_global))), rep("58590", nrow(rasterToPoints(res_RDA_proj58590$Proj_offset_global)))))
write.csv(RDA_proj_offset,file="RDA_proj_offset.csv",quote=F)


RDA_proj_offset_58550 <- read.csv("RDA_proj_offset_58550.csv")
raster_data <- rasterFromXYZ(RDA_proj_offset_58550[, c("x", "y", "Global_offset")])
asc_file <- "RDA_proj_offset_58550.asc"
writeRaster(raster_data, filename = asc_file, format = "ascii")
 

RDA_proj_offset_58590 <- read.csv("RDA_proj_offset_58590.csv")
raster_data <- rasterFromXYZ(RDA_proj_offset_58590[, c("x", "y", "Global_offset")])
asc_file <- "RDA_proj_offset_58590.asc"
writeRaster(raster_data, filename = asc_file, format = "ascii")



res_RDA_proj12690 <- genomic_offset(RDA_outliers, K = 3, env_pres = ras, env_fut = ras_12690, range = NULL, method = "loadings", scale_env = scale_env, center_env = center_env)
res_RDA_proj12650 <- genomic_offset(RDA_outliers, K = 3, env_pres = ras, env_fut = ras_12650, range = NULL, method = "loadings", scale_env = scale_env, center_env = center_env)
RDA_proj_offset <- data.frame(rbind(rasterToPoints(res_RDA_proj12650$Proj_offset_global), rasterToPoints(res_RDA_proj12690$Proj_offset_global)), Date = c(rep("12650", nrow(rasterToPoints(res_RDA_proj12650$Proj_offset_global))), rep("12690", nrow(rasterToPoints(res_RDA_proj12690$Proj_offset_global)))))
write.csv(RDA_proj_offset, "RDA_proj_offset_12650_90.csv", row.names = FALSE)

RDA_proj_offset_12650 <- read.csv("RDA_proj_offset_12650.csv")
raster_data <- rasterFromXYZ(RDA_proj_offset_12650[, c("x", "y", "Global_offset")])
asc_file <- "RDA_proj_offset_12650.asc"
writeRaster(raster_data, filename = asc_file, format = "ascii")
 
RDA_proj_offset_12690 <- read.csv("RDA_proj_offset_12690.csv")
raster_data <- rasterFromXYZ(RDA_proj_offset_12690[, c("x", "y", "Global_offset")])
asc_file <- "RDA_proj_offset_12690.asc"
writeRaster(raster_data, filename = asc_file, format = "ascii")


##########43pop_future_climate_GO
library(dplyr)
library(fuzzyjoin)
# Read the latitude and longitude data file
df_pop <- read.csv("15pop_longlat.csv", header = TRUE)
# Read the GO value data file
df_go <- read.csv("RDA_proj_offset_12650.csv", header = TRUE)
merged_data <- df_pop %>%
     fuzzy_inner_join(df_go, 
                      by = c("longitude" = "x", "latitude" = "y"),
                      match_fun = list(`>=`, `<=`))  # Use greater than or equal and less than or equal for fuzzy matching					  
merged_data <- merged_data %>%
     mutate(distance = sqrt((longitude - x)^2 + (latitude - y)^2))
write.csv(merged_data, "43_12650.csv", row.names = FALSE)
# Find the row with the minimum distance for each population and select the corresponding GO value
min_distance_rows <- merged_data %>%
  group_by(pop) %>%
  slice(which.min(distance)) %>%
  select(pop, Global_offset)
write.csv(min_distance_rows , "43_12650_min.csv", row.names = FALSE)


# Read the GO value data file
df_go <- read.csv("RDA_proj_offset_12690.csv", header = TRUE)
merged_data <- df_pop %>%
     fuzzy_inner_join(df_go, 
                      by = c("longitude" = "x", "latitude" = "y"),
                      match_fun = list(`>=`, `<=`))  # Use greater than or equal and less than or equal for fuzzy matching					  
merged_data <- merged_data %>%
     mutate(distance = sqrt((longitude - x)^2 + (latitude - y)^2))
write.csv(merged_data, "43_12690.csv", row.names = FALSE)
# Find the row with the minimum distance for each population and select the corresponding GO value
min_distance_rows <- merged_data %>%
  group_by(pop) %>%
  slice(which.min(distance)) %>%
  select(pop, Global_offset)
write.csv(min_distance_rows , "43_12690_min.csv", row.names = FALSE)


# Read the GO value data file
df_go <- read.csv("RDA_proj_offset_58550.csv", header = TRUE)
merged_data <- df_pop %>%
     fuzzy_inner_join(df_go, 
                      by = c("longitude" = "x", "latitude" = "y"),
                      match_fun = list(`>=`, `<=`))  # Use greater than or equal and less than or equal for fuzzy matching					  
merged_data <- merged_data %>%
     mutate(distance = sqrt((longitude - x)^2 + (latitude - y)^2))
write.csv(merged_data, "43_58550.csv", row.names = FALSE)
# Find the row with the minimum distance for each population and select the corresponding GO value
min_distance_rows <- merged_data %>%
  group_by(pop) %>%
  slice(which.min(distance)) %>%
  select(pop, Global_offset)
write.csv(min_distance_rows , "43_58550_min.csv", row.names = FALSE)


# Read the GO value data file
df_go <- read.csv("RDA_proj_offset_58590.csv", header = TRUE)
merged_data <- df_pop %>%
     fuzzy_inner_join(df_go, 
                      by = c("longitude" = "x", "latitude" = "y"),
                      match_fun = list(`>=`, `<=`))  # Use greater than or equal and less than or equal for fuzzy matching					  
merged_data <- merged_data %>%
     mutate(distance = sqrt((longitude - x)^2 + (latitude - y)^2))
write.csv(merged_data, "43_58590.csv", row.names = FALSE)
# Find the row with the minimum distance for each population and select the corresponding GO value
min_distance_rows <- merged_data %>%
  group_by(pop) %>%
  slice(which.min(distance)) %>%
  select(pop, Global_offset)
write.csv(min_distance_rows , "43_58590_min.csv", row.names = FALSE)

library(sp)
library(raster)
DH5026_raster <- raster("RDA_proj_offset_12650.asc")
DH9026_raster <- raster("RDA_proj_offset_12690.asc")
DH5085_raster <- raster("RDA_proj_offset_58550.asc")
DH9085_raster <- raster("RDA_proj_offset_58590.asc")

pdf(file = "corerdaDH_all_GO.pdf")
par(mfrow = c(2, 2), mar = c(2, 2, 2, 1), oma = c(1, 1, 1, 2))
# First subplot
plot(DH5026_raster, zlim = c(0.004, 2.1), legend = FALSE)  # Plot DH5026_raster, set z-axis range, hide legend
points_data <- read.csv("corerda_go.csv")  # Read GO values from CSV file
library(sp)  # Load spatial data processing library
# Create spatial points dataframe
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object
# Set coordinate reference system (assuming WGS84 for longitude and latitude)
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
# Calculate the maximum and minimum GO values
max_GO <- max(points_data$GO_5026)  # Get the maximum GO value
min_GO <- min(points_data$GO_5026)  # Get the minimum GO value
# Calculate symbol radius
points_data$radius <- 0.05 + 0.3 * (points_data$GO_5026 - min_GO) / (max_GO - min_GO)
# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(a)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label
# Second subplot
plot(DH9026_raster, zlim = c(0.004, 2.1), legend = FALSE)  # Plot DH9026_raster
points_data <- read.csv("corerda_go.csv")  # Read GO values from CSV file
# Create spatial points dataframe
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")  # Set coordinate reference system
# Calculate the maximum and minimum GO values
max_GO <- max(points_data$GO_9026)
min_GO <- min(points_data$GO_9026)
points_data$radius <- 0.05 + 0.3 * (points_data$GO_9026 - min_GO) / (max_GO - min_GO)
# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(b)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label
# Third subplot
plot(DH5085_raster, zlim = c(0.004, 2.1), legend = FALSE)  # Plot DH5085_raster
points_data <- read.csv("corerda_go.csv")  # Read GO values from CSV file
# Create spatial points dataframe
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")  # Set coordinate reference system
# Calculate the maximum and minimum GO values
max_GO <- max(points_data$GO_5085)
min_GO <- min(points_data$GO_5085)
points_data$radius <- 0.05 + 0.3 * (points_data$GO_5085 - min_GO) / (max_GO - min_GO)
# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(c)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label
# Fourth subplot
plot(DH9085_raster, zlim = c(0.004, 2.1))  # Plot DH9085_raster and display the legend
points_data <- read.csv("corerda_go.csv")  # Read GO values from CSV file
# Create spatial points dataframe
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")  # Set coordinate reference system
# Calculate the maximum and minimum GO values
max_GO <- max(points_data$GO_9085)
min_GO <- min(points_data$GO_9085)
points_data$radius <- 0.05 + 0.3 * (points_data$GO_9085 - min_GO) / (max_GO - min_GO)
# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(d)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label
dev.off()





###########################################################################
####pan east
AllFreq <- read.table("pan_east_AF.csv", sep = ",", header = T, row.names=1)
## Environment Variable Processing
# Historical Environment
library(raster)

# Get file list
file_list <- list.files("F:\\王聪颖\\景观_feng\\RDA_GO\\RDA_GO_171_3144\\历史\\适生区_asc\\east", pattern = ".tif$", full.names = TRUE)

# Load and check each file
raster_list <- lapply(file_list, raster)

# Check the attributes of each raster
for (i in seq_along(raster_list)) {
  cat("File", file_list[i], "has", nrow(raster_list[[i]]), "rows and", ncol(raster_list[[i]]), "columns.\n")
}

# Attempt to stack the files
tryCatch({
  ras <- stack(raster_list)
}, error = function(e) {
  message("Error:", e)
})

# Function to remove NAs from raster stack
remove.NAs.stack <- function(rast.stack) {
  nom <- names(rast.stack)
  test1 <- calc(rast.stack, fun = sum)
  test1[!is.na(test1)] <- 1
  test2 <- rast.stack * test1
  test2 <- stack(test2)
  names(test2) <- nom
  return(test2)
}

# Apply the remove.NAs.stack function to the raster stack created earlier
ras <- remove.NAs.stack(ras)

### 585
### 585_50
# Get file list
file_list <- list.files("F:\\王聪颖\\景观_feng\\RDA_GO\\RDA_GO_171_3144\\ssp585_50\\58550_适生区_asc\\east", pattern = ".tif$", full.names = TRUE)

# Load and check each file
raster_list <- lapply(file_list, raster)

# Check the attributes of each raster
for (i in seq_along(raster_list)) {
  cat("File", file_list[i], "has", nrow(raster_list[[i]]), "rows and", ncol(raster_list[[i]]), "columns.\n")
}

# Attempt to stack the files
tryCatch({
  ras_58550 <- stack(raster_list)
}, error = function(e) {
  message("Error:", e)
})

names(ras_58550)

# Apply the remove.NAs.stack function to the new raster stack
ras_58550 <- remove.NAs.stack(ras_58550)

#### 58590
# Get file list
file_list <- list.files("F:\\王聪颖\\景观_feng\\RDA_GO\\RDA_GO_171_3144\\ssp585_90\\58590_适生区_asc\\east", pattern = ".tif$", full.names = TRUE)

# Load and check each file
raster_list <- lapply(file_list, raster)

# Check the attributes of each raster
for (i in seq_along(raster_list)) {
cat("File", file_list[i], "has", nrow(raster_list[[i]]), "rows and", ncol(raster_list[[i]]), "columns.\n")
}
# Attempt to stack files
tryCatch({
  ras_58590 <- stack(raster_list)
}, error = function(e) {
  message("Error:", e)
})

# Get the names of the raster stack
names(ras_58590)

# Apply the remove.NAs.stack function to handle missing values
ras_58590 <- remove.NAs.stack(ras_58590)

# Read the coordinates data
Coordinates <- read.table("15pop_longlat.csv", sep = ",", header = TRUE, row.names = 1)

# Match coordinates data with row names of AllFreq
Coordinates <- Coordinates[match(row.names(AllFreq), Coordinates$id2, nomatch = 0),]

# Rename the columns of the coordinates data
colnames(Coordinates) <- c("Population", "Longitude", "Latitude", "Elevation")

# Extract environmental variables corresponding to the coordinates from the raster 'ras' and create a data frame
Env <- data.frame(extract(ras, Coordinates[, 2:3]))

# Standardize the variables
Env <- scale(Env, center = TRUE, scale = TRUE) # center=TRUE and scale=TRUE are the default settings of the scale() function

# Recover the scaling factors
scale_env <- attr(Env, 'scaled:scale') # Get the scaling factors
center_env <- attr(Env, 'scaled:center') # Get the center values

# Convert the standardized environmental data into a data frame
Env <- as.data.frame(Env)

# Set row names of the data frame as population names
row.names(Env) <- c(Coordinates$Population) 

# Import other R scripts (not specified in the provided code)
source("F:\\王聪颖\\景观_feng\\RDA_GO\\RDA_GO_171_3144\\输入脚本\\genomic_offset.R")
source("F:\\王聪颖\\景观_feng\\RDA_GO\\RDA_GO_171_3144\\输入脚本\\provgar_offset.R")
# Perform Principal Component Analysis (PCA)
pca <- rda(AllFreq, scale = TRUE)

# Plot the Eigenvalue chart for PCA
screeplot(pca, type = "barplot", npcs = 10, main = "PCA Eigenvalues")

### Extract the first few principal components (PCs)
PCs <- scores(pca, choices = c(1:3), display = "sites", scaling = 0)

# Create a data frame containing populations and principal components
PopStruct <- data.frame(POP = row.names(AllFreq), PCs)
colnames(PopStruct) <- c("POP", "PC1", "PC2", "PC3") # Rename columns
write.csv(PopStruct, file = "PopStruct.csv", quote = FALSE)

# Merge coordinates, population structure, and environmental variables into one data frame
Variables <- data.frame(Coordinates, PopStruct[,-1], Env)

# Perform Redundancy Analysis (RDA)
RDA_outliers <- rda(AllFreq ~ BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15, Variables)

## Run the genomic_offset function for specific populations (58550 and 58590)
res_RDA_proj58590 <- genomic_offset(RDA_outliers, K = 3, env_pres = ras, env_fut = ras_58590, range = NULL, method = "loadings", scale_env = scale_env, center_env = center_env)
res_RDA_proj58550 <- genomic_offset(RDA_outliers, K = 3, env_pres = ras, env_fut = ras_58550, range = NULL, method = "loadings", scale_env = scale_env, center_env = center_env)
## Table global genetic offset predicted for 58550 and 58590
RDA_proj_offset <- data.frame(rbind(rasterToPoints(res_RDA_proj58550$Proj_offset_global), rasterToPoints(res_RDA_proj58590$Proj_offset_global)), Date = c(rep("58550", nrow(rasterToPoints(res_RDA_proj58550$Proj_offset_global))), rep("58590", nrow(rasterToPoints(res_RDA_proj58590$Proj_offset_global)))))
write.csv(RDA_proj_offset,file="RDA_proj_offset.csv",quote=F)



RDA_proj_offset_58550 <- read.csv("RDA_proj_offset_58550.csv")
raster_data <- rasterFromXYZ(RDA_proj_offset_58550[, c("x", "y", "Global_offset")])
asc_file <- "RDA_proj_offset_58550.asc"
writeRaster(raster_data, filename = asc_file, format = "ascii")
 

RDA_proj_offset_58590 <- read.csv("RDA_proj_offset_58590.csv")
raster_data <- rasterFromXYZ(RDA_proj_offset_58590[, c("x", "y", "Global_offset")])
asc_file <- "RDA_proj_offset_58590.asc"
writeRaster(raster_data, filename = asc_file, format = "ascii")



####126
## Climate data for 12690
# Get the list of files
file_list <- list.files("F:\\王聪颖\\景观_feng\\RDA_GO\\RDA_GO_171_3144\\ssp126_90\\12690_适生区_asc\\east", pattern = ".tif$", full.names = TRUE)
# Load and check each file
raster_list <- lapply(file_list, raster)
# Check the attributes of each raster
for (i in seq_along(raster_list)) {
  cat("File", file_list[i], "has", nrow(raster_list[[i]]), "rows and", ncol(raster_list[[i]]), "columns.\n")
}
# Attempt to stack the files
tryCatch({
  ras_12690 <- stack(raster_list)
}, error = function(e) {
  message("Error:", e)
})
names(ras_12690)
ras_12690 <- remove.NAs.stack(ras_12690)

####12650
## Climate data for 12650
# Get the list of files
file_list <- list.files("F:\\王聪颖\\景观feng\\RDA_GO\\RDA_GO_171_3144\\ssp126_50\\12650_适生区_asc\\east", pattern = ".tif$", full.names = TRUE)
# Load and check each file
raster_list <- lapply(file_list, raster)
# Check the attributes of each raster
for (i in seq_along(raster_list)) {
  cat("File", file_list[i], "has", nrow(raster_list[[i]]), "rows and", ncol(raster_list[[i]]), "columns.\n")
}
# Attempt to stack the files
tryCatch({
  ras_12650 <- stack(raster_list)
}, error = function(e) {
  message("Error:", e)
})
names(ras_12650)
ras_12650 <- remove.NAs.stack(ras_12650)

# Perform genomic offset projections for future climate scenarios
res_RDA_proj12690 <- genomic_offset(RDA_outliers, K = 3, env_pres = ras, env_fut = ras_12690, range = NULL, method = "loadings", scale_env = scale_env, center_env = center_env)
res_RDA_proj12650 <- genomic_offset(RDA_outliers, K = 3, env_pres = ras, env_fut = ras_12650, range = NULL, method = "loadings", scale_env = scale_env, center_env = center_env)

# Combine results and write to CSV
RDA_proj_offset <- data.frame(rbind(rasterToPoints(res_RDA_proj12650$Proj_offset_global), rasterToPoints(res_RDA_proj12690$Proj_offset_global)),
                              Date = c(rep("12650", nrow(rasterToPoints(res_RDA_proj12650$Proj_offset_global))), 
                                       rep("12690", nrow(rasterToPoints(res_RDA_proj12690$Proj_offset_global)))))
write.csv(RDA_proj_offset, "RDA_proj_offset_12650_90.csv", row.names = FALSE)

# Convert results to raster and save as ASC file
RDA_proj_offset_12650 <- read.csv("RDA_proj_offset_12650.csv")
raster_data <- rasterFromXYZ(RDA_proj_offset_12650[, c("x", "y", "Global_offset")])
asc_file <- "RDA_proj_offset_12650.asc"
writeRaster(raster_data, filename = asc_file, format = "ascii")

RDA_proj_offset_12690 <- read.csv("RDA_proj_offset_12690.csv")
raster_data <- rasterFromXYZ(RDA_proj_offset_12690[, c("x", "y", "Global_offset")])
asc_file <- "RDA_proj_offset_12690.asc"
writeRaster(raster_data, filename = asc_file, format = "ascii")


########## 15pop_ Future Climate GO
library(dplyr)
library(fuzzyjoin)

# Read latitude and longitude data file
df_pop <- read.csv("15poplonglat.csv", header = TRUE)

# Read GO values data file
df_go <- read.csv("RDA_proj_offset_12650.csv", header = TRUE)

# Perform fuzzy join based on latitude and longitude
merged_data <- df_pop %>%
  fuzzy_inner_join(df_go, 
                   by = c("longitude" = "x", "latitude" = "y"),
                   match_fun = list(`>=`, `<=`))  # Use fuzzy matching with >= and <=

# Calculate distance between each population point and the raster points
merged_data <- merged_data %>%
  mutate(distance = sqrt((longitude - x)^2 + (latitude - y)^2))

# Write merged data to CSV
write.csv(merged_data, "15_12650.csv", row.names = FALSE)

# Find the row with the minimum distance for each population and select the corresponding GO value
min_distance_rows <- merged_data %>%
  group_by(pop) %>%
  slice(which.min(distance)) %>%
  select(pop, Global_offset)

# Write minimum distance rows to CSV
write.csv(min_distance_rows, "15_12650_min.csv", row.names = FALSE)


# Read GO values data file for 12690
df_go <- read.csv("RDA_proj_offset_12690.csv", header = TRUE)

# Perform fuzzy join based on latitude and longitude
merged_data <- df_pop %>%
  fuzzy_inner_join(df_go, 
                   by = c("longitude" = "x", "latitude" = "y"),
                   match_fun = list(`>=`, `<=`))  # Use fuzzy matching with >= and <=

# Calculate distance between each population point and the raster points
merged_data <- merged_data %>%
  mutate(distance = sqrt((longitude - x)^2 + (latitude - y)^2))

# Write merged data to CSV
write.csv(merged_data, "15_12690.csv", row.names = FALSE)

# Find the row with the minimum distance for each population and select the corresponding GO value
min_distance_rows <- merged_data %>%
  group_by(pop) %>%
  slice(which.min(distance)) %>%
  select(pop, Global_offset)

# Write minimum distance rows to CSV
write.csv(min_distance_rows, "15_12690_min.csv", row.names = FALSE)




# Read GO values data file
df_go <- read.csv("RDA_proj_offset_58550.csv", header = TRUE)

# Perform fuzzy join based on longitude and latitude
merged_data <- df_pop %>%
     fuzzy_inner_join(df_go, 
                      by = c("longitude" = "x", "latitude" = "y"),
                      match_fun = list(`>=`, `<=`))  # Use fuzzy matching with >= and <=

# Calculate the distance between each population point and the raster points
merged_data <- merged_data %>%
     mutate(distance = sqrt((longitude - x)^2 + (latitude - y)^2))

# Write the merged data to a CSV file
write.csv(merged_data, "15_58550.csv", row.names = FALSE)

# Find the row with the minimum distance for each population and select the corresponding GO value
min_distance_rows <- merged_data %>%
  group_by(pop) %>%
  slice(which.min(distance)) %>%
  select(pop, Global_offset)

# Write the minimum distance rows to a CSV file
write.csv(min_distance_rows , "15_58550_min.csv", row.names = FALSE)


# Read GO values data file
df_go <- read.csv("RDA_proj_offset_58590.csv", header = TRUE)

# Perform fuzzy join based on longitude and latitude
merged_data <- df_pop %>%
     fuzzy_inner_join(df_go, 
                      by = c("longitude" = "x", "latitude" = "y"),
                      match_fun = list(`>=`, `<=`))  # Use fuzzy matching with >= and <=

# Calculate the distance between each population point and the raster points
merged_data <- merged_data %>%
     mutate(distance = sqrt((longitude - x)^2 + (latitude - y)^2))

# Write the merged data to a CSV file
write.csv(merged_data, "15_58590.csv", row.names = FALSE)

# Find the row with the minimum distance for each population and select the corresponding GO value
min_distance_rows <- merged_data %>%
  group_by(pop) %>%
  slice(which.min(distance)) %>%
  select(pop, Global_offset)

# Write the minimum distance rows to a CSV file
write.csv(min_distance_rows , "15_58590_min.csv", row.names = FALSE)



library(sp)
library(raster)
DH5026_raster <- raster("RDA_proj_offset_12650.asc")
DH9026_raster <- raster("RDA_proj_offset_12690.asc")
DH5085_raster <- raster("RDA_proj_offset_58550.asc")
DH9085_raster <- raster("RDA_proj_offset_58590.asc")

pdf(file = "paneastrdaDH_all_GO.pdf")
par(mfrow = c(2, 2), mar = c(2, 2, 2, 1), oma = c(1, 1, 1, 2))

# First subplot
plot(DH5026_raster, zlim = c(0.009, 1), legend = FALSE)  # Plot DH5026_raster with z-axis limits and hide legend
points_data <- read.csv("paneast_go.csv")  # Read GO values data from CSV
library(sp)  # Load spatial data processing library
# Create spatial points data frame
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object
# Set coordinate reference system (assuming longitude and latitude use WGS84)
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
# Calculate max and min GO values
max_GO <- max(points_data$GO_5026)  # Get the maximum GO value
min_GO <- min(points_data$GO_5026)  # Get the minimum GO value
# Calculate symbol radius
points_data$radius <- 0.05 + 0.3 * (points_data$GO_5026 - min_GO) / (max_GO - min_GO)
# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(a)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label

# Second subplot
plot(DH9026_raster, zlim = c(0.009, 1), legend = FALSE)  # Plot DH9026_raster
points_data <- read.csv("paneast_go.csv")  # Read GO values data from CSV
# Create spatial points data frame
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")  # Set coordinate reference system
# Calculate max and min GO values
max_GO <- max(points_data$GO_9026)
min_GO <- min(points_data$GO_9026)
points_data$radius <- 0.05 + 0.3 * (points_data$GO_9026 - min_GO) / (max_GO - min_GO)
# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(b)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label

# Third subplot
plot(DH5085_raster, zlim = c(0.009, 1), legend = FALSE)  # Plot DH5085_raster
points_data <- read.csv("paneast_go.csv")  # Read GO values data from CSV
# Create spatial points data frame
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")  # Set coordinate reference system
# Calculate max and min GO values
max_GO <- max(points_data$GO_5085)
min_GO <- min(points_data$GO_5085)
points_data$radius <- 0.05 + 0.3 * (points_data$GO_5085 - min_GO) / (max_GO - min_GO)
# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(c)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label

# Fourth subplot
plot(DH9085_raster, zlim = c(0.009, 1))  # Plot DH9085_raster with legend
points_data <- read.csv("paneast_go.csv")  # Read GO values data from CSV
# Create spatial points data frame
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")  # Set coordinate reference system
# Calculate max and min GO values
max_GO <- max(points_data$GO_9085)
min_GO <- min(points_data$GO_9085)
points_data$radius <- 0.05 + 0.3 * (points_data$GO_9085 - min_GO) / (max_GO - min_GO)
# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(d)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label

dev.off()



################################
## Core East
AllFreq <- read.table("core_east_AF.csv", sep = ",", header = TRUE, row.names = 1)

# Perform Principal Component Analysis (PCA)
pca <- rda(AllFreq, scale = TRUE)

# Plot the PCA eigenvalue graph
screeplot(pca, type = "barplot", npcs = 10, main = "PCA Eigenvalues")

### Extract the first few principal components (PCs)
PCs <- scores(pca, choices = c(1:3), display = "sites", scaling = 0)

# Create a data frame containing populations and principal components
PopStruct <- data.frame(POP = row.names(AllFreq), PCs)
colnames(PopStruct) <- c("POP", "PC1", "PC2", "PC3") # Rename columns
write.csv(PopStruct, file = "PopStruct.csv", quote = FALSE)

# Merge coordinates, population structure, and environmental variables into one data frame
Variables <- data.frame(Coordinates, PopStruct[,-1], Env)

# Perform Redundancy Analysis (RDA)
RDA_outliers <- rda(AllFreq ~ BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15, Variables)

## Run the genomic_offset function for specified populations (58550 and 58590)
res_RDA_proj58590 <- genomic_offset(RDA_outliers, K = 3, env_pres = ras, env_fut = ras_58590, range = NULL, method = "loadings", scale_env = scale_env, center_env = center_env)
res_RDA_proj58550 <- genomic_offset(RDA_outliers, K = 3, env_pres = ras, env_fut = ras_58550, range = NULL, method = "loadings", scale_env = scale_env, center_env = center_env)
## Table global genetic offset predicted for 58550 and 58590
RDA_proj_offset <- data.frame(rbind(rasterToPoints(res_RDA_proj58550$Proj_offset_global), rasterToPoints(res_RDA_proj58590$Proj_offset_global)), Date = c(rep("58550", nrow(rasterToPoints(res_RDA_proj58550$Proj_offset_global))), rep("58590", nrow(rasterToPoints(res_RDA_proj58590$Proj_offset_global)))))
write.csv(RDA_proj_offset,file="RDA_proj_offset.csv",quote=F)



RDA_proj_offset_58550 <- read.csv("RDA_proj_offset_58550.csv")
raster_data <- rasterFromXYZ(RDA_proj_offset_58550[, c("x", "y", "Global_offset")])
asc_file <- "RDA_proj_offset_58550.asc"
writeRaster(raster_data, filename = asc_file, format = "ascii")
 

RDA_proj_offset_58590 <- read.csv("RDA_proj_offset_58590.csv")
raster_data <- rasterFromXYZ(RDA_proj_offset_58590[, c("x", "y", "Global_offset")])
asc_file <- "RDA_proj_offset_58590.asc"
writeRaster(raster_data, filename = asc_file, format = "ascii")

res_RDA_proj12690 <- genomic_offset(RDA_outliers, K = 3, env_pres = ras, env_fut = ras_12690, range = NULL, method = "loadings", scale_env = scale_env, center_env = center_env)
res_RDA_proj12650 <- genomic_offset(RDA_outliers, K = 3, env_pres = ras, env_fut = ras_12650, range = NULL, method = "loadings", scale_env = scale_env, center_env = center_env)
RDA_proj_offset <- data.frame(rbind(rasterToPoints(res_RDA_proj12650$Proj_offset_global), rasterToPoints(res_RDA_proj12690$Proj_offset_global)), Date = c(rep("12650", nrow(rasterToPoints(res_RDA_proj12650$Proj_offset_global))), rep("12690", nrow(rasterToPoints(res_RDA_proj12690$Proj_offset_global)))))
write.csv(RDA_proj_offset, "RDA_proj_offset_12650_90.csv", row.names = FALSE)

RDA_proj_offset_12650 <- read.csv("RDA_proj_offset_12650.csv")
raster_data <- rasterFromXYZ(RDA_proj_offset_12650[, c("x", "y", "Global_offset")])
asc_file <- "RDA_proj_offset_12650.asc"
writeRaster(raster_data, filename = asc_file, format = "ascii")
 
RDA_proj_offset_12690 <- read.csv("RDA_proj_offset_12690.csv")
raster_data <- rasterFromXYZ(RDA_proj_offset_12690[, c("x", "y", "Global_offset")])
asc_file <- "RDA_proj_offset_12690.asc"
writeRaster(raster_data, filename = asc_file, format = "ascii")


########## 15pop Future Climate GO
library(dplyr)
library(fuzzyjoin)

# Read the latitude and longitude data file
df_pop <- read.csv("15poplonglat.csv", header = TRUE)

# Read the GO value data file
df_go <- read.csv("RDA_proj_offset_12650.csv", header = TRUE)
merged_data <- df_pop %>%
     fuzzy_inner_join(df_go, 
                      by = c("longitude" = "x", "latitude" = "y"),
                      match_fun = list(`>=`, `<=`))  # Use greater than or equal to and less than or equal to for fuzzy matching
merged_data <- merged_data %>%
     mutate(distance = sqrt((longitude - x)^2 + (latitude - y)^2))
write.csv(merged_data, "15_12650.csv", row.names = FALSE)

# Find the row with the smallest distance for each population and select the corresponding GO value
min_distance_rows <- merged_data %>%
  group_by(pop) %>%
  slice(which.min(distance)) %>%
  select(pop, Global_offset)
write.csv(min_distance_rows , "15_12650_min.csv", row.names = FALSE)


# Read the GO value data file
df_go <- read.csv("RDA_proj_offset_12690.csv", header = TRUE)
merged_data <- df_pop %>%
     fuzzy_inner_join(df_go, 
                      by = c("longitude" = "x", "latitude" = "y"),
                      match_fun = list(`>=`, `<=`))  # Use greater than or equal to and less than or equal to for fuzzy matching
merged_data <- merged_data %>%
     mutate(distance = sqrt((longitude - x)^2 + (latitude - y)^2))
write.csv(merged_data, "15_12690.csv", row.names = FALSE)

# Find the row with the smallest distance for each population and select the corresponding GO value
min_distance_rows <- merged_data %>%
  group_by(pop) %>%
  slice(which.min(distance)) %>%
  select(pop, Global_offset)
write.csv(min_distance_rows , "15_12690_min.csv", row.names = FALSE)


# Read the GO value data file
df_go <- read.csv("RDA_proj_offset_58550.csv", header = TRUE)
merged_data <- df_pop %>%
     fuzzy_inner_join(df_go, 
                      by = c("longitude" = "x", "latitude" = "y"),
                      match_fun = list(`>=`, `<=`))  # Use greater than or equal to and less than or equal to for fuzzy matching
merged_data <- merged_data %>%
     mutate(distance = sqrt((longitude - x)^2 + (latitude - y)^2))
write.csv(merged_data, "15_58550.csv", row.names = FALSE)

# Find the row with the smallest distance for each population and select the corresponding GO value
min_distance_rows <- merged_data %>%
  group_by(pop) %>%
  slice(which.min(distance)) %>%
  select(pop, Global_offset)
write.csv(min_distance_rows , "15_58550_min.csv", row.names = FALSE)


# Read the GO value data file
df_go <- read.csv("RDA_proj_offset_58590.csv", header = TRUE)
merged_data <- df_pop %>%
     fuzzy_inner_join(df_go, 
                      by = c("longitude" = "x", "latitude" = "y"),
                      match_fun = list(`>=`, `<=`))  # Use greater than or equal to and less than or equal to for fuzzy matching
merged_data <- merged_data %>%
     mutate(distance = sqrt((longitude - x)^2 + (latitude - y)^2))
write.csv(merged_data, "15_58590.csv", row.names = FALSE)

# Find the row with the smallest distance for each population and select the corresponding GO value
min_distance_rows <- merged_data %>%
  group_by(pop) %>%
  slice(which.min(distance)) %>%
  select(pop, Global_offset)
write.csv(min_distance_rows , "15_58590_min.csv", row.names = FALSE)


library(sp)
library(raster)
DH5026_raster <- raster("RDA_proj_offset_12650.asc")
DH9026_raster <- raster("RDA_proj_offset_12690.asc")
DH5085_raster <- raster("RDA_proj_offset_58550.asc")
DH9085_raster <- raster("RDA_proj_offset_58590.asc")

pdf(file = "paneastrdaDH_all_GO.pdf")
par(mfrow = c(2, 2), mar = c(2, 2, 2, 1), oma = c(1, 1, 1, 2))
# First subplot
plot(DH5026_raster, zlim = c(0.01, 1), legend = FALSE)  # Plot DH5026_raster, set z-axis range, hide legend
points_data <- read.csv("paneast_go.csv")  # Read GO values data from CSV file
library(sp)  # Load spatial data processing library
# Create spatial points data frame
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object
# Set coordinate reference system (assuming latitude and longitude are using WGS84)
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
# Calculate the maximum and minimum GO values
max_GO <- max(points_data$GO_5026)  # Get maximum GO value
min_GO <- min(points_data$GO_5026)  # Get minimum GO value
# Calculate symbol radius
points_data$radius <- 0.05 + 0.3 * (points_data$GO_5026 - min_GO) / (max_GO - min_GO)
# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(a)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label
# Second subplot
plot(DH9026_raster, zlim = c(0.01, 1), legend = FALSE)  # Plot DH9026_raster
points_data <- read.csv("paneast_go.csv")  # Read GO values data from CSV file
# Create spatial points data frame
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")  # Set coordinate reference system
# Calculate the maximum and minimum GO values
max_GO <- max(points_data$GO_9026)
min_GO <- min(points_data$GO_9026)
points_data$radius <- 0.05 + 0.3 * (points_data$GO_9026 - min_GO) / (max_GO - min_GO)
# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(b)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label
# Third subplot
plot(DH5085_raster, zlim = c(0.01, 1), legend = FALSE)  # Plot DH5085_raster
points_data <- read.csv("paneast_go.csv")  # Read GO values data from CSV file
# Create spatial points data frame
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")  # Set coordinate reference system
# Calculate the maximum and minimum GO values
max_GO <- max(points_data$GO_5085)
min_GO <- min(points_data$GO_5085)
points_data$radius <- 0.05 + 0.3 * (points_data$GO_5085 - min_GO) / (max_GO - min_GO)
# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(c)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label
# Fourth subplot
plot(DH9085_raster, zlim = c(0.01, 1))  # Plot DH9085_raster and show legend
points_data <- read.csv("paneast_go.csv")  # Read GO values data from CSV file
# Create spatial points data frame
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")  # Set coordinate reference system
# Calculate the maximum and minimum GO values
max_GO <- max(points_data$GO_9085)
min_GO <- min(points_data$GO_9085)
points_data$radius <- 0.05 + 0.3 * (points_data$GO_9085 - min_GO) / (max_GO - min_GO)
# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(d)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label
dev.off()









####pan west
AllFreq <- read.table("pan_west_AF.csv", sep = ",", header = T, row.names=1)

## Environmental variable processing
# Historical Environment
library(raster)

# Get the list of files
file_list <- list.files("F:\\王聪颖\\景观_feng\\RDA_GO\\RDA_GO_171_3144\\历史\\适生区 asc\\west", pattern = ".tif$", full.names = TRUE)

# Load and check each file
raster_list <- lapply(file_list, raster)

# Check the attributes of each raster
for (i in seq_along(raster_list)) {
  cat("File", file_list[i], "has", nrow(raster_list[[i]]), "rows and", ncol(raster_list[[i]]), "columns.\n")
}

# Try stacking the files
tryCatch({
  ras <- stack(raster_list)
}, error = function(e) {
  message("Error:", e)
})

# Function to remove NAs from the stack
remove.NAs.stack <- function(rast.stack) {
  nom <- names(rast.stack)
  test1 <- calc(rast.stack, fun = sum)
  test1[!is.na(test1)] <- 1
  test2 <- rast.stack * test1
  test2 <- stack(test2)
  names(test2) <- nom
  return(test2)
}

ras <- remove.NAs.stack(ras) # Apply remove.NAs.stack function to the previously created raster stack

### 585
### 585_50
# Get the list of files
file_list <- list.files("F:\\王聪颖\\景观_feng\\RDA_GO\\RDA_GO_171_3144\\ssp585_50\\58550适生区asc\\west", pattern = ".tif$", full.names = TRUE)

# Load and check each file
raster_list <- lapply(file_list, raster)

# Check the attributes of each raster
for (i in seq_along(raster_list)) {
  cat("File", file_list[i], "has", nrow(raster_list[[i]]), "rows and", ncol(raster_list[[i]]), "columns.\n")
}

# Try stacking the files
tryCatch({
  ras_58550 <- stack(raster_list)
}, error = function(e) {
  message("Error:", e)
})

names(ras_58550)
ras_58550 <- remove.NAs.stack(ras_58550)

#### 58590
# Get the list of files
file_list <- list.files("F:\\王聪颖\\景观_feng\\RDA_GO\\RDA_GO_171_3144\\ssp585_90\\58590适生区asc\\west", pattern = ".tif$", full.names = TRUE)

# Load and check each file
raster_list <- lapply(file_list, raster)

# Check the attributes of each raster
for (i in seq_along(raster_list)) {
  cat("File", file_list[i], "has", nrow(raster_list[[i]]), "rows and", ncol(raster_list[[i]]), "columns.\n")
}

# Try stacking the files
tryCatch({
ras_58590 <- stack(raster_list)
}, error = function(e) {
message("Error:", e)
})
names(ras_58590)
ras_58590 <- remove.NAs.stack(ras_58590)

Coordinates <- read.table("28pop_longlat.csv", sep = ",", header = TRUE, row.names = 1)
# Match the coordinate data based on the row names of AllFreq
Coordinates <- Coordinates[match(row.names(AllFreq), Coordinates$id2, nomatch = 0),]
# Rename the columns of the coordinate data
colnames(Coordinates) <- c("Population", "Longitude", "Latitude", "Elevation")
# Extract environmental variables from the raster 'ras' based on the coordinates and create a data frame
Env <- data.frame(extract(ras, Coordinates[, 2:3]))
## Standardize the variables
Env <- scale(Env, center = TRUE, scale = TRUE) # center=TRUE and scale=TRUE are the default settings for the scale() function
## Retrieve the scaling coefficients
scale_env <- attr(Env, 'scaled:scale') # Get scaling factors
center_env <- attr(Env, 'scaled:center') # Get centering values
## Convert the standardized environmental data into a data frame
Env <- as.data.frame(Env)
row.names(Env) <- c(Coordinates$Population) # Set row names to the population names

# Import additional R scripts
source("F:\\王聪颖\\景观\\RDA_GO\\RDA_GO_171_3144\\输入脚本\\genomic_offset.R")
source("F:\\王聪颖\\景观_feng\\RDA_GO\\RDA_GO_171_3144\\输入脚本\\provgar_offset.R")

# Perform Principal Component Analysis (PCA)
pca <- rda(AllFreq, scale = TRUE)

# Plot the PCA eigenvalue bar plot
screeplot(pca, type = "barplot", npcs = 10, main = "PCA Eigenvalues")

### Extract the first few principal components (PCs)
PCs <- scores(pca, choices = c(1:3), display = "sites", scaling = 0)

# Create a data frame with populations and the principal components
PopStruct <- data.frame(POP = row.names(AllFreq), PCs)
colnames(PopStruct) <- c("POP", "PC1", "PC2", "PC3") # Rename columns
write.csv(PopStruct, file = "PopStruct.csv", quote = FALSE)

# Combine coordinates, population structure, and environmental variables into one data frame
Variables <- data.frame(Coordinates, PopStruct[,-1], Env)

# Perform Redundancy Analysis (RDA)
RDA_outliers <- rda(AllFreq ~ BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15, Variables)

## Run the genomic_offset function for specified populations (58550 and 58590)
res_RDA_proj58590 <- genomic_offset(RDA_outliers, K = 3, env_pres = ras, env_fut = ras_58590, range = NULL, method = "loadings", scale_env = scale_env, center_env = center_env)
res_RDA_proj58550 <- genomic_offset(RDA_outliers, K = 3, env_pres = ras, env_fut = ras_58550, range = NULL, method = "loadings", scale_env = scale_env, center_env = center_env)
## Table global genetic offset predicted for 58550 and 58590
RDA_proj_offset <- data.frame(rbind(rasterToPoints(res_RDA_proj58550$Proj_offset_global), rasterToPoints(res_RDA_proj58590$Proj_offset_global)), Date = c(rep("58550", nrow(rasterToPoints(res_RDA_proj58550$Proj_offset_global))), rep("58590", nrow(rasterToPoints(res_RDA_proj58590$Proj_offset_global)))))
write.csv(RDA_proj_offset,file="RDA_proj_offset.csv",quote=F)



RDA_proj_offset_58550 <- read.csv("RDA_proj_offset_58550.csv")
raster_data <- rasterFromXYZ(RDA_proj_offset_58550[, c("x", "y", "Global_offset")])
asc_file <- "RDA_proj_offset_58550.asc"
writeRaster(raster_data, filename = asc_file, format = "ascii")
 

RDA_proj_offset_58590 <- read.csv("RDA_proj_offset_58590.csv")
raster_data <- rasterFromXYZ(RDA_proj_offset_58590[, c("x", "y", "Global_offset")])
asc_file <- "RDA_proj_offset_58590.asc"
writeRaster(raster_data, filename = asc_file, format = "ascii")



####126
## Climate data for 12690
# Get the file list
file_list <- list.files("F:\\王聪颖\\景观_feng\\RDA_GO\\RDA_GO_171_3144\\ssp126_90\\12690_适生区_asc\\west", pattern = ".tif$", full.names = TRUE)
# Load each file and check
raster_list <- lapply(file_list, raster)
# Check the attributes of each raster
for (i in seq_along(raster_list)) {
  cat("File", file_list[i], "has", nrow(raster_list[[i]]), "rows and", ncol(raster_list[[i]]), "columns.\n")
}
# Try stacking the files
tryCatch({
  ras_12690 <- stack(raster_list)
}, error = function(e) {
  message("Error:", e)
})
names(ras_12690)
ras_12690 <- remove.NAs.stack(ras_12690)

####12650
## Climate data for 12650
# Get the file list
file_list <- list.files("F:\\王聪颖\\景观_feng\\RDA_GO\\RDA_GO_171_3144\\ssp126_50\\12650_适生区_asc\\west", pattern = ".tif$", full.names = TRUE)
# Load each file and check
raster_list <- lapply(file_list, raster)
# Check the attributes of each raster
for (i in seq_along(raster_list)) {
  cat("File", file_list[i], "has", nrow(raster_list[[i]]), "rows and", ncol(raster_list[[i]]), "columns.\n")
}
# Try stacking the files
tryCatch({
  ras_12650 <- stack(raster_list)
}, error = function(e) {
  message("Error:", e)
})
names(ras_12650)
ras_12650 <- remove.NAs.stack(ras_12650)

# Perform genomic offset analysis for 12690 and 12650
res_RDA_proj12690 <- genomic_offset(RDA_outliers, K = 3, env_pres = ras, env_fut = ras_12690, range = NULL, method = "loadings", scale_env = scale_env, center_env = center_env)
res_RDA_proj12650 <- genomic_offset(RDA_outliers, K = 3, env_pres = ras, env_fut = ras_12650, range = NULL, method = "loadings", scale_env = scale_env, center_env = center_env)

# Combine results for 12650 and 12690 into one data frame
RDA_proj_offset <- data.frame(
  rbind(
    rasterToPoints(res_RDA_proj12650$Proj_offset_global),
    rasterToPoints(res_RDA_proj12690$Proj_offset_global)
  ),
  Date = c(rep("12650", nrow(rasterToPoints(res_RDA_proj12650$Proj_offset_global))),
           rep("12690", nrow(rasterToPoints(res_RDA_proj12690$Proj_offset_global))))
)

# Write the combined result to a CSV file
write.csv(RDA_proj_offset, "RDA_proj_offset_12650_90.csv", row.names = FALSE)

# Convert the results for 12650 to raster and save as ASC file
RDA_proj_offset_12650 <- read.csv("RDA_proj_offset_12650.csv")
raster_data <- rasterFromXYZ(RDA_proj_offset_12650[, c("x", "y", "Global_offset")])
asc_file <- "RDA_proj_offset_12650.asc"
writeRaster(raster_data, filename = asc_file, format = "ascii")

# Convert the results for 12690 to raster and save as ASC file
RDA_proj_offset_12690 <- read.csv("RDA_proj_offset_12690.csv")
raster_data <- rasterFromXYZ(RDA_proj_offset_12690[, c("x", "y", "Global_offset")])
asc_file <- "RDA_proj_offset_12690.asc"
writeRaster(raster_data, filename = asc_file, format = "ascii")


##########28pop_Future Climate GO
library(dplyr)
library(fuzzyjoin)

# Read population longitude and latitude data
df_pop <- read.csv("28poplonglat.csv", header = TRUE)

# Read the GO offset values for 12650
df_go <- read.csv("RDA_proj_offset_12650.csv", header = TRUE)
merged_data <- df_pop %>%
  fuzzy_inner_join(df_go, 
                   by = c("longitude" = "x", "latitude" = "y"),
                   match_fun = list(`>=`, `<=`))  # Use fuzzy matching with >= and <=
merged_data <- merged_data %>%
  mutate(distance = sqrt((longitude - x)^2 + (latitude - y)^2))  # Calculate distance
write.csv(merged_data, "28_12650.csv", row.names = FALSE)

# Find the row with the minimum distance for each population and select corresponding GO value
min_distance_rows <- merged_data %>%
  group_by(pop) %>%
  slice(which.min(distance)) %>%
  select(pop, Global_offset)
write.csv(min_distance_rows , "28_12650_min.csv", row.names = FALSE)


# Repeat the process for 12690
df_go <- read.csv("RDA_proj_offset_12690.csv", header = TRUE)
merged_data <- df_pop %>%
  fuzzy_inner_join(df_go, 
                   by = c("longitude" = "x", "latitude" = "y"),
                   match_fun = list(`>=`, `<=`))  # Use fuzzy matching with >= and <=
merged_data <- merged_data %>%
  mutate(distance = sqrt((longitude - x)^2 + (latitude - y)^2))  # Calculate distance
write.csv(merged_data, "28_12690.csv", row.names = FALSE)

# Find the row with the minimum distance for each population and select corresponding GO value
min_distance_rows <- merged_data %>%
  group_by(pop) %>%
  slice(which.min(distance)) %>%
  select(pop, Global_offset)
write.csv(min_distance_rows , "28_12690_min.csv", row.names = FALSE)


# Repeat the process for population 58550
df_go <- read.csv("RDA_proj_offset_58550.csv", header = TRUE)
merged_data <- df_pop %>%
  fuzzy_inner_join(df_go, 
                   by = c("longitude" = "x", "latitude" = "y"),
                   match_fun = list(`>=`, `<=`))  # Use fuzzy matching with >= and <=
merged_data <- merged_data %>%
  mutate(distance = sqrt((longitude - x)^2 + (latitude - y)^2))  # Calculate distance
write.csv(merged_data, "28_58550.csv", row.names = FALSE)

# Find the row with the minimum distance for each population and select corresponding GO value
min_distance_rows <- merged_data %>%
  group_by(pop) %>%
  slice(which.min(distance)) %>%
  select(pop, Global_offset)
write.csv(min_distance_rows , "28_58550_min.csv", row.names = FALSE)


# Repeat the process for population 58590
df_go <- read.csv("RDA_proj_offset_58590.csv", header = TRUE)
merged_data <- df_pop %>%
  fuzzy_inner_join(df_go, 
                   by = c("longitude" = "x", "latitude" = "y"),
                   match_fun = list(`>=`, `<=`))  # Use fuzzy matching with >= and <=
merged_data <- merged_data %>%
  mutate(distance = sqrt((longitude - x)^2 + (latitude - y)^2))  # Calculate distance
write.csv(merged_data, "28_58590.csv", row.names = FALSE)

# Find the row with the minimum distance for each population and select corresponding GO value
min_distance_rows <- merged_data %>%
  group_by(pop) %>%
  slice(which.min(distance)) %>%
  select(pop, Global_offset)
write.csv(min_distance_rows , "28_58590_min.csv", row.names = FALSE)





################################
##core west

# Read the frequency data
AllFreq <- read.table("core_west_AF.csv", sep = ",", header = TRUE, row.names = 1)

# Perform Principal Component Analysis (PCA)
pca <- rda(AllFreq, scale = TRUE)

# Plot the scree plot for PCA eigenvalues
screeplot(pca, type = "barplot", npcs = 10, main = "PCA Eigenvalues")

# Extract the first three principal components (PCs)
PCs <- scores(pca, choices = c(1:3), display = "sites", scaling = 0)

# Create a data frame containing the populations and their principal components
PopStruct <- data.frame(POP = row.names(AllFreq), PCs)

# Rename the columns of the data frame
colnames(PopStruct) <- c("POP", "PC1", "PC2", "PC3")

# Write the population structure data to a CSV file
write.csv(PopStruct, file = "PopStruct.csv", quote = FALSE)

# Combine coordinates, population structure, and environmental variables into one data frame
Variables <- data.frame(Coordinates, PopStruct[,-1], Env)

# Perform Redundancy Analysis (RDA)
RDA_outliers <- rda(AllFreq ~ BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15, Variables)

# Run the genomic_offset function for specified populations (58550 and 58590)
res_RDA_proj58590 <- genomic_offset(RDA_outliers, K = 3, env_pres = ras, env_fut = ras_58590, range = NULL, method = "loadings", scale_env = scale_env, center_env = center_env)
res_RDA_proj58550 <- genomic_offset(RDA_outliers, K = 3, env_pres = ras, env_fut = ras_58550, range = NULL, method = "loadings", scale_env = scale_env, center_env = center_env)
## Table global genetic offset predicted for 58550 and 58590
RDA_proj_offset <- data.frame(rbind(rasterToPoints(res_RDA_proj58550$Proj_offset_global), rasterToPoints(res_RDA_proj58590$Proj_offset_global)), Date = c(rep("58550", nrow(rasterToPoints(res_RDA_proj58550$Proj_offset_global))), rep("58590", nrow(rasterToPoints(res_RDA_proj58590$Proj_offset_global)))))
write.csv(RDA_proj_offset,file="RDA_proj_offset.csv",quote=F)



RDA_proj_offset_58550 <- read.csv("RDA_proj_offset_58550.csv")
raster_data <- rasterFromXYZ(RDA_proj_offset_58550[, c("x", "y", "Global_offset")])
asc_file <- "RDA_proj_offset_58550.asc"
writeRaster(raster_data, filename = asc_file, format = "ascii")
 

RDA_proj_offset_58590 <- read.csv("RDA_proj_offset_58590.csv")
raster_data <- rasterFromXYZ(RDA_proj_offset_58590[, c("x", "y", "Global_offset")])
asc_file <- "RDA_proj_offset_58590.asc"
writeRaster(raster_data, filename = asc_file, format = "ascii")

res_RDA_proj12690 <- genomic_offset(RDA_outliers, K = 3, env_pres = ras, env_fut = ras_12690, range = NULL, method = "loadings", scale_env = scale_env, center_env = center_env)
res_RDA_proj12650 <- genomic_offset(RDA_outliers, K = 3, env_pres = ras, env_fut = ras_12650, range = NULL, method = "loadings", scale_env = scale_env, center_env = center_env)
RDA_proj_offset <- data.frame(rbind(rasterToPoints(res_RDA_proj12650$Proj_offset_global), rasterToPoints(res_RDA_proj12690$Proj_offset_global)), Date = c(rep("12650", nrow(rasterToPoints(res_RDA_proj12650$Proj_offset_global))), rep("12690", nrow(rasterToPoints(res_RDA_proj12690$Proj_offset_global)))))
write.csv(RDA_proj_offset, "RDA_proj_offset_12650_90.csv", row.names = FALSE)

RDA_proj_offset_12650 <- read.csv("RDA_proj_offset_12650.csv")
raster_data <- rasterFromXYZ(RDA_proj_offset_12650[, c("x", "y", "Global_offset")])
asc_file <- "RDA_proj_offset_12650.asc"
writeRaster(raster_data, filename = asc_file, format = "ascii")
 
RDA_proj_offset_12690 <- read.csv("RDA_proj_offset_12690.csv")
raster_data <- rasterFromXYZ(RDA_proj_offset_12690[, c("x", "y", "Global_offset")])
asc_file <- "RDA_proj_offset_12690.asc"
writeRaster(raster_data, filename = asc_file, format = "ascii")


##########28pop_Future Climate GO
library(dplyr)
library(fuzzyjoin)

# Read longitude and latitude data file
df_pop <- read.csv("28poplonglat.csv", header = TRUE)

# Read GO values data file for 12650
df_go <- read.csv("RDA_proj_offset_12650.csv", header = TRUE)
merged_data <- df_pop %>%
     fuzzy_inner_join(df_go, 
                      by = c("longitude" = "x", "latitude" = "y"),
                      match_fun = list(`>=`, `<=`))  # Use 'greater than or equal' and 'less than or equal' for fuzzy matching					  
merged_data <- merged_data %>%
     mutate(distance = sqrt((longitude - x)^2 + (latitude - y)^2))  # Calculate distance
write.csv(merged_data, "28_12650.csv", row.names = FALSE)

# Find the row with the minimum distance for each population and select the corresponding GO value
min_distance_rows <- merged_data %>%
  group_by(pop) %>%
  slice(which.min(distance)) %>%
  select(pop, Global_offset)
write.csv(min_distance_rows , "28_12650_min.csv", row.names = FALSE)


# Read GO values data file for 12690
df_go <- read.csv("RDA_proj_offset_12690.csv", header = TRUE)
merged_data <- df_pop %>%
     fuzzy_inner_join(df_go, 
                      by = c("longitude" = "x", "latitude" = "y"),
                      match_fun = list(`>=`, `<=`))  # Use 'greater than or equal' and 'less than or equal' for fuzzy matching					  
merged_data <- merged_data %>%
     mutate(distance = sqrt((longitude - x)^2 + (latitude - y)^2))  # Calculate distance
write.csv(merged_data, "28_12690.csv", row.names = FALSE)

# Find the row with the minimum distance for each population and select the corresponding GO value
min_distance_rows <- merged_data %>%
  group_by(pop) %>%
  slice(which.min(distance)) %>%
  select(pop, Global_offset)
write.csv(min_distance_rows , "28_12690_min.csv", row.names = FALSE)



# Read GO values data file for 58550
df_go <- read.csv("RDA_proj_offset_58550.csv", header = TRUE)
merged_data <- df_pop %>%
     fuzzy_inner_join(df_go, 
                      by = c("longitude" = "x", "latitude" = "y"),
                      match_fun = list(`>=`, `<=`))  # Use 'greater than or equal' and 'less than or equal' for fuzzy matching					  
merged_data <- merged_data %>%
     mutate(distance = sqrt((longitude - x)^2 + (latitude - y)^2))  # Calculate distance
write.csv(merged_data, "28_58550.csv", row.names = FALSE)

# Find the row with the minimum distance for each population and select the corresponding GO value
min_distance_rows <- merged_data %>%
  group_by(pop) %>%
  slice(which.min(distance)) %>%
  select(pop, Global_offset)
write.csv(min_distance_rows , "28_58550_min.csv", row.names = FALSE)


# Read GO values data file for 58590
df_go <- read.csv("RDA_proj_offset_58590.csv", header = TRUE)
merged_data <- df_pop %>%
     fuzzy_inner_join(df_go, 
                      by = c("longitude" = "x", "latitude" = "y"),
                      match_fun = list(`>=`, `<=`))  # Use 'greater than or equal' and 'less than or equal' for fuzzy matching					  
merged_data <- merged_data %>%
     mutate(distance = sqrt((longitude - x)^2 + (latitude - y)^2))  # Calculate distance
write.csv(merged_data, "28_58590.csv", row.names = FALSE)

# Find the row with the minimum distance for each population and select the corresponding GO value
min_distance_rows <- merged_data %>%
  group_by(pop) %>%
  slice(which.min(distance)) %>%
  select(pop, Global_offset)
write.csv(min_distance_rows , "28_58590_min.csv", row.names = FALSE)




#####plot
library(sp)
library(raster)
DH5026_raster <- raster("RDA_proj_offset_12650.asc")
DH9026_raster <- raster("RDA_proj_offset_12690.asc")
DH5085_raster <- raster("RDA_proj_offset_58550.asc")
DH9085_raster <- raster("RDA_proj_offset_58590.asc")

pdf(file = "coreeastrdaDH_all_GO.pdf")
par(mfrow = c(2, 2), mar = c(2, 2, 2, 1), oma = c(1, 1, 1, 2))

# First subplot
plot(DH5026_raster, zlim = c(0.01, 7.5), legend = FALSE)  # Plot DH5026_raster, set z-axis range, hide the legend
points_data <- read.csv("coreeast_go.csv")  # Read GO value data from CSV file
library(sp)  # Load spatial data processing library

# Create spatial points data frame
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object
# Set coordinate reference system (assuming longitude and latitude use WGS84)
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# Calculate the maximum and minimum GO values
max_GO <- max(points_data$GO_5026)  # Get the maximum GO value
min_GO <- min(points_data$GO_5026)  # Get the minimum GO value

# Calculate symbol radius
points_data$radius <- 0.05 + 0.3 * (points_data$GO_5026 - min_GO) / (max_GO - min_GO)

# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(a)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label

# Second subplot
plot(DH9026_raster, zlim = c(0.01, 7.5), legend = FALSE)  # Plot DH9026_raster
points_data <- read.csv("coreeast_go.csv")  # Read GO value data from CSV file

# Create spatial points data frame
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")  # Set coordinate reference system

# Calculate the maximum and minimum GO values
max_GO <- max(points_data$GO_9026)
min_GO <- min(points_data$GO_9026)
points_data$radius <- 0.05 + 0.3 * (points_data$GO_9026 - min_GO) / (max_GO - min_GO)

# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(b)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label

# Third subplot
plot(DH5085_raster, zlim = c(0.01, 7.5), legend = FALSE)  # Plot DH5085_raster
points_data <- read.csv("coreeast_go.csv")  # Read GO value data from CSV file

# Create spatial points data frame
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")  # Set coordinate reference system

# Calculate the maximum and minimum GO values
max_GO <- max(points_data$GO_5085)
min_GO <- min(points_data$GO_5085)
points_data$radius <- 0.05 + 0.3 * (points_data$GO_5085 - min_GO) / (max_GO - min_GO)

# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(c)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label

# Fourth subplot
plot(DH9085_raster, zlim = c(0.01, 7.5))  # Plot DH9085_raster, and show the legend
points_data <- read.csv("coreeast_go.csv")  # Read GO value data from CSV file

# Create spatial points data frame
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")  # Set coordinate reference system

# Calculate the maximum and minimum GO values
max_GO <- max(points_data$GO_9085)
min_GO <- min(points_data$GO_9085)
points_data$radius <- 0.05 + 0.3 * (points_data$GO_9085 - min_GO) / (max_GO - min_GO)

# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(d)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label

dev.off()




library(sp)
library(raster)
DH5026_raster <- raster("RDA_proj_offset_12650.asc")
DH9026_raster <- raster("RDA_proj_offset_12690.asc")
DH5085_raster <- raster("RDA_proj_offset_58550.asc")
DH9085_raster <- raster("RDA_proj_offset_58590.asc")

pdf(file = "corewestrdaDH_all_GO.pdf")
par(mfrow = c(2, 2), mar = c(2, 2, 2, 1), oma = c(1, 1, 1, 2))

# First subplot
plot(DH5026_raster, zlim = c(0.01, 7.5), legend = FALSE)  # Plot DH5026_raster, set z-axis range, hide the legend
points_data <- read.csv("coewwest.csv")  # Read GO value data from the CSV file
library(sp)  # Load spatial data processing library

# Create spatial points data frame
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object

# Set coordinate reference system (assuming longitude and latitude use WGS84)
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# Calculate the maximum and minimum GO values
max_GO <- max(points_data$GO_5026)  # Get the maximum GO value
min_GO <- min(points_data$GO_5026)  # Get the minimum GO value

# Calculate symbol radius
points_data$radius <- 0.05 + 0.3 * (points_data$GO_5026 - min_GO) / (max_GO - min_GO)

# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(a)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label

# Second subplot
plot(DH9026_raster, zlim = c(0.01, 7.5), legend = FALSE)  # Plot DH9026_raster
points_data <- read.csv("coewwest.csv")  # Read GO value data from the CSV file

# Create spatial points data frame
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")  # Set coordinate reference system

# Calculate the maximum and minimum GO values
max_GO <- max(points_data$GO_9026)
min_GO <- min(points_data$GO_9026)
points_data$radius <- 0.05 + 0.3 * (points_data$GO_9026 - min_GO) / (max_GO - min_GO)

# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(b)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label

# Third subplot
plot(DH5085_raster, zlim = c(0.01, 7.5), legend = FALSE)  # Plot DH5085_raster
points_data <- read.csv("coewwest.csv")  # Read GO value data from the CSV file

# Create spatial points data frame
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")  # Set coordinate reference system

# Calculate the maximum and minimum GO values
max_GO <- max(points_data$GO_5085)
min_GO <- min(points_data$GO_5085)
points_data$radius <- 0.05 + 0.3 * (points_data$GO_5085 - min_GO) / (max_GO - min_GO)

# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(c)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label

# Fourth subplot
plot(DH9085_raster, zlim = c(0.01, 7.5))  # Plot DH9085_raster, and show the legend
points_data <- read.csv("coewwest.csv")  # Read GO value data from the CSV file

# Create spatial points data frame
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")  # Set coordinate reference system

# Calculate the maximum and minimum GO values
max_GO <- max(points_data$GO_9085)
min_GO <- min(points_data$GO_9085)
points_data$radius <- 0.05 + 0.3 * (points_data$GO_9085 - min_GO) / (max_GO - min_GO)

# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(d)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label

dev.off()



##pan east
DH5026_raster <- raster("RDA_proj_offset_12650.asc")
DH9026_raster <- raster("RDA_proj_offset_12690.asc")
DH5085_raster <- raster("RDA_proj_offset_58550.asc")
DH9085_raster <- raster("RDA_proj_offset_58590.asc")

pdf(file = "paneastrdaDH_all_GO.pdf")
par(mfrow = c(2, 2), mar = c(2, 2, 2, 1), oma = c(1, 1, 1, 2))

# First subplot
plot(DH5026_raster, zlim = c(0.009, 5.9), legend = FALSE)  # Plot DH5026_raster, set z-axis range, hide legend
points_data <- read.csv("paneast_go.csv")  # Read GO value data from CSV file
library(sp)  # Load spatial data processing library

# Create spatial points data frame
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object

# Set coordinate reference system (assuming longitude and latitude use WGS84)
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# Calculate the maximum and minimum GO values
max_GO <- max(points_data$GO_5026)  # Get the maximum GO value
min_GO <- min(points_data$GO_5026)  # Get the minimum GO value

# Calculate symbol radius
points_data$radius <- 0.05 + 0.3 * (points_data$GO_5026 - min_GO) / (max_GO - min_GO)

# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(a)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label

# Second subplot
plot(DH9026_raster, zlim = c(0.009, 5.9), legend = FALSE)  # Plot DH9026_raster
points_data <- read.csv("paneast_go.csv")  # Read GO value data from CSV file

# Create spatial points data frame
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")  # Set coordinate reference system

# Calculate the maximum and minimum GO values
max_GO <- max(points_data$GO_9026)
min_GO <- min(points_data$GO_9026)
points_data$radius <- 0.05 + 0.3 * (points_data$GO_9026 - min_GO) / (max_GO - min_GO)

# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(b)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label

# Third subplot
plot(DH5085_raster, zlim = c(0.009, 5.9), legend = FALSE)  # Plot DH5085_raster
points_data <- read.csv("paneast_go.csv")  # Read GO value data from CSV file

# Create spatial points data frame
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")  # Set coordinate reference system

# Calculate the maximum and minimum GO values
max_GO <- max(points_data$GO_5085)
min_GO <- min(points_data$GO_5085)
points_data$radius <- 0.05 + 0.3 * (points_data$GO_5085 - min_GO) / (max_GO - min_GO)

# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(c)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label

# Fourth subplot
plot(DH9085_raster, zlim = c(0.009, 5.9))  # Plot DH9085_raster, and show legend
points_data <- read.csv("paneast_go.csv")  # Read GO value data from CSV file

# Create spatial points data frame
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")  # Set coordinate reference system

# Calculate the maximum and minimum GO values
max_GO <- max(points_data$GO_9085)
min_GO <- min(points_data$GO_9085)
points_data$radius <- 0.05 + 0.3 * (points_data$GO_9085 - min_GO) / (max_GO - min_GO)

# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(d)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label

dev.off()


##pan west
DH5026_raster <- raster("RDA_proj_offset_12650.asc")
DH9026_raster <- raster("RDA_proj_offset_12690.asc")
DH5085_raster <- raster("RDA_proj_offset_58550.asc")
DH9085_raster <- raster("RDA_proj_offset_58590.asc")

pdf(file = "panwestrdaDH_all_GO.pdf")
par(mfrow = c(2, 2), mar = c(2, 2, 2, 1), oma = c(1, 1, 1, 2))

# First subplot
plot(DH5026_raster, zlim = c(0.009, 5.9), legend = FALSE)  # Plot DH5026_raster, set z-axis range, hide legend
points_data <- read.csv("panwest_go.csv")  # Read GO value data from CSV file
library(sp)  # Load spatial data processing library

# Create spatial points data frame
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object

# Set coordinate reference system (assuming longitude and latitude use WGS84)
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# Calculate the maximum and minimum GO values
max_GO <- max(points_data$GO_5026)  # Get the maximum GO value
min_GO <- min(points_data$GO_5026)  # Get the minimum GO value

# Calculate symbol radius
points_data$radius <- 0.05 + 0.3 * (points_data$GO_5026 - min_GO) / (max_GO - min_GO)

# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(a)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label

# Second subplot
plot(DH9026_raster, zlim = c(0.009, 5.9), legend = FALSE)  # Plot DH9026_raster
points_data <- read.csv("panwest_go.csv")  # Read GO value data from CSV file

# Create spatial points data frame
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")  # Set coordinate reference system

# Calculate the maximum and minimum GO values
max_GO <- max(points_data$GO_9026)
min_GO <- min(points_data$GO_9026)
points_data$radius <- 0.05 + 0.3 * (points_data$GO_9026 - min_GO) / (max_GO - min_GO)

# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(b)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label

# Third subplot
plot(DH5085_raster, zlim = c(0.009, 5.9), legend = FALSE)  # Plot DH5085_raster
points_data <- read.csv("panwest_go.csv")  # Read GO value data from CSV file

# Create spatial points data frame
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")  # Set coordinate reference system

# Calculate the maximum and minimum GO values
max_GO <- max(points_data$GO_5085)
min_GO <- min(points_data$GO_5085)
points_data$radius <- 0.05 + 0.3 * (points_data$GO_5085 - min_GO) / (max_GO - min_GO)

# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(c)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label

# Fourth subplot
plot(DH9085_raster, zlim = c(0.009, 5.9))  # Plot DH9085_raster, and show legend
points_data <- read.csv("panwest_go.csv")  # Read GO value data from CSV file

# Create spatial points data frame
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")  # Set coordinate reference system

# Calculate the maximum and minimum GO values
max_GO <- max(points_data$GO_9085)
min_GO <- min(points_data$GO_9085)
points_data$radius <- 0.05 + 0.3 * (points_data$GO_9085 - min_GO) / (max_GO - min_GO)

# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(d)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label

dev.off()



#######################################################################################################################
###Contemporary Adaptive Landscapes

source("C:\\Users\\feng\\Desktop\\景观\\RDA_GO\\RDA_GO_171_3144\\输入脚本\\adaptive_index.R")
library(pegas)
library(ggplot2)
library(raster)
library(rgdal)
library(LEA)
library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)
library(ggpubr)
library(vegan)
library(qvalue)
library(robust)
library(WMDB)
library(ggVennDiagram)
library(cowplot)
library(corrplot)
library(rgeos)
summary(RDA_outliers)
eigens <- eigenvals(RDA_outliers)
eigens <- RDA_outliers$CA$eig
total_variance <- sum(eigens)
variance_explained <- eigens / total_variance * 100
# View the proportion of variance explained by each principal component
variance_explained
total_variance <- sum(eigens)
# Calculate the proportion of variance explained by each principal component
variance_explained <- eigens / total_variance * 100
# Output the results
print(variance_explained)

# Extract species and environmental variable scores from RDA results
TAB_loci <- as.data.frame(scores(RDA_outliers, choices=c(1:2), display="species", scaling="none"))
TAB_var <- as.data.frame(scores(RDA_outliers, choices=c(1:2), display="bp"))

ggplot() +
  # Add horizontal and vertical reference lines
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  # Plot species points
  geom_point(data = TAB_loci, aes(x=RDA1*3, y=RDA2*3), colour = "#EB8055FF", size = 2, alpha = 0.8) +
  # Plot environmental variable arrows
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, 
               arrow=arrow(length = unit(0.02, "npc"))) + 
  # Add environmental variable labels
  geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  xlab("RDA 1 (Corresponding Value)") + ylab("RDA 2 (Corresponding Value)") +
  facet_wrap(~"Adaptively enriched RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), 
        strip.text = element_text(size=11))

# Calculate the adaptive index and extract RDA results
res_RDA_proj_current <- adaptive_index(RDA = RDA_outliers, K = 3, env_pres = ras, range = NULL, method = "loadings", scale_env = scale_env, center_env = center_env)
RDA_proj <- list(res_RDA_proj_current$RDA1, res_RDA_proj_current$RDA2)
# Convert RDA results into point data
RDA_proj <- lapply(RDA_proj, function(x) rasterToPoints(x))
# Normalize values
for(i in 1:length(RDA_proj)){
  RDA_proj[[i]][,3] <- (RDA_proj[[i]][,3]-min(RDA_proj[[i]][,3]))/(max(RDA_proj[[i]][,3])-min(RDA_proj[[i]][,3]))
}
# Combine RDA results into a data frame
TAB_RDA <- as.data.frame(do.call(rbind, RDA_proj[1:2]))
colnames(TAB_RDA)[3] <- "value"
TAB_RDA$variable <- factor(c(rep("RDA1", nrow(RDA_proj[[1]])), rep("RDA2", nrow(RDA_proj[[2]]))), levels = c("RDA1","RDA2"))

ggplot(data = TAB_RDA) + 
  geom_raster(aes(x = x, y = y, fill = cut(value, breaks=seq(0, 1, length.out=10), include.lowest = TRUE))) + 
  scale_fill_viridis_d(alpha = 0.8, direction = -1, option = "A", 
                       labels = c("Negative scores", "", "", "", "Intermediate scores", "", "", "", "Positive scores")) +
  coord_sf(xlim = c(85, 115), ylim = c(20, 45), expand = FALSE) +  
  xlab("Longitude") + 
  ylab("Latitude") +
  guides(fill = guide_legend(title = "Adaptive index")) +
  facet_grid(~ variable) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.grid = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        strip.text = element_text(size = 11))
















####### Mask Environmental Variables ######
library(raster)
library(rgdal)

# Read shapefile for masking --- east
shp <- raster("east.tif")
polygons <- rasterToPolygons(shp, dissolve = TRUE)
# Ensure the projection of the shapefile and raster are the same
proj4string(polygons) <- proj4string(shp)
# Get paths for all .asc files in the directory
asc_files <- list.files(pattern = "\\.asc$", full.names = TRUE)

# Loop through each .asc file and perform masking extraction
for (asc_file in asc_files) {
  # Read each .asc file
  r <- raster(asc_file)
  # Ensure the projection of the shapefile and raster are the same
  proj4string(polygons) <- proj4string(r)
  # Mask the raster using the shapefile
  masked_raster <- mask(r, polygons)
  # Create output filename (name based on input file)
  output_file <- gsub(".asc$", "_east.tif", basename(asc_file))
  # Save the masked result
  writeRaster(masked_raster, filename=output_file, format="GTiff", overwrite=TRUE)
  # Optional: Print processing status
  cat("Processed:", output_file, "\n")
}
# Output completion message
cat("All .asc files have been processed and saved as corresponding .tif files.\n")

# Read shapefile for masking --- west
shp <- raster("west.tif")
polygons <- rasterToPolygons(shp, dissolve = TRUE)
# Ensure the projection of the shapefile and raster are the same
proj4string(polygons) <- proj4string(shp)
# Get paths for all .asc files in the directory
asc_files <- list.files(pattern = "\\.asc$", full.names = TRUE)

# Loop through each .asc file and perform masking extraction
for (asc_file in asc_files) {
  # Read each .asc file
  r <- raster(asc_file)
  # Ensure the projection of the shapefile and raster are the same
  proj4string(polygons) <- proj4string(r)
  # Mask the raster using the shapefile
  masked_raster <- mask(r, polygons)
  # Create output filename (name based on input file)
  output_file <- gsub(".asc$", "_west.tif", basename(asc_file))
  # Save the masked result
  writeRaster(masked_raster, filename=output_file, format="GTiff", overwrite=TRUE)
  # Optional: Print processing status
  cat("Processed:", output_file, "\n")
}
# Output completion message
cat("All .asc files have been processed and saved as corresponding .tif files.\n")
