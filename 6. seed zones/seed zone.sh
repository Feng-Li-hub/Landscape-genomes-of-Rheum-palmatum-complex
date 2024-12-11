
###### Pan Zone
library(gradientForest)
gfData <- read.csv("pan_loci2_AF.csv")
envGF <- gfData[, 4:9]  # Read climate variable data from the file
snp <- gfData[, grep("RpChr", colnames(gfData))]  # Extract allele frequency data from columns with "RpChr" in the header
length(snp)  # Check the number of SNP sites extracted
#[1] 16709
maxLevel <- log2(0.368 * nrow(envGF) / 2)  # Account for correlations, see ?gradientForest
# Build the GF model
gfsnp_final <- gradientForest(cbind(envGF, snp), predictor.vars = colnames(envGF), response.vars = colnames(snp), ntree = 500, maxLevel = maxLevel, trace = TRUE, corr.threshold = 0.50)
DH_current_bio <- read.csv("test_115000_current+43pop.csv")
most_important <- c("BIO02", "BIO03", "BIO04", "BIO10", "BIO12", "BIO15")  ## Output the 6 most important BIO variables
imp.vars <- most_important
df_pca <- cbind(DH_current_bio[, c("long", "lat")], 
                predict(gfsnp_final, DH_current_bio[, imp.vars]))  # Combine coordinates and prediction results

# Principal Component Analysis (PCA)
PC <- prcomp(df_pca[, imp.vars])  
pcx <- PC$x  # Extract principal component scores

# Assign principal components
pc1 <- pcx[, 1]  
pc2 <- pcx[, 2]
pc3 <- pcx[, 3]

# Define RGB color palette
r <- pc2
g <- pc3 + pc1 - pc2
b <- pc3 - pc2
r <- (r - min(r)) / (max(r) - min(r))  # Normalize
g <- (g - min(g)) / (max(g) - min(g))
b <- (b - min(b)) / (max(b) - min(b))

# Plot bivariate plot
plot(pcx[, 1:2], pch = ".", cex = 1, col = rgb(r, g, b), asp = 1)

# Add arrows representing important climate variables
vec <- c("BIO02", "BIO03", "BIO04", "BIO10", "BIO12", "BIO15")
lv <- length(vec)
vind <- rownames(PC$rotation) %in% vec  # Select important variables from the principal components
class(vind) 
arrow1 <- PC$rotation[c(1, 2, 3, 4, 5, 6), 1]  # Direction for principal component 1
arrow2 <- PC$rotation[c(1, 2, 3, 4, 5, 6), 2]  # Direction for principal component 2
arrow_scale <- 10  # Set arrow length scale

# Redraw the plot with arrows
plot(pcx[, 1:2], pch = ".", cex = 2, col = rgb(r, g, b), asp = 1)
arrows(rep(0, lv), rep(0, lv), arrow1 / arrow_scale, arrow2 / arrow_scale, length = 0.0625)
jit <- 0.0012  # Distance between text and arrows
text(arrow1 / arrow_scale + jit * sign(arrow1), arrow2 / arrow_scale + jit * sign(arrow2), labels = vec)

# Plot spatial map showing GF prediction results
plot(df_pca[, c("long", "lat")], pch = ".", cex = 0.1, asp = 1, col = rgb(r, g, b))  # Plot spatial map

library("knitr")
library("ggplot2")
library("factoextra")
library("cluster")

# Perform clustering analysis using Clara algorithm and plot elbow plot to determine the optimal number of clusters
f <- fviz_nbclust(pcx, clara, method = "wss", k.max = 10)
# Extract variation data
variation <- f$data$y
write.csv(variation, file = "variation_all.csv", row.names = FALSE)
variation <- f$data$y  ## Find the elbow point

# Bar plot + line plot
k_values <- 1:length(variation)   # Define the number of clusters

# Set plot area
par(mar = c(5, 4, 4, 2) + 0.1)  # Adjust margins to fit both plots

# Draw bar plot for elbow method
barplot_heights <- barplot(variation, names.arg = k_values,
                            xlab = "Number of Clusters (K)", 
                            ylab = "Variation",
                            main = "Elbow Method",
                            col = rgb(0.2, 0.5, 0.8, 0.7),
                            ylim = c(0, max(variation) * 1.2))  # Set y-axis range

# Add line plot to the same graph
par(new = TRUE)  # Allow adding new plot on the same figure
plot(k_values, variation, type = "b", pch = 19, col = "black", lwd = 2,
     axes = FALSE, xlab = "", ylab = "", ylim = c(0, max(variation) * 1.2))

library(cluster)
# Cluster points into zones
#-------------------------
ncl <- 3  # The number of seed zones depends on the previous result
clPCs <- clara(pcx, ncl, sampsize = 10000)

# Set up the medoid color palette
medcolR <- r[clPCs$i.med]
medcolG <- g[clPCs$i.med]
medcolB <- b[clPCs$i.med]

summary(medcolR)
summary(medcolG)
summary(medcolB)

# Color adjustment
# Define colors for each zone
zone_colors <- c("coral1", "olivedrab3", "goldenrod1")

# PCA biplot into groups
plot(pcx[, 1:2], pch = ".", cex = 1, 
     col = zone_colors[clPCs$clustering], asp = 1)  # Use the defined colors for plotting
arrows(rep(0, lv), rep(0, lv), arrow1 / arrow_scale, arrow2 / arrow_scale, length = 0.0625)
text(arrow1 / arrow_scale + jit * sign(arrow1), arrow2 / arrow_scale + jit * sign(arrow2), labels = vec)



# Plot the PCA chart
plot(df_pca[, c("long", "lat")], 
     pch = ".", 
     cex = 0.1, 
     asp = 1, 
     col = zone_colors[clPCs$clustering],  # Use defined colors
     main = "", 
     xlim = range(df_pca$long),             # Set to the range of the data
     ylim = range(df_pca$lat))              # Set to the range of the data

# Add a legend
legend("bottomleft", 
       as.character(seq(1, ncl)), 
       pch = 15, 
       cex = 1, 
       col = zone_colors)                    # Use defined colors

# Mark the centroid points of the clusters
points(df_pca[clPCs$i.med, c("long", "lat")], 
       pch = as.character(seq(1, ncl)), 
       col = zone_colors)        

# Generate the RGB color mix for each location based on clustering results
mix <- rgb(medcolR[clPCs$clustering], medcolG[clPCs$clustering], medcolB[clPCs$clustering]) 
summary(mix)
str(mix)

# Get unique colors from the color vector
unique(mix)

# Replace specific color codes with numeric labels
mix2 <- gsub("#37ADC1", "1", mix)  # Replace color #37ADC1 with "1"
mix3 <- gsub("#84C1A0", "2", mix2)  # Replace color #84C1A0 with "2"
mix4 <- gsub("#965B98", "3", mix3)  # Replace color #965B98 with "3"

# Convert the modified character vector into numeric values
mix9 <- as.numeric(mix4)
mix9
unique(mix9)

# Extract latitude and longitude from the PCA dataframe
mixcluster <- df_pca[, c("long", "lat")]

# Add the color values to the mixcluster dataframe
mixcluster$color <- mix9

# Save the dataframe to a CSV file
write.csv(mixcluster, "rzone3_pan.csv")


####################################################################################################################################
####pan Seed zone migration
 gf <- gfsnp_final
 write.csv(df_pca, file="df_pca.csv", row.names=FALSE)
 vars <- most_important
 df_5026 <- read.csv("test_115000_ssp126_50_mean+43pop.csv")
 df_5026 <- df_5026[,2:9] 
names(df_5026)[1] <- "Longitude"
names(df_5026)[2] <- "Latitude"
df_5085 <- read.csv("test_115000_ssp585_50_mean+43pop.csv")
df_5085 <- df_5085[,2:9]
names(df_5085)[1] <- "Longitude"
names(df_5085)[2] <- "Latitude"
library(gradientForest)
##50-126
Trns_grid_5026 <- cbind(df_5026[,c("Longitude","Latitude")], predict(gf,df_5026[,vars]))
write.csv(Trns_grid_5026, file="Trns_grid_5026.csv", row.names=FALSE)
##50-585
Trns_grid_5085 <- cbind(df_5085[,c("Longitude","Latitude")], predict(gf,df_5085[,vars]))
write.csv(Trns_grid_5085, file="Trns_grid_5085.csv", row.names=FALSE)
#install.packages("FNN")
library(FNN)
# Convert data frames to matrices for FNN
pca_matrix <- as.matrix(df_pca[, vars])
future_matrix <- as.matrix(Trns_grid_5085[, vars])
#future_matrix26 <- as.matrix(Trns_grid_5026[, vars]) #2050 ssp126

# Finding the nearest neighbor
knn_result <- get.knnx(data = future_matrix, query = pca_matrix, k = 1)
#knn_result <- get.knnx(data = future_matrix26, query = pca_matrix, k = 1)

# Extract minimum distances and their indices
min_distances <- knn_result$nn.dist
min_indices <- knn_result$nn.index
# Create a dataframe with results
#SSP585
results_df <- data.frame(
  Original_Longitude = df_pca$long,
  Original_Latitude = df_pca$lat,
  Min_Distance = min_distances,
  Future_Longitude = Trns_grid_5085$Longitude[min_indices],
  Future_Latitude = Trns_grid_5085$Latitude[min_indices]
)

write.csv(results_df, "Seed_Zone_Shifts_under_climate_change_2050_ssp585_KNN_method_extended.csv", row.names = FALSE)


#ssp126
future_matrix26 <- as.matrix(Trns_grid_5026[, vars]) #2050 ssp126
knn_result <- get.knnx(data = future_matrix26, query = pca_matrix, k = 1)

# Extract minimum distances and their indices
min_distances <- knn_result$nn.dist
min_indices <- knn_result$nn.index
results26_df <- data.frame(
  Original_Longitude = df_pca$long,
  Original_Latitude = df_pca$lat,
  Min_Distance = min_distances,
  Future_Longitude = Trns_grid_5026$Longitude[min_indices],
  Future_Latitude = Trns_grid_5026$Latitude[min_indices]
)

write.csv(results26_df, "Seed_Zone_Shifts_under_climate_change_2050_ssp126_KNN_method_extended.csv", row.names = FALSE)

########## Plotting #########
library(readxl)
library(dplyr)
library(ggplot2)
library(geosphere)
library(rnaturalearth)
library(mapdata)
library(purrr)
library(sf)
library(raster)

#ssp585_50
two_zones_df <- read.csv("Seed_Zone_Shifts_under_climate_change_2050_ssp585_KNN_method_extended.csv")
centroids <- two_zones_df %>%
  group_by(Seed_Zone) %>%
  summarise(
    Orig_Centroid_Lon = mean(Original_Longitude),
    Orig_Centroid_Lat = mean(Original_Latitude),
    Future_Centroid_Lon = mean(Future_Longitude),
    Future_Centroid_Lat = mean(Future_Latitude),
    mean_offset = mean(Min_Distance)
  )
centroids$Distance_km <- mapply(function(lon1, lat1, lon2, lat2) {
  distHaversine(c(lon1, lat1), c(lon2, lat2)) / 1000  # convert from meters to kilometers
}, centroids$Orig_Centroid_Lon, centroids$Orig_Centroid_Lat, centroids$Future_Centroid_Lon, centroids$Future_Centroid_Lat)
### The mapply function is used to apply the given function (in this case, an anonymous function) to corresponding elements of multiple lists or vectors. 
## It takes four parameters: Orig_Centroid_Lon, Orig_Centroid_Lat, Future_Centroid_Lon, and Future_Centroid_Lat, representing the start and end coordinates (longitude and latitude).
## The anonymous function uses the distHaversine function to calculate the great-circle distance between two latitude and longitude coordinates (in meters), 
## and then the result is divided by 1000 to convert the distance to kilometers.
## centroids$Distance_km is assigned the great-circle distance between each pair of start and end points (in kilometers).

arrow_data_50585 <- centroids %>%
  transmute(
    x = Orig_Centroid_Lon, y = Orig_Centroid_Lat,
    xend = Future_Centroid_Lon, yend = Future_Centroid_Lat,
    Seed_Zone)
## centroids %>%: The pipe operator %>% passes the centroids dataframe to the next function.
## transmute function is used to select and rename columns in the dataframe. Here, it selects four columns: Orig_Centroid_Lon is renamed as x, Orig_Centroid_Lat as y, 
## Future_Centroid_Lon as xend, and Future_Centroid_Lat as yend. It also includes the Seed_Zone column.
write.csv(arrow_data_50585,"arrow_data_50585.csv")

# Get geographical data for China
china <- ne_countries(country = "china", returnclass = "sf")
# Get provincial boundary data for China
china_provinces <- ne_states(country = "china", returnclass = "sf")

arrow_data <- arrow_data_50585
#"coral1", "olivedrab3", "goldenrod1"

ggplot() +
  geom_sf(data = china, fill = "grey90", color = "black") +  # Plot the background of China's map
  geom_sf(data = china_provinces, fill = NA, color = "black", lwd = 0.5) +  # Plot provincial boundaries
  geom_point(aes(x = two_zones_df$Future_Longitude, y = two_zones_df$Future_Latitude, col = as.factor(two_zones_df$Seed_Zone)), pch = ".", cex = 0.01) +
  scale_color_manual(values = c("1" = "coral1", "2" = "olivedrab3", "3" = "goldenrod1")) +
  geom_point(aes(x = centroids[1,]$Orig_Centroid_Lon, y = centroids[1,]$Orig_Centroid_Lat), 
             pch = 21, size = 7, fill = "coral1", color = "black", stroke = 1.2) +
  geom_point(aes(x = centroids[2,]$Orig_Centroid_Lon, y = centroids[2,]$Orig_Centroid_Lat), 
             pch = 21, size = 7, fill = "olivedrab3", color = "black", stroke = 1.2) +
  geom_point(aes(x = centroids[3,]$Orig_Centroid_Lon, y = centroids[3,]$Orig_Centroid_Lat), 
             pch = 21, size = 7, fill = "goldenrod1", color = "black", stroke = 1.2) +
  geom_point(aes(x = centroids[1,]$Future_Centroid_Lon, y = centroids[1,]$Future_Centroid_Lat), 
             pch = 21, size = 7, fill = "white", color = "black", stroke = 1.2) +
  geom_point(aes(x = centroids[2,]$Future_Centroid_Lon, y = centroids[2,]$Future_Centroid_Lat), 
             pch = 21, size = 7, fill = "white", color = "black", stroke = 1.2) +
  geom_point(aes(x = centroids[3,]$Future_Centroid_Lon, y = centroids[3,]$Future_Centroid_Lat), 
             pch = 21, size = 7, fill = "white", color = "black", stroke = 1.2) + 
  geom_point(aes(x = centroids$Orig_Centroid_Lon, y = centroids$Orig_Centroid_Lat), pch = c("1", "2", "3"), size = 5, col = "white") +
  geom_point(aes(x = centroids[1,]$Future_Centroid_Lon, y = centroids[1,]$Future_Centroid_Lat), pch = c("1"), size = 5, col = "coral1") +
  geom_point(aes(x = centroids[2,]$Future_Centroid_Lon, y = centroids[2,]$Future_Centroid_Lat), pch = c("2"), size = 5, col = "olivedrab3") +
  geom_point(aes(x = centroids[3,]$Future_Centroid_Lon, y = centroids[3,]$Future_Centroid_Lat), pch = c("3"), size = 5, col = "goldenrod1") +
  geom_segment(data = arrow_data, aes(x = x + 0.53, y = y + 0.23, xend = (xend - 0.5), yend = (yend - 0.5)), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black", size = 0.8) +
  coord_sf(xlim = c(85, 115), ylim = c(20, 45), expand = FALSE) +
  xlab("\nLongitude") + ylab("Latitude\n") +
  theme_bw(base_size = 11) +
  theme(plot.background = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        panel.grid.major = element_blank(), legend.position = "none")

#ssp126_50
# Load data for SSP126 scenario
two_zones_df <- read.csv("Seed_Zone_Shifts_under_climate_change_2050_ssp126_KNN_method_extended.csv")

# Group by Seed_Zone and calculate centroids and distance for each seed zone
centroids <- two_zones_df %>%
  group_by(Seed_Zone) %>%
  summarise(
    Orig_Centroid_Lon = mean(Original_Longitude),  # Calculate the mean of the original longitude
    Orig_Centroid_Lat = mean(Original_Latitude),   # Calculate the mean of the original latitude
    Future_Centroid_Lon = mean(Future_Longitude),  # Calculate the mean of the future longitude
    Future_Centroid_Lat = mean(Future_Latitude),   # Calculate the mean of the future latitude
    mean_offset = mean(Min_Distance)               # Calculate the average of the minimum distance
  )

# Calculate the distance (in kilometers) between the original and future centroids
centroids$Distance_km <- mapply(function(lon1, lat1, lon2, lat2) {
  distHaversine(c(lon1, lat1), c(lon2, lat2)) / 1000  # Using Haversine formula to calculate distance, converted to kilometers
}, centroids$Orig_Centroid_Lon, centroids$Orig_Centroid_Lat, centroids$Future_Centroid_Lon, centroids$Future_Centroid_Lat)

# Create a new data frame for arrow plotting data
arrow_data_50126 <- centroids %>%
  transmute(
    x = Orig_Centroid_Lon,            # Original Longitude
    y = Orig_Centroid_Lat,            # Original Latitude
    xend = Future_Centroid_Lon,       # Future Longitude
    yend = Future_Centroid_Lat,       # Future Latitude
    Seed_Zone                         # Seed Zone Information
)

# Save the arrow data to a CSV file
write.csv(arrow_data_50126, "arrow_data_50126.csv")

# Assign the arrow data to a new variable for further use
arrow_data <- arrow_data_50126

# Plotting the data using ggplot2
ggplot() +
  geom_sf(data = china, fill = "grey90", color = "black") +  # Plot the background of China's map
  geom_sf(data = china_provinces, fill = NA, color = "black", lwd = 0.5) +  # Plot provincial boundaries
  geom_point(aes(x = two_zones_df$Future_Longitude, y = two_zones_df$Future_Latitude, col = as.factor(two_zones_df$Seed_Zone)), pch = ".", cex = 0.01) +  # Plot points for future locations, colored by seed zone
  scale_color_manual(values = c("1" = "coral1", "2" = "olivedrab3", "3" = "goldenrod1")) +  # Assign colors to seed zones
  geom_point(aes(x = centroids[1,]$Orig_Centroid_Lon, y = centroids[1,]$Orig_Centroid_Lat), 
             pch = 21, size = 7, fill = "coral1", color = "black", stroke = 1.2) +  # Plot the original centroid for Seed Zone 1
  geom_point(aes(x = centroids[2,]$Orig_Centroid_Lon, y = centroids[2,]$Orig_Centroid_Lat), 
             pch = 21, size = 7, fill = "olivedrab3", color = "black", stroke = 1.2) +  # Plot the original centroid for Seed Zone 2
  geom_point(aes(x = centroids[3,]$Orig_Centroid_Lon, y = centroids[3,]$Orig_Centroid_Lat), 
             pch = 21, size = 7, fill = "goldenrod1", color = "black", stroke = 1.2) +  # Plot the original centroid for Seed Zone 3
  geom_point(aes(x = centroids[1,]$Future_Centroid_Lon, y = centroids[1,]$Future_Centroid_Lat), 
             pch = 21, size = 7, fill = "white", color = "black", stroke = 1.2) +  # Plot the future centroid for Seed Zone 1
  geom_point(aes(x = centroids[2,]$Future_Centroid_Lon, y = centroids[2,]$Future_Centroid_Lat), 
             pch = 21, size = 7, fill = "white", color = "black", stroke = 1.2) +  # Plot the future centroid for Seed Zone 2
  geom_point(aes(x = centroids[3,]$Future_Centroid_Lon, y = centroids[3,]$Future_Centroid_Lat), 
             pch = 21, size = 7, fill = "white", color = "black", stroke = 1.2) +  # Plot the future centroid for Seed Zone 3
  geom_point(aes(x = centroids$Orig_Centroid_Lon, y = centroids$Orig_Centroid_Lat), pch = c("1", "2", "3"), size = 5, col = "white") +  # Highlight original centroids
  geom_point(aes(x = centroids[1,]$Future_Centroid_Lon, y = centroids[1,]$Future_Centroid_Lat), pch = c("1"), size = 5, col = "coral1") +  # Highlight future centroid for Seed Zone 1
  geom_point(aes(x = centroids[2,]$Future_Centroid_Lon, y = centroids[2,]$Future_Centroid_Lat), pch = c("2"), size = 5, col = "olivedrab3") +  # Highlight future centroid for Seed Zone 2
  geom_point(aes(x = centroids[3,]$Future_Centroid_Lon, y = centroids[3,]$Future_Centroid_Lat), pch = c("3"), size = 5, col = "goldenrod1") +  # Highlight future centroid for Seed Zone 3
  geom_segment(data = arrow_data, aes(x = x + 0.53, y = y + 0.23, xend = (xend - 0.5), yend = (yend - 0.5)), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black", size = 0.8) +  # Plot arrows between original and future centroids
  coord_sf(xlim = c(85, 115), ylim = c(20, 45), expand = FALSE) +  # Set the map boundaries (longitude and latitude limits)
  xlab("\nLongitude") + ylab("Latitude\n") +  # Label the axes
  theme_bw(base_size = 11) +  # Set theme with a clean white background
  theme(plot.background = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        panel.grid.major = element_blank(), legend.position = "none")  # Further customize the theme

######90_126、585
df_9026 <- read.csv("test_115000_ssp126_90_mean+43pop.csv")
df_9026 <- df_9026[,2:9] 
names(df_9026)[1] <- "Longitude"
names(df_9026)[2] <- "Latitude"
df_9085 <- read.csv("test_115000_ssp585_90_mean+43pop.csv")
df_9085 <- df_9085[,2:9]
names(df_9085)[1] <- "Longitude"
names(df_9085)[2] <- "Latitude"
library(gradientForest)
##90-126
Trns_grid_9026 <- cbind(df_9026[,c("Longitude","Latitude")], predict(gf,df_9026[,vars]))
write.csv(Trns_grid_9026, file="Trns_grid_9026.csv", row.names=FALSE)
##90-585
Trns_grid_9085 <- cbind(df_9085[,c("Longitude","Latitude")], predict(gf,df_9085[,vars]))
write.csv(Trns_grid_9085, file="Trns_grid_9085.csv", row.names=FALSE)

#install.packages("FNN")
library(FNN)
# Convert data frames to matrices for FNN
pca_matrix <- as.matrix(df_pca[, vars])
future_matrix <- as.matrix(Trns_grid_9085[, vars])
future_matrix26 <- as.matrix(Trns_grid_9026[, vars]) #2090 ssp126

# Finding the nearest neighbor
knn_result <- get.knnx(data = future_matrix, query = pca_matrix, k = 1)
#knn_result <- get.knnx(data = future_matrix26, query = pca_matrix, k = 1)

# Extract minimum distances and their indices
min_distances <- knn_result$nn.dist
min_indices <- knn_result$nn.index

results_df <- data.frame(
  Original_Longitude = df_pca$long,
  Original_Latitude = df_pca$lat,
  Min_Distance = min_distances,
  Future_Longitude = Trns_grid_9085$Longitude[min_indices],
  Future_Latitude = Trns_grid_9085$Latitude[min_indices]
)

write.csv(results_df, "Seed_Zone_Shifts_under_climate_change_2090_ssp585_KNN_method_extended.csv", row.names = FALSE)

#ssp126

knn_result <- get.knnx(data = future_matrix26, query = pca_matrix, k = 1)

# Extract minimum distances and their indices
min_distances <- knn_result$nn.dist
min_indices <- knn_result$nn.index
results26_df <- data.frame(
  Original_Longitude = df_pca$long,
  Original_Latitude = df_pca$lat,
  Min_Distance = min_distances,
  Future_Longitude = Trns_grid_9026$Longitude[min_indices],
  Future_Latitude = Trns_grid_9026$Latitude[min_indices]
)

write.csv(results26_df, "Seed_Zone_Shifts_under_climate_change_2090_ssp126_KNN_method_extended.csv", row.names = FALSE)

#ssp126_90

# Load the CSV file for 2090 projections under SSP126 scenario using the KNN method
two_zones_df <- read.csv("Seed_Zone_Shifts_under_climate_change_2090_ssp126_KNN_method_extended.csv")

# Group by Seed_Zone and calculate centroids and distance metrics for each seed zone
centroids <- two_zones_df %>%
  group_by(Seed_Zone) %>%
  summarise(
    Orig_Centroid_Lon = mean(Original_Longitude),  # Calculate the mean original longitude
    Orig_Centroid_Lat = mean(Original_Latitude),   # Calculate the mean original latitude
    Future_Centroid_Lon = mean(Future_Longitude),   # Calculate the mean future longitude
    Future_Centroid_Lat = mean(Future_Latitude),    # Calculate the mean future latitude
    mean_offset = mean(Min_Distance)                # Calculate the mean of the minimum distances
  )

# Calculate the distance (in km) between the original and future centroids
centroids$Distance_km <- mapply(function(lon1, lat1, lon2, lat2) {
  distHaversine(c(lon1, lat1), c(lon2, lat2)) / 1000  # Use Haversine formula for distance calculation (in kilometers)
}, centroids$Orig_Centroid_Lon, centroids$Orig_Centroid_Lat, centroids$Future_Centroid_Lon, centroids$Future_Centroid_Lat)

# Create a new dataframe containing the coordinates for arrows between original and future centroids
arrow_data_90126 <- centroids %>%
  transmute(
    x = Orig_Centroid_Lon,            # Original centroid longitude
    y = Orig_Centroid_Lat,             # Original centroid latitude
    xend = Future_Centroid_Lon,       # Future centroid longitude
    yend = Future_Centroid_Lat,       # Future centroid latitude
    Seed_Zone                         # Seed zone information
  )

# Save the arrow data to a CSV file
write.csv(arrow_data_90126, "arrow_data_90126.csv")
arrow_data <- arrow_data_90126

# Create the plot
ggplot() +
  geom_sf(data = china, fill = "grey90", color = "black") +  # Plot the map of China
  geom_sf(data = china_provinces, fill = NA, color = "black", lwd = 0.5) +  # Plot the provincial borders
  geom_point(aes(x = two_zones_df$Future_Longitude, y = two_zones_df$Future_Latitude, col = as.factor(two_zones_df$Seed_Zone)), pch = ".", cex = 0.01) +  # Plot the points for future seed zone locations
  scale_color_manual(values = c("1" = "coral1", "2" = "olivedrab3", "3" = "goldenrod1")) +  # Customize color scale for seed zones
  # Plot original centroids for each seed zone
  geom_point(aes(x = centroids[1,]$Orig_Centroid_Lon, y = centroids[1,]$Orig_Centroid_Lat), pch = 21, size = 7, fill = "coral1", color = "black", stroke = 1.2) +
  geom_point(aes(x = centroids[2,]$Orig_Centroid_Lon, y = centroids[2,]$Orig_Centroid_Lat), pch = 21, size = 7, fill = "olivedrab3", color = "black", stroke = 1.2) +
  geom_point(aes(x = centroids[3,]$Orig_Centroid_Lon, y = centroids[3,]$Orig_Centroid_Lat), pch = 21, size = 7, fill = "goldenrod1", color = "black", stroke = 1.2) +
  # Plot future centroids for each seed zone
  geom_point(aes(x = centroids[1,]$Future_Centroid_Lon, y = centroids[1,]$Future_Centroid_Lat), pch = 21, size = 7, fill = "white", color = "black", stroke = 1.2) +
  geom_point(aes(x = centroids[2,]$Future_Centroid_Lon, y = centroids[2,]$Future_Centroid_Lat), pch = 21, size = 7, fill = "white", color = "black", stroke = 1.2) +
  geom_point(aes(x = centroids[3,]$Future_Centroid_Lon, y = centroids[3,]$Future_Centroid_Lat), pch = 21, size = 7, fill = "white", color = "black", stroke = 1.2) +
  # Add points for seed zone identifiers
  geom_point(aes(x = centroids$Orig_Centroid_Lon, y = centroids$Orig_Centroid_Lat), pch = c("1", "2", "3"), size = 5, col = "white") +
  geom_point(aes(x = centroids[1,]$Future_Centroid_Lon, y = centroids[1,]$Future_Centroid_Lat), pch = c("1"), size = 5, col = "coral1") +
  geom_point(aes(x = centroids[2,]$Future_Centroid_Lon, y = centroids[2,]$Future_Centroid_Lat), pch = c("2"), size = 5, col = "olivedrab3") +
  geom_point(aes(x = centroids[3,]$Future_Centroid_Lon, y = centroids[3,]$Future_Centroid_Lat), pch = c("3"), size = 5, col = "goldenrod1") +
  # Add arrows showing the shift from original to future centroids
  geom_segment(data = arrow_data, aes(x = x + 0.53, y = y + 0.23, xend = (xend - 0.5), yend = (yend - 0.5)),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", size = 0.8) +
  # Set the map view limits and labels
  coord_sf(xlim = c(85, 115), ylim = c(20, 45), expand = FALSE) +
  xlab("\nLongitude") + ylab("Latitude\n") +
  # Customize the theme of the plot
  theme_bw(base_size = 11) +
  theme(plot.background = element_blank(), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14),
        panel.grid.major = element_blank(), 
        legend.position = "none")

 ##58590

# Load the CSV file for 2090 projections under SSP585 scenario using the KNN method
two_zones_df <- read.csv("Seed_Zone_Shifts_under_climate_change_2090_ssp585_KNN_method_extended.csv")

# Group by Seed_Zone and calculate centroids and distance metrics for each seed zone
centroids <- two_zones_df %>%
  group_by(Seed_Zone) %>%
  summarise(
    Orig_Centroid_Lon = mean(Original_Longitude),  # Calculate the mean of original longitude
    Orig_Centroid_Lat = mean(Original_Latitude),   # Calculate the mean of original latitude
    Future_Centroid_Lon = mean(Future_Longitude),   # Calculate the mean of future longitude
    Future_Centroid_Lat = mean(Future_Latitude),    # Calculate the mean of future latitude
    mean_offset = mean(Min_Distance)                # Calculate the mean of the minimum distances
  )

# Calculate the distance (in km) between the original and future centroids
centroids$Distance_km <- mapply(function(lon1, lat1, lon2, lat2) {
  distHaversine(c(lon1, lat1), c(lon2, lat2)) / 1000  # Use Haversine formula to calculate distance in kilometers
}, centroids$Orig_Centroid_Lon, centroids$Orig_Centroid_Lat, centroids$Future_Centroid_Lon, centroids$Future_Centroid_Lat)

# Create a new dataframe containing the coordinates for the arrows between original and future centroids
arrow_data_90585 <- centroids %>%
  transmute(
    x = Orig_Centroid_Lon,            # Original longitude
    y = Orig_Centroid_Lat,             # Original latitude
    xend = Future_Centroid_Lon,       # Future longitude
    yend = Future_Centroid_Lat,       # Future latitude
    Seed_Zone                         # Seed zone information
  )

# Save the arrow data to a CSV file
write.csv(arrow_data_90585, "arrow_data_90585.csv")
arrow_data <- arrow_data_90585

# Create the plot
ggplot() +
  geom_sf(data = china, fill = "grey90", color = "black") +  # Plot the map of China
  geom_sf(data = china_provinces, fill = NA, color = "black", lwd = 0.5) +  # Plot the provincial borders
  geom_point(aes(x = two_zones_df$Future_Longitude, y = two_zones_df$Future_Latitude, col = as.factor(two_zones_df$Seed_Zone)), pch = ".", cex = 0.01) +  # Plot points for future seed zone locations
  scale_color_manual(values = c("1" = "coral1", "2" = "olivedrab3", "3" = "goldenrod1")) +  # Set colors for different seed zones
  # Plot original centroids for each seed zone
  geom_point(aes(x = centroids[1,]$Orig_Centroid_Lon, y = centroids[1,]$Orig_Centroid_Lat), pch = 21, size = 7, fill = "coral1", color = "black", stroke = 1.2) +
  geom_point(aes(x = centroids[2,]$Orig_Centroid_Lon, y = centroids[2,]$Orig_Centroid_Lat), pch = 21, size = 7, fill = "olivedrab3", color = "black", stroke = 1.2) +
  geom_point(aes(x = centroids[3,]$Orig_Centroid_Lon, y = centroids[3,]$Orig_Centroid_Lat), pch = 21, size = 7, fill = "goldenrod1", color = "black", stroke = 1.2) +
  # Plot future centroids for each seed zone
  geom_point(aes(x = centroids[1,]$Future_Centroid_Lon, y = centroids[1,]$Future_Centroid_Lat), pch = 21, size = 7, fill = "white", color = "black", stroke = 1.2) +
  geom_point(aes(x = centroids[2,]$Future_Centroid_Lon, y = centroids[2,]$Future_Centroid_Lat), pch = 21, size = 7, fill = "white", color = "black", stroke = 1.2) +
  geom_point(aes(x = centroids[3,]$Future_Centroid_Lon, y = centroids[3,]$Future_Centroid_Lat), pch = 21, size = 7, fill = "white", color = "black", stroke = 1.2) +
  # Add points for seed zone identifiers
  geom_point(aes(x = centroids$Orig_Centroid_Lon, y = centroids$Orig_Centroid_Lat), pch = c("1", "2", "3"), size = 5, col = "white") +
  geom_point(aes(x = centroids[1,]$Future_Centroid_Lon, y = centroids[1,]$Future_Centroid_Lat), pch = c("1"), size = 5, col = "coral1") +
  geom_point(aes(x = centroids[2,]$Future_Centroid_Lon, y = centroids[2,]$Future_Centroid_Lat), pch = c("2"), size = 5, col = "olivedrab3") +
  geom_point(aes(x = centroids[3,]$Future_Centroid_Lon, y = centroids[3,]$Future_Centroid_Lat), pch = c("3"), size = 5, col = "goldenrod1") +
  # Add arrows showing the shift from original to future centroids
  geom_segment(data = arrow_data, aes(x = x + 0.53, y = y + 0.23, xend = (xend - 0.5), yend = (yend - 0.5)),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", size = 0.8) +
  # Set the map view limits and labels
  coord_sf(xlim = c(85, 115), ylim = c(20, 45), expand = FALSE) +
  xlab("\nLongitude") + ylab("Latitude\n") +
  # Customize the theme of the plot
  theme_bw(base_size = 11) +
  theme(plot.background = element_blank(), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14),
        panel.grid.major = element_blank(), 
        legend.position = "none")






##############################################################################
###core
library(gradientForest)
gfData <- read.csv("core_loci2_AF.csv")
envGF <- gfData[, 4:9] # Read climate variables from the file
snp <- gfData[, grep("RpChr", colnames(gfData))] # Read allele frequency data based on column names
length(snp) # Check the number of SNPs extracted
#[1] 16709
maxLevel <- log2(0.368 * nrow(envGF) / 2) # Account for correlations, see ?gradientForest

# Build the GF model
gfsnp_final <- gradientForest(cbind(envGF, snp), predictor.vars = colnames(envGF), 
                              response.vars = colnames(snp), ntree = 500, maxLevel = maxLevel, 
                              trace = T, corr.threshold = 0.50)

DH_current_bio <- read.csv("test_115000_current+43pop.csv")
most_important <- c("BIO02", "BIO03", "BIO04", "BIO10", "BIO12", "BIO15") # Output the top 6 important BIO variables
imp.vars <- most_important
df_pca <- cbind(DH_current_bio[, c("long", "lat")], 
                predict(gfsnp_final, DH_current_bio[, imp.vars]))  # Combine coordinates with prediction results

# Principal Component Analysis (PCA)
PC <- prcomp(df_pca[, imp.vars])  
pcx <- PC$x  # Extract principal component scores

# Assign principal components
pc1 <- pcx[, 1]  
pc2 <- pcx[, 2]
pc3 <- pcx[, 3]

# Define RGB color palette
r <- pc2
g <- pc3 + pc1 - pc2
b <- pc3 - pc2
r <- (r - min(r)) / (max(r) - min(r))  # Normalize
g <- (g - min(g)) / (max(g) - min(g))
b <- (b - min(b)) / (max(b) - min(b))

# Plot bivariate graph
plot(pcx[, 1:2], pch = ".", cex = 1, col = rgb(r, g, b), asp = 1)

# Add arrows to indicate important climate variables
vec <- c("BIO02", "BIO03", "BIO04", "BIO10", "BIO12", "BIO15")
lv <- length(vec)
vind <- rownames(PC$rotation) %in% vec  # Select corresponding important variables in the principal components
class(vind)
arrow1 <- PC$rotation[c(1, 2, 3, 4, 5, 6), 1]  # Direction of principal component 1
arrow2 <- PC$rotation[c(1, 2, 3, 4, 5, 6), 2]  # Direction of principal component 2
arrow_scale <- 10  # Set arrow length scale

# Plot again with arrows
plot(pcx[, 1:2], pch = ".", cex = 2, col = rgb(r, g, b), asp = 1)
arrows(rep(0, lv), rep(0, lv), arrow1 / arrow_scale, arrow2 / arrow_scale, length = 0.0625)
jit <- 0.0012  # Distance between text and arrows
text(arrow1 / arrow_scale + jit * sign(arrow1), arrow2 / arrow_scale + jit * sign(arrow2), labels = vec)

# Plot spatial map to show GF prediction results
plot(df_pca[, c("long", "lat")], pch = ".", cex = 0.1, asp = 1, col = rgb(r, g, b))  # Spatial map

library("knitr")
library("ggplot2")
library("factoextra")
library("cluster")

# Use Clara algorithm for clustering analysis and plot the elbow method to determine the optimal number of clusters
f <- fviz_nbclust(pcx, clara, method = "wss", k.max = 10)
# Extract variation data
variation <- f$data$y
write.csv(variation, file = "variation_all.csv", row.names = FALSE)
variation <- f$data$y ## Find the elbow point

# Bar plot + line plot
k_values <- 1:length(variation)   # Define the number of clusters

# Set up the plotting area
par(mar = c(5, 4, 4, 2) + 0.1)  # Adjust margins to fit both plots

# Bar plot for elbow method
barplot_heights <- barplot(variation, names.arg = k_values,
                            xlab = "Number of Clusters (K)", 
                            ylab = "Variation",
                            main = "Elbow Method",
                            col = rgb(0.2, 0.5, 0.8, 0.7),
                            ylim = c(0, max(variation) * 1.2))  # Set y-axis range

# Add line plot on top
par(new = TRUE)  # Allow adding new plot to the same figure
plot(k_values, variation, type = "b", pch = 19, col = "black", lwd = 2,
     axes = FALSE, xlab = "", ylab = "", ylim = c(0, max(variation) * 1.2))

library(cluster)
# Cluster points into zones
#-------------------------
ncl <- 3 # The number of seed zones depends on the previous result
clPCs <- clara(pcx, ncl, sampsize = 10000)

# Set up the medoid color palette
medcolR <- r[clPCs$i.med]
medcolG <- g[clPCs$i.med]
medcolB <- b[clPCs$i.med]

summary(medcolR)
summary(medcolG)
summary(medcolB)

# Define colors for each zone
zone_colors <- c("coral1", "olivedrab3", "goldenrod1")

# PCA biplot into groups
plot(pcx[, 1:2], pch = ".", cex = 1, 
     col = zone_colors[clPCs$clustering], asp = 1)  # Use the defined colors for plotting
arrows(rep(0, lv), rep(0, lv), arrow1 / arrow_scale, arrow2 / arrow_scale, length = 0.0625)
text(arrow1 / arrow_scale + jit * sign(arrow1), arrow2 / arrow_scale + jit * sign(arrow2), labels = vec)



# Plot PCA graph
plot(df_pca[, c("long", "lat")], 
     pch = ".", 
     cex = 0.1, 
     asp = 1, 
     col = zone_colors[clPCs$clustering],  # Use the defined colors
     main = "", 
     xlim = range(df_pca$long),             # Set the range to data range
     ylim = range(df_pca$lat))              # Set the range to data range

# Add legend
legend("bottomleft", 
       as.character(seq(1, ncl)), 
       pch = 15, 
       cex = 1, 
       col = zone_colors)                    # Use the defined colors

# Mark cluster medoid points
points(df_pca[clPCs$i.med, c("long", "lat")], 
       pch = as.character(seq(1, ncl)), 
       col = zone_colors)        

# Generate RGB color mix for each location based on clustering results
mix <- rgb(medcolR[clPCs$clustering], medcolG[clPCs$clustering], medcolB[clPCs$clustering]) 
summary(mix)
str(mix)

# Get unique color values from the color vector
unique(mix)

# Replace specific color codes with numeric labels
mix2 <- gsub("#BD3635", "1", mix)  # Replace color #BD3635 with "1"
mix3 <- gsub("#693852", "2", mix2)  # Replace color #693852 with "2"
mix4 <- gsub("#72B14E", "3", mix3)  # Replace color #72B14E with "3"

# Convert the replaced character vector to numeric
mix9 <- as.numeric(mix4)
mix9
unique(mix9)

# Extract longitude and latitude from PCA data frame
mixcluster <- df_pca[, c("long", "lat")]

# Add color values to the mixcluster data frame
mixcluster$color <- mix9

# Write the result to a CSV file
write.csv(mixcluster, "rzone3_core.csv")


####################################################################################################################################
####core Seed zone migration
 gf <- gfsnp_final
 write.csv(df_pca, file="df_pca.csv", row.names=FALSE)
 vars <- most_important
 df_5026 <- read.csv("test_115000_ssp126_50_mean+43pop.csv")
 df_5026 <- df_5026[,2:9] 
names(df_5026)[1] <- "Longitude"
names(df_5026)[2] <- "Latitude"
df_5085 <- read.csv("test_115000_ssp585_50_mean+43pop.csv")
df_5085 <- df_5085[,2:9]
names(df_5085)[1] <- "Longitude"
names(df_5085)[2] <- "Latitude"
library(gradientForest)
##50-126
Trns_grid_5026 <- cbind(df_5026[,c("Longitude","Latitude")], predict(gf,df_5026[,vars]))
write.csv(Trns_grid_5026, file="Trns_grid_5026.csv", row.names=FALSE)
##50-585
Trns_grid_5085 <- cbind(df_5085[,c("Longitude","Latitude")], predict(gf,df_5085[,vars]))
write.csv(Trns_grid_5085, file="Trns_grid_5085.csv", row.names=FALSE)
#install.packages("FNN")
library(FNN)
# Convert data frames to matrices for FNN
pca_matrix <- as.matrix(df_pca[, vars])
future_matrix <- as.matrix(Trns_grid_5085[, vars])
#future_matrix26 <- as.matrix(Trns_grid_5026[, vars]) #2050 ssp126

# Finding the nearest neighbor
knn_result <- get.knnx(data = future_matrix, query = pca_matrix, k = 1)
#knn_result <- get.knnx(data = future_matrix26, query = pca_matrix, k = 1)

# Extract minimum distances and their indices
min_distances <- knn_result$nn.dist
min_indices <- knn_result$nn.index
# Create a dataframe with results
#SSP585
results_df <- data.frame(
  Original_Longitude = df_pca$long,
  Original_Latitude = df_pca$lat,
  Min_Distance = min_distances,
  Future_Longitude = Trns_grid_5085$Longitude[min_indices],
  Future_Latitude = Trns_grid_5085$Latitude[min_indices]
)

write.csv(results_df, "Seed_Zone_Shifts_under_climate_change_2050_ssp585_KNN_method_extended.csv", row.names = FALSE)


#ssp126
future_matrix26 <- as.matrix(Trns_grid_5026[, vars]) #2050 ssp126
knn_result <- get.knnx(data = future_matrix26, query = pca_matrix, k = 1)

# Extract minimum distances and their indices
min_distances <- knn_result$nn.dist
min_indices <- knn_result$nn.index
results26_df <- data.frame(
  Original_Longitude = df_pca$long,
  Original_Latitude = df_pca$lat,
  Min_Distance = min_distances,
  Future_Longitude = Trns_grid_5026$Longitude[min_indices],
  Future_Latitude = Trns_grid_5026$Latitude[min_indices]
)

write.csv(results26_df, "Seed_Zone_Shifts_under_climate_change_2050_ssp126_KNN_method_extended.csv", row.names = FALSE)

########## Plotting #########
library(readxl)
library(dplyr)
library(ggplot2)
library(geosphere)
library(rnaturalearth)
library(mapdata)
library(purrr)
library(sf)
library(raster)

#ssp585_50
two_zones_df <- read.csv("Seed_Zone_Shifts_under_climate_change_2050_ssp585_KNN_method_extended.csv")
centroids <- two_zones_df %>%
  group_by(Seed_Zone) %>%
  summarise(
    Orig_Centroid_Lon = mean(Original_Longitude),
    Orig_Centroid_Lat = mean(Original_Latitude),
    Future_Centroid_Lon = mean(Future_Longitude),
    Future_Centroid_Lat = mean(Future_Latitude),
    mean_offset = mean(Min_Distance)
  )
centroids$Distance_km <- mapply(function(lon1, lat1, lon2, lat2) {
  distHaversine(c(lon1, lat1), c(lon2, lat2)) / 1000  # convert from meters to kilometers
}, centroids$Orig_Centroid_Lon, centroids$Orig_Centroid_Lat, centroids$Future_Centroid_Lon, centroids$Future_Centroid_Lat)
### The mapply function is used to apply a given function (in this case, an anonymous function) to corresponding elements from multiple lists or vectors. It takes four parameters: Orig_Centroid_Lon, Orig_Centroid_Lat, Future_Centroid_Lon, and Future_Centroid_Lat, which represent the longitude and latitude coordinates of the starting and ending points.
## The anonymous function uses the distHaversine function to calculate the great-circle distance between the two latitude-longitude pairs (in meters), then divides the result by 1000 to convert the distance to kilometers.
## centroids$Distance_km is assigned the great-circle distance (in kilometers) for each pair of starting and ending points.

arrow_data_50585 <- centroids %>%
  transmute(
    x = Orig_Centroid_Lon, y = Orig_Centroid_Lat,
    xend = Future_Centroid_Lon, yend = Future_Centroid_Lat,
    Seed_Zone)
## The centroids %>% pipe operator passes the centroids data frame to the next function.
## The transmute function selects and renames columns in the data frame. In this case, it selects four columns: Orig_Centroid_Lon renamed to x, Orig_Centroid_Lat renamed to y, Future_Centroid_Lon renamed to xend, and Future_Centroid_Lat renamed to yend, and includes the Seed_Zone column.
write.csv(arrow_data_50585,"arrow_data_50585.csv")

# Get geographic data for China
china <- ne_countries(country = "china", returnclass = "sf")
# Get provincial boundary data for China
china_provinces <- ne_states(country = "china", returnclass = "sf")
arrow_data <- arrow_data_50585 
#"coral1", "olivedrab3", "goldenrod1"

ggplot() +
  geom_sf(data = china, fill = "grey90", color = "black") +  # Draw the map of China with a grey background
  geom_sf(data = china_provinces, fill = NA, color = "black", lwd = 0.5) +  # Draw provincial boundaries
  geom_point(aes(x = two_zones_df$Future_Longitude, y = two_zones_df$Future_Latitude, col = as.factor(two_zones_df$Seed_Zone)), pch = ".", cex = 0.01) +
  scale_color_manual(values = c("1" = "coral1", "2" = "olivedrab3", "3" = "goldenrod1")) +
  geom_point(aes(x = centroids[1,]$Orig_Centroid_Lon, y = centroids[1,]$Orig_Centroid_Lat), 
             pch = 21, size = 7, fill = "coral1", color = "black", stroke = 1.2) +
  geom_point(aes(x = centroids[2,]$Orig_Centroid_Lon, y = centroids[2,]$Orig_Centroid_Lat), 
             pch = 21, size = 7, fill = "olivedrab3", color = "black", stroke = 1.2) +
  geom_point(aes(x = centroids[3,]$Orig_Centroid_Lon, y = centroids[3,]$Orig_Centroid_Lat), 
             pch = 21, size = 7, fill = "goldenrod1", color = "black", stroke = 1.2) +
  geom_point(aes(x = centroids[1,]$Future_Centroid_Lon, y = centroids[1,]$Future_Centroid_Lat), 
             pch = 21, size = 7, fill = "white", color = "black", stroke = 1.2) +
  geom_point(aes(x = centroids[2,]$Future_Centroid_Lon, y = centroids[2,]$Future_Centroid_Lat), 
             pch = 21, size = 7, fill = "white", color = "black", stroke = 1.2) +
  geom_point(aes(x = centroids[3,]$Future_Centroid_Lon, y = centroids[3,]$Future_Centroid_Lat), 
             pch = 21, size = 7, fill = "white", color = "black", stroke = 1.2) +		 
  geom_point(aes(x = centroids$Orig_Centroid_Lon, y = centroids$Orig_Centroid_Lat), pch = c("1", "2", "3"), size = 5, col = "white") +
  geom_point(aes(x = centroids[1,]$Future_Centroid_Lon, y = centroids[1,]$Future_Centroid_Lat), pch = c("1"), size = 5, col = "coral1") +
  geom_point(aes(x = centroids[2,]$Future_Centroid_Lon, y = centroids[2,]$Future_Centroid_Lat), pch = c("2"), size = 5, col = "olivedrab3") +
  geom_point(aes(x = centroids[3,]$Future_Centroid_Lon, y = centroids[3,]$Future_Centroid_Lat), pch = c("3"), size = 5, col = "goldenrod1") +
  geom_segment(data = arrow_data, aes(x = x + 0.53, y = y + 0.23, xend = (xend - 0.5), yend = (yend - 0.5)), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black", size = 0.8) +
  coord_sf(xlim = c(85, 115), ylim = c(20, 45), expand = FALSE) +
  xlab("\nLongitude") + ylab("Latitude\n") +
  theme_bw(base_size = 11) +
  theme(plot.background = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 14), panel.grid.major = element_blank(), legend.position = "none")

#ssp126_50
two_zones_df <- read.csv("Seed_Zone_Shifts_under_climate_change_2050_ssp126_KNN_method_extended.csv")

# Group by Seed_Zone and calculate centroids and distances for each seed zone
centroids <- two_zones_df %>%
  group_by(Seed_Zone) %>%
  summarise(
    Orig_Centroid_Lon = mean(Original_Longitude),  # Calculate the average of the original longitude
    Orig_Centroid_Lat = mean(Original_Latitude),   # Calculate the average of the original latitude
    Future_Centroid_Lon = mean(Future_Longitude),  # Calculate the average of the future longitude
    Future_Centroid_Lat = mean(Future_Latitude),   # Calculate the average of the future latitude
    mean_offset = mean(Min_Distance)               # Calculate the average of the minimum distance
  )

# Calculate the distance between the original and future centroids (in kilometers)
centroids$Distance_km <- mapply(function(lon1, lat1, lon2, lat2) {
  distHaversine(c(lon1, lat1), c(lon2, lat2)) / 1000  # Use the Haversine formula to calculate distance, and convert to kilometers
}, centroids$Orig_Centroid_Lon, centroids$Orig_Centroid_Lat, centroids$Future_Centroid_Lon, centroids$Future_Centroid_Lat)

# Create a new data frame with data required for arrow plotting
arrow_data_50126 <- centroids %>%
  transmute(
    x = Orig_Centroid_Lon,            # Original longitude
    y = Orig_Centroid_Lat,            # Original latitude
    xend = Future_Centroid_Lon,       # Future longitude
    yend = Future_Centroid_Lat,       # Future latitude
    Seed_Zone                         # Seed zone information
)

write.csv(arrow_data_50126, "arrow_data_50126.csv")
arrow_data <- arrow_data_50126

ggplot() +
   geom_sf(data = china, fill = "grey90", color = "black") +  # Draw the background of China's map
   geom_sf(data = china_provinces, fill = NA, color = "black", lwd = 0.5) +  # Draw provincial boundaries
   geom_point(aes(x = two_zones_df$Future_Longitude, y = two_zones_df$Future_Latitude, col = as.factor(two_zones_df$Seed_Zone)), pch = ".", cex = 0.01) +
   scale_color_manual(values = c("1" = "coral1", "2" = "olivedrab3", "3" = "goldenrod1")) +
   geom_point(aes(x = centroids[1,]$Orig_Centroid_Lon, y = centroids[1,]$Orig_Centroid_Lat), 
             pch = 21, size = 7, fill = "coral1", color = "black", stroke = 1.2) +
   geom_point(aes(x = centroids[2,]$Orig_Centroid_Lon, y = centroids[2,]$Orig_Centroid_Lat), 
             pch = 21, size = 7, fill = "olivedrab3", color = "black", stroke = 1.2) +
   geom_point(aes(x = centroids[3,]$Orig_Centroid_Lon, y = centroids[3,]$Orig_Centroid_Lat), 
             pch = 21, size = 7, fill = "goldenrod1", color = "black", stroke = 1.2) +
   geom_point(aes(x = centroids[1,]$Future_Centroid_Lon, y = centroids[1,]$Future_Centroid_Lat), 
             pch = 21, size = 7, fill = "white", color = "black", stroke = 1.2) +
   geom_point(aes(x = centroids[2,]$Future_Centroid_Lon, y = centroids[2,]$Future_Centroid_Lat), 
             pch = 21, size = 7, fill = "white", color = "black", stroke = 1.2) +
   geom_point(aes(x = centroids[3,]$Future_Centroid_Lon, y = centroids[3,]$Future_Centroid_Lat), 
             pch = 21, size = 7, fill = "white", color = "black", stroke = 1.2) +		 
   geom_point(aes(x = centroids$Orig_Centroid_Lon, y = centroids$Orig_Centroid_Lat), pch = c("1", "2", "3"), size = 5, col = "white") +
   geom_point(aes(x = centroids[1,]$Future_Centroid_Lon, y = centroids[1,]$Future_Centroid_Lat), pch = c("1"), size = 5, col = "coral1") +
   geom_point(aes(x = centroids[2,]$Future_Centroid_Lon, y = centroids[2,]$Future_Centroid_Lat), pch = c("2"), size = 5, col = "olivedrab3") +
   geom_point(aes(x = centroids[3,]$Future_Centroid_Lon, y = centroids[3,]$Future_Centroid_Lat), pch = c("3"), size = 5, col = "goldenrod1") +
   geom_segment(data = arrow_data, aes(x = x + 0.53, y = y + 0.23, xend = (xend - 0.5), yend = (yend - 0.5)), 
                arrow = arrow(length = unit(0.2, "cm")), color = "black", size = 0.8) +
   coord_sf(xlim = c(85, 115), ylim = c(20, 45), expand = FALSE) +
   xlab("\nLongitude") + ylab("Latitude\n") +
   theme_bw(base_size = 11) +
   theme(plot.background = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 14), 
         panel.grid.major = element_blank(), legend.position = "none")

######90_126、585
df_9026 <- read.csv("test_115000_ssp126_90_mean+43pop.csv")
df_9026 <- df_9026[,2:9] 
names(df_9026)[1] <- "Longitude"
names(df_9026)[2] <- "Latitude"
df_9085 <- read.csv("test_115000_ssp585_90_mean+43pop.csv")
df_9085 <- df_9085[,2:9]
names(df_9085)[1] <- "Longitude"
names(df_9085)[2] <- "Latitude"
library(gradientForest)
##90-126
Trns_grid_9026 <- cbind(df_9026[,c("Longitude","Latitude")], predict(gf,df_9026[,vars]))
write.csv(Trns_grid_9026, file="Trns_grid_9026.csv", row.names=FALSE)
##90-585
Trns_grid_9085 <- cbind(df_9085[,c("Longitude","Latitude")], predict(gf,df_9085[,vars]))
write.csv(Trns_grid_9085, file="Trns_grid_9085.csv", row.names=FALSE)

#install.packages("FNN")
library(FNN)
# Convert data frames to matrices for FNN
pca_matrix <- as.matrix(df_pca[, vars])
future_matrix <- as.matrix(Trns_grid_9085[, vars])
future_matrix26 <- as.matrix(Trns_grid_9026[, vars]) #2090 ssp126

# Finding the nearest neighbor
knn_result <- get.knnx(data = future_matrix, query = pca_matrix, k = 1)
#knn_result <- get.knnx(data = future_matrix26, query = pca_matrix, k = 1)

# Extract minimum distances and their indices
min_distances <- knn_result$nn.dist
min_indices <- knn_result$nn.index

results_df <- data.frame(
  Original_Longitude = df_pca$long,
  Original_Latitude = df_pca$lat,
  Min_Distance = min_distances,
  Future_Longitude = Trns_grid_9085$Longitude[min_indices],
  Future_Latitude = Trns_grid_9085$Latitude[min_indices]
)

write.csv(results_df, "Seed_Zone_Shifts_under_climate_change_2090_ssp585_KNN_method_extended.csv", row.names = FALSE)

#ssp126

knn_result <- get.knnx(data = future_matrix26, query = pca_matrix, k = 1)

# Extract minimum distances and their indices
min_distances <- knn_result$nn.dist
min_indices <- knn_result$nn.index
results26_df <- data.frame(
  Original_Longitude = df_pca$long,
  Original_Latitude = df_pca$lat,
  Min_Distance = min_distances,
  Future_Longitude = Trns_grid_9026$Longitude[min_indices],
  Future_Latitude = Trns_grid_9026$Latitude[min_indices]
)

write.csv(results26_df, "Seed_Zone_Shifts_under_climate_change_2090_ssp126_KNN_method_extended.csv", row.names = FALSE)

#ssp126_90

two_zones_df <- read.csv("Seed_Zone_Shifts_under_climate_change_2090_ssp126_KNN_method_extended.csv")

# Group by Seed_Zone and calculate the centroid and distance for each seed zone
centroids <- two_zones_df %>%
  group_by(Seed_Zone) %>%
  summarise(
    Orig_Centroid_Lon = mean(Original_Longitude),  # Calculate the mean of the original longitude
    Orig_Centroid_Lat = mean(Original_Latitude),   # Calculate the mean of the original latitude
    Future_Centroid_Lon = mean(Future_Longitude),  # Calculate the mean of the future longitude
    Future_Centroid_Lat = mean(Future_Latitude),   # Calculate the mean of the future latitude
    mean_offset = mean(Min_Distance)               # Calculate the mean of the minimum distance
  )

# Calculate the distance between the original and future centroids (in kilometers)
centroids$Distance_km <- mapply(function(lon1, lat1, lon2, lat2) {
  distHaversine(c(lon1, lat1), c(lon2, lat2)) / 1000  # Use the Haversine formula to calculate distance and convert to kilometers
}, centroids$Orig_Centroid_Lon, centroids$Orig_Centroid_Lat, centroids$Future_Centroid_Lon, centroids$Future_Centroid_Lat)

# Create a new data frame with data required for arrow plotting
arrow_data_90126 <- centroids %>%
  transmute(
    x = Orig_Centroid_Lon,            # Original longitude
    y = Orig_Centroid_Lat,            # Original latitude
    xend = Future_Centroid_Lon,       # Future longitude
    yend = Future_Centroid_Lat,       # Future latitude
    Seed_Zone                         # Seed zone information
  )

write.csv(arrow_data_90126, "arrow_data_90126.csv")
arrow_data <- arrow_data_90126

# Create the plot
ggplot() +
   geom_sf(data = china, fill = "grey90", color = "black") +  # Draw the background of China's map
   geom_sf(data = china_provinces, fill = NA, color = "black", lwd = 0.5) +  # Draw provincial boundaries
   geom_point(aes(x = two_zones_df$Future_Longitude, y = two_zones_df$Future_Latitude, col = as.factor(two_zones_df$Seed_Zone)), pch = ".", cex = 0.01) +
   scale_color_manual(values = c("1" = "coral1", "2" = "olivedrab3", "3" = "goldenrod1")) +
   geom_point(aes(x = centroids[1,]$Orig_Centroid_Lon, y = centroids[1,]$Orig_Centroid_Lat), 
              pch = 21, size = 7, fill = "coral1", color = "black", stroke = 1.2) +
   geom_point(aes(x = centroids[2,]$Orig_Centroid_Lon, y = centroids[2,]$Orig_Centroid_Lat), 
              pch = 21, size = 7, fill = "olivedrab3", color = "black", stroke = 1.2) +
   geom_point(aes(x = centroids[3,]$Orig_Centroid_Lon, y = centroids[3,]$Orig_Centroid_Lat), 
              pch = 21, size = 7, fill = "goldenrod1", color = "black", stroke = 1.2) +
   geom_point(aes(x = centroids[1,]$Future_Centroid_Lon, y = centroids[1,]$Future_Centroid_Lat), 
              pch = 21, size = 7, fill = "white", color = "black", stroke = 1.2) +
   geom_point(aes(x = centroids[2,]$Future_Centroid_Lon, y = centroids[2,]$Future_Centroid_Lat), 
              pch = 21, size = 7, fill = "white", color = "black", stroke = 1.2) +
   geom_point(aes(x = centroids[3,]$Future_Centroid_Lon, y = centroids[3,]$Future_Centroid_Lat), 
              pch = 21, size = 7, fill = "white", color = "black", stroke = 1.2) +		 
   geom_point(aes(x = centroids$Orig_Centroid_Lon, y = centroids$Orig_Centroid_Lat), pch = c("1", "2", "3"), size = 5, col = "white") +
   geom_point(aes(x = centroids[1,]$Future_Centroid_Lon, y = centroids[1,]$Future_Centroid_Lat), pch = c("1"), size = 5, col = "coral1") +
   geom_point(aes(x = centroids[2,]$Future_Centroid_Lon, y = centroids[2,]$Future_Centroid_Lat), pch = c("2"), size = 5, col = "olivedrab3") +
   geom_point(aes(x = centroids[3,]$Future_Centroid_Lon, y = centroids[3,]$Future_Centroid_Lat), pch = c("3"), size = 5, col = "goldenrod1") +
   geom_segment(data = arrow_data, aes(x = x + 0.53, y = y + 0.23, xend = (xend - 0.5), yend = (yend - 0.5)), 
                arrow = arrow(length = unit(0.2, "cm")), color = "black", size = 0.8) +
   coord_sf(xlim = c(85, 115), ylim = c(20, 45), expand = FALSE) +
   xlab("\nLongitude") + ylab("Latitude\n") +
   theme_bw(base_size = 11) +
   theme(plot.background = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 14), 
         panel.grid.major = element_blank(), legend.position = "none")

 ##58590

two_zones_df <- read.csv("Seed_Zone_Shifts_under_climate_change_2090_ssp585_KNN_method_extended.csv")

# Group by Seed_Zone and calculate the centroid and distance for each seed zone
centroids <- two_zones_df %>%
  group_by(Seed_Zone) %>%
  summarise(
    Orig_Centroid_Lon = mean(Original_Longitude),  # Calculate the average of the original longitude
    Orig_Centroid_Lat = mean(Original_Latitude),   # Calculate the average of the original latitude
    Future_Centroid_Lon = mean(Future_Longitude),  # Calculate the average of the future longitude
    Future_Centroid_Lat = mean(Future_Latitude),   # Calculate the average of the future latitude
    mean_offset = mean(Min_Distance)               # Calculate the average of the minimum distance
  )

# Calculate the distance between the original and future centroids (in kilometers)
centroids$Distance_km <- mapply(function(lon1, lat1, lon2, lat2) {
  distHaversine(c(lon1, lat1), c(lon2, lat2)) / 1000  # Use the Haversine formula to calculate distance and convert to kilometers
}, centroids$Orig_Centroid_Lon, centroids$Orig_Centroid_Lat, centroids$Future_Centroid_Lon, centroids$Future_Centroid_Lat)

# Create a new data frame with data required for arrow plotting
arrow_data_90585 <- centroids %>%
  transmute(
    x = Orig_Centroid_Lon,            # Original longitude
    y = Orig_Centroid_Lat,            # Original latitude
    xend = Future_Centroid_Lon,       # Future longitude
    yend = Future_Centroid_Lat,       # Future latitude
    Seed_Zone                         # Seed zone information
  )

write.csv(arrow_data_90585, "arrow_data_90585.csv")
arrow_data <- arrow_data_90585

# Create the plot
ggplot() +
   geom_sf(data = china, fill = "grey90", color = "black") +  # Draw the background of China's map
   geom_sf(data = china_provinces, fill = NA, color = "black", lwd = 0.5) +  # Draw provincial boundaries
   geom_point(aes(x = two_zones_df$Future_Longitude, y = two_zones_df$Future_Latitude, col = as.factor(two_zones_df$Seed_Zone)), pch = ".", cex = 0.01) +
   scale_color_manual(values = c("1" = "coral1", "2" = "olivedrab3", "3" = "goldenrod1")) +
   geom_point(aes(x = centroids[1,]$Orig_Centroid_Lon, y = centroids[1,]$Orig_Centroid_Lat), 
              pch = 21, size = 7, fill = "coral1", color = "black", stroke = 1.2) +
   geom_point(aes(x = centroids[2,]$Orig_Centroid_Lon, y = centroids[2,]$Orig_Centroid_Lat), 
              pch = 21, size = 7, fill = "olivedrab3", color = "black", stroke = 1.2) +
   geom_point(aes(x = centroids[3,]$Orig_Centroid_Lon, y = centroids[3,]$Orig_Centroid_Lat), 
              pch = 21, size = 7, fill = "goldenrod1", color = "black", stroke = 1.2) +
   geom_point(aes(x = centroids[1,]$Future_Centroid_Lon, y = centroids[1,]$Future_Centroid_Lat), 
              pch = 21, size = 7, fill = "white", color = "black", stroke = 1.2) +
   geom_point(aes(x = centroids[2,]$Future_Centroid_Lon, y = centroids[2,]$Future_Centroid_Lat), 
              pch = 21, size = 7, fill = "white", color = "black", stroke = 1.2) +
   geom_point(aes(x = centroids[3,]$Future_Centroid_Lon, y = centroids[3,]$Future_Centroid_Lat), 
              pch = 21, size = 7, fill = "white", color = "black", stroke = 1.2) +		 
   geom_point(aes(x = centroids$Orig_Centroid_Lon, y = centroids$Orig_Centroid_Lat), pch = c("1", "2", "3"), size = 5, col = "white") +
   geom_point(aes(x = centroids[1,]$Future_Centroid_Lon, y = centroids[1,]$Future_Centroid_Lat), pch = c("1"), size = 5, col = "coral1") +
   geom_point(aes(x = centroids[2,]$Future_Centroid_Lon, y = centroids[2,]$Future_Centroid_Lat), pch = c("2"), size = 5, col = "olivedrab3") +
   geom_point(aes(x = centroids[3,]$Future_Centroid_Lon, y = centroids[3,]$Future_Centroid_Lat), pch = c("3"), size = 5, col = "goldenrod1") +
   geom_segment(data = arrow_data, aes(x = x + 0.53, y = y + 0.23, xend = (xend - 0.5), yend = (yend - 0.5)), 
                arrow = arrow(length = unit(0.2, "cm")), color = "black", size = 0.8) +
   coord_sf(xlim = c(85, 115), ylim = c(20, 45), expand = FALSE) +
   xlab("\nLongitude") + ylab("Latitude\n") +
   theme_bw(base_size = 11) +
   theme(plot.background = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 14), 
         panel.grid.major = element_blank(), legend.position = "none")



# Sample
plink --bfile core_adapt --export vcf --out core_adapt --allow-extra-chr
vcftools --vcf core_adapt.vcf --keep corezone1.txt --out corezone1 --recode
vcftools --vcf core_adapt.vcf --keep corezone2.txt --out corezone2 --recode
vcftools --vcf core_adapt.vcf --keep corezone3.txt --out corezone3 --recode

# Run the script to calculate the optimal sample size
nohup python capturing_diversity3.py --vcf corezone1.recode.vcf > zone1.log 2>&1 &
nohup python capturing_diversity3.py --vcf corezone2.recode.vcf > zone2.log 2>&1 &
nohup python capturing_diversity3.py --vcf corezone3.recode.vcf > zone3.log 2>&1 &

# Set different thresholds with -d
nohup python capturing_diversity3.py --vcf corezone1.recode.vcf -d 0.95 > zone1_0.95.log 2>&1 &
nohup python capturing_diversity3.py --vcf corezone2.recode.vcf -d 0.95 > zone2_0.95.log 2>&1 &
nohup python capturing_diversity3.py --vcf corezone3.recode.vcf -d 0.95 > zone3_0.95.log 2>&1 &

nohup python capturing_diversity3.py --vcf corezone1.recode.vcf -d 0.99 > zone1_0.99.log 2>&1 &
nohup python capturing_diversity3.py --vcf corezone2.recode.vcf -d 0.99 > zone2_0.99.log 2>&1 &
nohup python capturing_diversity3.py --vcf corezone3.recode.vcf -d 0.99 > zone3_0.99.log 2>&1 &

# Process the pan-adaptation dataset
plink --bfile pan_adapt --export vcf --out pan_adapt --allow-extra-chr
vcftools --vcf pan_adapt.vcf --keep corezone1.txt --out panzone1 --recode
vcftools --vcf pan_adapt.vcf --keep corezone2.txt --out panzone2 --recode
vcftools --vcf pan_adapt.vcf --keep corezone3.txt --out panzone3 --recode

# Run the script to calculate the optimal sample size for pan-adapt dataset
nohup python capturing_diversity3.py --vcf panzone1.recode.vcf > zone1.log 2>&1 &
nohup python capturing_diversity3.py --vcf panzone2.recode.vcf > zone2.log 2>&1 &
nohup python capturing_diversity3.py --vcf panzone3.recode.vcf > zone3.log 2>&1 &

# Set different thresholds with -d for pan-adapt dataset
nohup python capturing_diversity3.py --vcf panzone1.recode.vcf -d 0.95 > zone1_0.95.log 2>&1 &
nohup python capturing_diversity3.py --vcf panzone2.recode.vcf -d 0.95 > zone2_0.95.log 2>&1 &
nohup python capturing_diversity3.py --vcf panzone3.recode.vcf -d 0.95 > zone3_0.95.log 2>&1 &

nohup python capturing_diversity3.py --vcf panzone1.recode.vcf -d 0.99 > zone1_0.99.log 2>&1 &
nohup python capturing_diversity3.py --vcf panzone2.recode.vcf -d 0.99 > zone2_0.99.log 2>&1 &
nohup python capturing_diversity3.py --vcf panzone3.recode.vcf -d 0.99 > zone3_0.99.log 2>&1 &



