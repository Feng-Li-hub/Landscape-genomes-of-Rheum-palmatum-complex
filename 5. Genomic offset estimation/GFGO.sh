#### Adaptive loci AF file processing

# For each row in core_loci.csv, concatenate the first and second columns with an underscore, 
# and save the result to core_loci2.csv
awk -F',' '{print \$1"_"\$2}' core_loci.csv > core_loci2.csv

# For each row in pan_loci.csv, concatenate the first and second columns with an underscore, 
# and save the result to pan_loci2.csv
awk -F',' '{print \$1"_"\$2}' pan_loci.csv > pan_loci2.csv

library(dplyr)       
library(data.table)  

# Read the allele frequency file
allele_freq <- fread("F:\\王聪颖\\景观_2\\AF\\allele_freq_AF2.csv")
# Read the core loci file, assuming no column names and only SNP names in the first column
core_loci2 <- fread("core_loci2.csv", header = FALSE)
# Extract SNP names from the first column of core_loci2
snp_names <- core_loci2$V1

# Extract the relevant columns from allele_freq corresponding to SNP names
# Filter out the ID column and select columns based on the SNP names
columns_to_extract <- c(snp_names)  # SNP names to extract, including the first column and the columns of interest
# Select the columns from allele_freq based on the SNP names
filtered_data_core <- allele_freq %>% select(all_of(columns_to_extract))
fwrite(filtered_data_core, "core_loci2_AF.csv", row.names = FALSE)

# Read the pan loci file, assuming no column names and only SNP names in the first column
pan_loci2 <- fread("pan_loci2.csv", header = FALSE)
# Extract SNP names from the first column of pan_loci2
snp_names <- pan_loci2$V1
# Extract the relevant columns from allele_freq corresponding to SNP names
# Filter out the ID column and select columns based on the SNP names
columns_to_extract <- c(snp_names)  # SNP names to extract, including the first column and the columns of interest
# Select the columns from allele_freq based on the SNP names
filtered_data_pan <- allele_freq %>% select(all_of(columns_to_extract))
fwrite(filtered_data_pan, "pan_loci2_AF.csv", row.names = FALSE)

#############################GFGO################################
library(gradientForest)  
# Read the allele frequency data for pan (the dataset containing SNPs and environmental data)
gfData <- read.csv("pan_loci2_AF.csv")
# Extract climate variables (columns 4 to 9) from the data
envGF <- gfData[,4:9] 
# Extract SNP allele frequency data based on the column names containing "RpChr"
snp <- gfData[,grep("RpChr", colnames(gfData))] 
# Check the number of SNPs (alleles) extracted
length(snp)  # Output: 16709 SNPs

# Calculate the maximum level parameter based on the number of rows in envGF and correlation threshold
maxLevel <- log2(0.368 * nrow(envGF) / 2)  # Account for correlations, reference in ?gradientForest

# Build the gradient forest (GF) model using the environmental variables and SNP allele frequencies
gfsnp_final <- gradientForest(
  cbind(envGF, snp),               # Combine environmental and SNP data
  predictor.vars = colnames(envGF), # Specify predictor variables (environmental variables)
  response.vars = colnames(snp),    # Specify response variables (SNP data)
  ntree = 500,                      # Number of trees in the random forest model
  maxLevel = maxLevel,              # Maximum level parameter, calculated above
  trace = TRUE,                     # Print progress as the model is being built
  corr.threshold = 0.50             # Correlation threshold for variable selection
)

# Read contemporary climate data for predictions
DH_current_bio <- read.csv("test_115000_current+43pop.csv")

# Define the most important climate variables for this analysis
most_important <- c("BIO02", "BIO03", "BIO04", "BIO10", "BIO12", "BIO15") 

# Select the important environmental variables
imp.vars <- most_important

# Calculate the allele frequencies for the selected environment variables and combine with geographical coordinates
Trns_grid <- cbind(
  DH_current_bio[, c("long", "lat")],                # Add longitude and latitude from the climate data
  predict(gfsnp_final, DH_current_bio[, imp.vars])   # Predict SNP allele frequencies based on selected environmental variables
)

# Save the result (genetic variation in relation to climate) as a CSV file
write.csv(Trns_grid, "Trns_grid.csv")

# Perform PCA (Principal Component Analysis) on the selected environmental variables
PCs <- prcomp(Trns_grid[, imp.vars])

# Extract the first three principal components for plotting
a1 <- PCs$x[,1]
a2 <- PCs$x[,2]
a3 <- PCs$x[,3]

# Generate RGB values for visualization based on the PCA results
r <- a1 + a2
g <- -a2
b <- a3 + a2 - a1

# Normalize the RGB values to be in the range of 0 to 255
r <- (r - min(r)) / (max(r) - min(r)) * 255
g <- (g - min(g)) / (max(g) - min(g)) * 255
b <- (b - min(b)) / (max(b) - min(b)) * 255

# Create a PDF for the map showing the current distribution of genetic variation
pdf(file = "current distribution of genetic variation.pdf")
plot(
  Trns_grid[, c("long", "lat")],           # Plot longitude and latitude
  pch = ".",                               # Set plot character as dot
  cex = 3,                                 # Set dot size
  asp = 1,                                 # Set aspect ratio to 1 (equal scaling of x and y axes)
  col = rgb(r, g, b, max = 255)            # Use RGB colors for points based on PCA results
)
dev.off()  # Close the PDF file after plotting

# Output PCA biplot for the selected important environmental variables
nvs <- dim(PCs$rotation)[1]  # Number of variables in the rotation matrix
vec <- c("BIO02", "BIO03", "BIO04", "BIO10", "BIO12", "BIO15")  # List of important variables
lv <- length(vec)  # Number of selected variables

# Set scaling factor for PCA vectors
scal <- 15

# Define the x and y ranges for the PCA plot
xrng <- range(PCs$x[, 1], PCs$rotation[, 1] / scal) * 1.1
yrng <- range(PCs$x[, 2], PCs$rotation[, 2] / scal) * 1.1

# Create a PDF for the PCA biplot
pdf(file = "PCA biplot.pdf")
plot(
  (PCs$x[, 1:2]),                           # Plot the first two PCA components
  xlim = xrng,                              # Set x-axis limits
  ylim = yrng,                              # Set y-axis limits
  pch = ".",                                # Plot character as dot
  cex = 4,                                  # Set dot size
  col = rgb(r, g, b, max = 255),            # Use RGB colors based on PCA results
  asp = 1                                    # Set aspect ratio to 1
)

# Add arrows for the important variables (BIO variables) in the PCA plot
arrows(
  rep(0, lv),                               # Starting point of the arrows (origin)
  rep(0, lv),                               # Starting point of the arrows (origin)
  PCs$rotation[vec, 1] / scal,              # x-coordinate of the arrow tips
  PCs$rotation[vec, 2] / scal,              # y-coordinate of the arrow tips
  length = 0.0625,                          # Length of the arrowheads
  lwd = 1                                    # Line width for the arrows
)

# Add labels to the arrows with some jitter for better visibility
jit <- 0.004
text(
  PCs$rotation[vec, 1] / scal + jit * sign(PCs$rotation[vec, 1]),
  PCs$rotation[vec, 2] / scal + jit * sign(PCs$rotation[vec, 2]),
  labels = vec,                             # Labels for the arrows
  cex = 1,                                  # Font size for the labels
  family = "serif"                          # Font family
)

dev.off()  # Close the PDF file after plotting

################# Calculate genetic components under future important variables #################
### ssp126_50
fut_5026 <- read.csv("test_115000_ssp126_50_mean+43pop.csv") 
# Predict allele frequencies based on the future climate data and selected important variables
Trns_grid_5026 <- cbind(fut_5026[,c("long","lat")], predict(gfsnp_final,fut_5026[,imp.vars]))
write.csv(Trns_grid_5026,"Trns_grid_5026.csv")
## Calculate GO (Genetic Offset) values using Euclidean distance
euclidean <- function(a, b) sqrt(sum((a - b)^2))  # Define the Euclidean distance function for genetic offset
df <- data.frame()  # Create an empty data frame to store GO values

# Loop over each row in the future dataset
for (i in 1:nrow(fut_5026)) {
    enc_dist <- euclidean(Trns_grid[i, imp.vars], Trns_grid_5026[i, imp.vars])  # Calculate genetic offset
    df <- rbind(df, enc_dist)  # Add the calculated distance to the data frame
    names(df) <- "GO_5026"  # Name the column as GO_5026
}

write.csv(df, "DH_GO_5026.csv")  # Save the calculated genetic offset values

DH_GO_5026_final <- cbind(fut_5026[,c("long", "lat")], df)  # Combine the geographic data with GO values
write.csv(DH_GO_5026_final, "DH_GO_5026_final.csv")  # Save the final data to a CSV
## After saving the result, remove the first 43 sample points from the dataset, ensure that each latitude appears at least twice, 
## and save the file as "test_DH_GO_5026_final.csv".

# Convert points of average value to raster
library(dplyr)  ## Load dplyr for data manipulation
library(raster)  ## Load raster package for converting data to raster format

test_DH_GO_5026 <- read.csv("test_DH_GO_5026_final.csv")  # Read the final dataset with latitude and longitude
test <- test_DH_GO_5026[,-1]  ## Remove the first column (cell number)

# Convert the data to a raster object
DH5026_raster <- test %>% rasterFromXYZ(crs="+proj=longlat +datum=WGS84 +no_defs")

pdf(file="DHGO5026.pdf")  
plot(DH5026_raster) 
dev.off() 
writeRaster(DH5026_raster, "DH5026_raster.asc", format="ascii") 

### ssp126_90
fut_9026 <- read.csv("test_115000_ssp126_90_mean+43pop.csv") 
Trns_grid_9026 <- cbind(fut_9026[,c("long", "lat")], predict(gfsnp_final, fut_9026[,imp.vars]))
write.csv(Trns_grid_9026, "Trns_grid_9026.csv")

## Calculate GO (Genetic Offset) values using Euclidean distance
df <- data.frame()  # Create an empty data frame to store GO values
# Loop over each row in the future dataset
for (i in 1:nrow(fut_9026)) {
    enc_dist <- euclidean(Trns_grid[i, imp.vars], Trns_grid_9026[i, imp.vars])  # Calculate genetic offset
    df <- rbind(df, enc_dist)  # Add the calculated distance to the data frame
    names(df) <- "GO_9026"  # Name the column as GO_9026
}
write.csv(df, "DH_GO_9026.csv") 

DH_GO_9026_final <- cbind(fut_9026[,c("long", "lat")], df)  # Combine geographic data with GO values
write.csv(DH_GO_9026_final, "DH_GO_9026_final.csv")  # Save the final data to a CSV
## After saving the result, remove the first 43 sample points from the dataset, ensure that each latitude appears at least twice,
## and save the file as "test_DH_GO_9026_final.csv".

# Convert points of average value to raster
test_DH_GO_9026 <- read.csv("test_DH_GO_9026_final.csv")  # Read the final dataset with latitude and longitude
test <- test_DH_GO_9026[,-1]  ## Remove the first column (cell number)

# Convert the data to a raster object
DH9026_raster <- test %>% rasterFromXYZ(crs="+proj=longlat +datum=WGS84 +no_defs")

pdf(file="DHGO9026.pdf") 
plot(DH9026_raster)  
dev.off() 
writeRaster(DH9026_raster, "DH9026_raster.asc", format="ascii")  

##### ssp585_50
fut_5085 <- read.csv("test_115000_ssp585_50_mean+43pop.csv") 
Trns_grid_5085 <- cbind(fut_5085[,c("long", "lat")], predict(gfsnp_final, fut_5085[,imp.vars]))
write.csv(Trns_grid_5085, "Trns_grid_5085.csv")
## Calculate GO (Genetic Offset) values using Euclidean distance
df <- data.frame()  # Create an empty data frame to store GO values

# Loop over each row in the future dataset
for (i in 1:nrow(fut_5085)) {
    enc_dist <- euclidean(Trns_grid[i, imp.vars], Trns_grid_5085[i, imp.vars])  # Calculate genetic offset
    df <- rbind(df, enc_dist)  # Add the calculated distance to the data frame
    names(df) <- "GO_5085"  # Name the column as GO_5085
}

write.csv(df, "DH_GO_5085.csv")  

DH_GO_5085_final <- cbind(fut_5085[,c("long", "lat")], df)  # Combine geographic data with GO values
write.csv(DH_GO_5085_final, "DH_GO_5085_final.csv") 
## After saving the result, remove the first 43 sample points from the dataset, ensure that each latitude appears at least twice, 
## and save the file as "test_DH_GO_5085_final.csv".

# Convert points of average value to raster
test_DH_GO_5085 <- read.csv("test_DH_GO_5085_final.csv")  # Read the final dataset with latitude and longitude
test <- test_DH_GO_5085[,-1]  ## Remove the first column (cell number)

# Convert the data to a raster object
DH5085_raster <- test %>% rasterFromXYZ(crs="+proj=longlat +datum=WGS84 +no_defs")

pdf(file="DHGO5085.pdf")  
plot(DH5085_raster)  
dev.off() 
writeRaster(DH5085_raster, "DH5085_raster.asc", format="ascii")  

### ssp585_90
fut_9085 <- read.csv("test_115000_ssp585_90_mean+43pop.csv") 
Trns_grid_9085 <- cbind(fut_9085[,c("long", "lat")], predict(gfsnp_final, fut_9085[,imp.vars]))
write.csv(Trns_grid_9085, "Trns_grid_9085.csv")
## Calculate GO (Genetic Offset) values using Euclidean distance
df <- data.frame()  # Create an empty data frame to store GO values

for (i in 1:nrow(fut_9085)){
  # Calculate the Euclidean distance between the predicted allele frequencies for the current row in 'Trns_grid' and 'Trns_grid_9085'
  enc_dist <- euclidean(Trns_grid[i, imp.vars], Trns_grid_9085[i, imp.vars])
  
  # Append the calculated genetic offset distance to the data frame 'df'
  df <- rbind(df, enc_dist)
  
  # Assign the column name 'GO_9085' to the data frame 'df'
  names(df) <- "GO_9085"
}
write.csv(df, "DH_GO_9085.csv")
DH_GO_9085_final <- cbind(fut_9085[, c("long", "lat")], df)
write.csv(DH_GO_9085_final, "DH_GO_9085_final.csv")
## After outputting the results, remove the first 43 sample points from the dataset,
## check that the latitude column has at least two occurrences for each latitude value,
## and then save the file as 'test_DH_GO_9085_final.csv'.


library(dplyr)  
test_DH_GO_9085 <- read.csv("test_DH_GO_9085_final.csv")
test <- test_DH_GO_9085[, -1] 

# Convert the dataset to a raster object using the 'rasterFromXYZ' function from the raster package
DH9085_raster <- test %>% rasterFromXYZ(crs = "+proj=longlat +datum=WGS84 +no_defs")


pdf(file = "DHGO9085.pdf")
plot(DH9085_raster)
dev.off()
writeRaster(DH9085_raster, "DH9085_raster.asc", format = "ascii")

 


###################################################################################################################
##########pan_east
library(gradientForest)
gfData <- read.csv("pan_east_AF.csv")
envGF <- gfData[, 4:9]
# Extract the SNP (allele frequency) data by finding columns with "RpChr" in their name
snp <- gfData[, grep("RpChr", colnames(gfData))]
# Check the number of SNPs extracted (should return 16709)
length(snp)
# Calculate maxLevel parameter to account for correlations (used in gradientForest)
maxLevel <- log2(0.368 * nrow(envGF) / 2)

# Build the Gradient Forest (GF) model
gfsnp_final <- gradientForest(
  cbind(envGF, snp),                # Combine the climate variables and SNP data
  predictor.vars = colnames(envGF),  # Define the predictor variables (climate variables)
  response.vars = colnames(snp),     # Define the response variables (SNP data)
  ntree = 500,                       # Number of trees to use in the forest
  maxLevel = maxLevel,               # Max level parameter to handle correlations
  trace = TRUE,                      # Show progress during model building
  corr.threshold = 0.50              # Correlation threshold to filter variables
)

# Load the current data for genetic variation analysis
DH_current_bio <- read.csv("test_115000_current+15pop.csv")

# Define the most important 6 BIO variables (climate variables)
most_important <- c("BIO02", "BIO03", "BIO04", "BIO10", "BIO12", "BIO15")

# Predict genetic variation (allele frequency) based on the most important environmental variables
imp.vars <- most_important
Trns_grid <- cbind(
  DH_current_bio[, c("long", "lat")],  # Combine longitude and latitude data
  predict(gfsnp_final, DH_current_bio[, imp.vars])  # Predict genetic variation
)
write.csv(Trns_grid, "Trns_grid.csv")

# Perform Principal Component Analysis (PCA) on the selected environmental variables
PCs <- prcomp(Trns_grid[, imp.vars])

# Extract the first three principal components (PC1, PC2, PC3)
a1 <- PCs$x[, 1]
a2 <- PCs$x[, 2]
a3 <- PCs$x[, 3]

# Create RGB values based on the PCA components for visualization
r <- a1 + a2
g <- -a2
b <- a3 + a2 - a1

# Normalize RGB values to a range of 0-255
r <- (r - min(r)) / (max(r) - min(r)) * 255
g <- (g - min(g)) / (max(g) - min(g)) * 255
b <- (b - min(b)) / (max(b) - min(b)) * 255

pdf(file = "current distribution of genetic variation.pdf")
plot(Trns_grid[, c("long", "lat")], pch = ".", cex = 3, asp = 1, col = rgb(r, g, b, max = 255))  # Plot using RGB colors
dev.off()

# Plot PCA biplot for the important environmental variables
nvs <- dim(PCs$rotation)[1]
vec <- c("BIO02", "BIO03", "BIO04", "BIO10", "BIO12", "BIO15")
lv <- length(vec)

# Scale and adjust the axis ranges for the PCA biplot
scal <- 15
xrng <- range(PCs$x[, 1], PCs$rotation[, 1] / scal) * 1.1
yrng <- range(PCs$x[, 2], PCs$rotation[, 2] / scal) * 1.1

# Save the PCA biplot to a PDF
pdf(file = "PCA biplot.pdf")
plot(
  (PCs$x[, 1:2]), 
  xlim = xrng, ylim = yrng, 
  pch = ".", cex = 4, col = rgb(r, g, b, max = 255), asp = 1
)
arrows(rep(0, lv), rep(0, lv), PCs$rotation[vec, 1] / scal, PCs$rotation[vec, 2] / scal, length = 0.0625, lwd = 1)
jit <- 0.004
text(PCs$rotation[vec, 1] / scal + jit * sign(PCs$rotation[vec, 1]), 
     PCs$rotation[vec, 2] / scal + jit * sign(PCs$rotation[vec, 2]), 
     labels = vec, cex = 1, family = "serif")
dev.off()

# Process genetic variation for future climate scenario (SSP126-50)
fut_5026 <- read.csv("test_115000_ssp126_50_mean+15pop.csv")

# Predict genetic variation for the future scenario using the GF model
Trns_grid_5026 <- cbind(fut_5026[, c("long", "lat")], predict(gfsnp_final, fut_5026[, imp.vars]))
write.csv(Trns_grid_5026, "Trns_grid_5026.csv")

# Calculate Genetic Offset (GO) values using Euclidean distance
euclidean <- function(a, b) sqrt(sum((a - b)^2))  # Function to compute Euclidean distance

# Calculate the GO values for each point
df <- data.frame()
for (i in 1:nrow(fut_5026)) {
  enc_dist <- euclidean(Trns_grid[i, imp.vars], Trns_grid_5026[i, imp.vars])
  df <- rbind(df, enc_dist)
  names(df) <- "GO_5026"
}
write.csv(df, "DH_GO_5026.csv")

# Combine the GO values with geographic coordinates and save the results
DH_GO_5026_final <- cbind(fut_5026[, c("long", "lat")], df)
write.csv(DH_GO_5026_final, "DH_GO_5026_final.csv")
# Clean up the data: remove the first 43 sample points, check for duplicate latitudes, and save the cleaned data
# Save the cleaned file as 'test_DH_GO_5026_final.csv'

# Convert the genetic offset data to raster format
library(dplyr)
library(raster)
test_DH_GO_5026 <- read.csv("test_DH_GO_5026_final.csv")
test <- test_DH_GO_5026[, -1]

# Convert the data to a raster object with specified coordinate reference system (WGS84)
DH5026_raster <- test %>% rasterFromXYZ(crs = "+proj=longlat +datum=WGS84 +no_defs")

# Plot the raster and save it as a PDF
pdf(file = "DHGO5026.pdf")
plot(DH5026_raster)
dev.off()
writeRaster(DH5026_raster, "DH5026_raster.asc", format = "ascii")


### ssp126_90
fut_9026 <- read.csv("test_115000_ssp126_90_mean+15pop.csv") 
Trns_grid_9026 <- cbind(fut_9026[, c("long", "lat")], predict(gfsnp_final, fut_9026[, imp.vars]))
# Predict genetic variation for future scenario and combine with coordinates
write.csv(Trns_grid_9026, "Trns_grid_9026.csv")

## Calculate Genetic Offset (GO) values using Euclidean distance
euclidean <- function(a, b) sqrt(sum((a - b)^2))  # Function to compute Euclidean distance
df <- data.frame()
for (i in 1:nrow(fut_9026)) {
  enc_dist <- euclidean(Trns_grid[i, imp.vars], Trns_grid_9026[i, imp.vars])
  df <- rbind(df, enc_dist)
  names(df) <- "GO_9026"
}
write.csv(df, "DH_GO_9026.csv")
DH_GO_9026_final <- cbind(fut_9026[, c("long", "lat")], df)
write.csv(DH_GO_9026_final, "DH_GO_9026_final.csv")
# Combine the geographic coordinates and GO values and save

# Convert the genetic offset data into a raster for visualization
library(dplyr)
test_DH_GO_9026 <- read.csv("test_DH_GO_9026_final.csv")
test <- test_DH_GO_9026[, -1] 
DH9026_raster <- test %>% rasterFromXYZ(crs = "+proj=longlat +datum=WGS84 +no_defs")
pdf(file = "DHGO9026.pdf")
plot(DH9026_raster)
dev.off()
writeRaster(DH9026_raster, "DH9026_raster.asc", format = "ascii")

##### ssp585_50 
fut_5085 <- read.csv("test_115000_ssp585_50_mean+15pop.csv") 
Trns_grid_5085 <- cbind(fut_5085[, c("long", "lat")], predict(gfsnp_final, fut_5085[, imp.vars]))
write.csv(Trns_grid_5085, "Trns_grid_5085.csv")
# Predict genetic variation for future scenario and combine with coordinates

## Calculate Genetic Offset (GO) values for the SSP585-50 scenario
df <- data.frame()
for (i in 1:nrow(fut_5085)) {
  enc_dist <- euclidean(Trns_grid[i, imp.vars], Trns_grid_5085[i, imp.vars])
  df <- rbind(df, enc_dist)
  names(df) <- "GO_5085"
}
write.csv(df, "DH_GO_5085.csv")
DH_GO_5085_final <- cbind(fut_5085[, c("long", "lat")], df)
write.csv(DH_GO_5085_final, "DH_GO_5085_final.csv")
# Combine geographic coordinates and GO values and save

# Convert the genetic offset data into a raster for visualization
test_DH_GO_5085 <- read.csv("test_DH_GO_5085_final.csv")
test <- test_DH_GO_5085[, -1]  
DH5085_raster <- test %>% rasterFromXYZ(crs = "+proj=longlat +datum=WGS84 +no_defs")
pdf(file = "DHGO5085.pdf")
plot(DH5085_raster)
dev.off()
writeRaster(DH5085_raster, "DH5085_raster.asc", format = "ascii")
# Create and save a raster plot for the SSP585-50 scenario

### ssp585_90 Scenario
fut_9085 <- read.csv("test_115000_ssp585_90_mean+15pop.csv") 
Trns_grid_9085 <- cbind(fut_9085[, c("long", "lat")], predict(gfsnp_final, fut_9085[, imp.vars]))
write.csv(Trns_grid_9085, "Trns_grid_9085.csv")
# Predict genetic variation for the future scenario and combine with coordinates

## Calculate Genetic Offset (GO) values for the SSP585-90 scenario
df <- data.frame()
for (i in 1:nrow(fut_9085)) {
  enc_dist <- euclidean(Trns_grid[i, imp.vars], Trns_grid_9085[i, imp.vars])
  df <- rbind(df, enc_dist)
  names(df) <- "GO_9085"
}
write.csv(df, "DH_GO_9085.csv")
DH_GO_9085_final <- cbind(fut_9085[, c("long", "lat")], df)
write.csv(DH_GO_9085_final, "DH_GO_9085_final.csv")
# Combine geographic coordinates and GO values and save

# Convert the genetic offset data into a raster for visualization
test_DH_GO_9085 <- read.csv("test_DH_GO_9085_final.csv")
test <- test_DH_GO_9085[, -1]  
DH9085_raster <- test %>% rasterFromXYZ(crs = "+proj=longlat +datum=WGS84 +no_defs")
pdf(file = "DHGO9085.pdf")
plot(DH9085_raster)
dev.off()
writeRaster(DH9085_raster, "DH9085_raster.asc", format = "ascii")

# Output variable importance results from the Gradient Forest model
pdf(file = "all_predictoroverallimportance.pdf")
plot(gfsnp_final, plot.type = "Overall.Importance")
dev.off()
# Plot and save the overall importance of environmental predictors used in the model
write.csv(gfsnp_final$result, file = "all_gf_result.csv")
write.csv(gfsnp_final$overall.imp, file = "all_gf_overall.imp.csv") 
write.csv(gfsnp_final$overall.imp2, file = "all_gf_overall.imp2.csv")  







########################################################
 #### pan west
library(gradientForest)
# Read the dataset containing allele frequencies and environmental variables
gfData <- read.csv("pan_west_AF.csv")
envGF <- gfData[, 4:9]
# Extract the allele frequency data for the SNPs (columns containing "RpChr" in the header)
snp <- gfData[, grep("RpChr", colnames(gfData))]
# Check the number of SNPs
length(snp)  # Output: [1] 16709
# Calculate the maximum level parameter for the Gradient Forest model, accounting for correlations
maxLevel <- log2(0.368 * nrow(envGF) / 2)
# Build the Gradient Forest (GF) model
gfsnp_final <- gradientForest(cbind(envGF, snp), predictor.vars = colnames(envGF), response.vars = colnames(snp), 
                              ntree = 500, maxLevel = maxLevel, trace = TRUE, corr.threshold = 0.50)
# Read current environmental data for the locations (e.g., climate data)
DH_current_bio <- read.csv("test_115000_current+28pop.csv")
# Select the most important 6 BIO variables for the model
most_important <- c("BIO02", "BIO03", "BIO04", "BIO10", "BIO12", "BIO15")
imp.vars <- most_important
# Predict allele frequencies based on the selected environmental variables (for current conditions)
Trns_grid <- cbind(DH_current_bio[, c("long", "lat")], predict(gfsnp_final, DH_current_bio[, imp.vars]))
write.csv(Trns_grid, "Trns_grid.csv")

# Perform Principal Component Analysis (PCA) on the selected environmental variables (BIO variables)
PCs <- prcomp(Trns_grid[, imp.vars])

# Extract the first three principal components (PCs)
a1 <- PCs$x[, 1]
a2 <- PCs$x[, 2]
a3 <- PCs$x[, 3]

# Define RGB colors based on PCA components
r <- a1 + a2
g <- -a2
b <- a3 + a2 - a1

# Normalize the RGB values to range from 0 to 255
r <- (r - min(r)) / (max(r) - min(r)) * 255
g <- (g - min(g)) / (max(g) - min(g)) * 255
b <- (b - min(b)) / (max(b) - min(b)) * 255

# Plot the geographic distribution of genetic variation using the calculated RGB colors
pdf(file = "current distribution of genetic variation.pdf")
plot(Trns_grid[, c("long", "lat")], pch = ".", cex = 3, asp = 1, col = rgb(r, g, b, max = 255)) 
dev.off()

# Plot the PCA biplot showing the relationship between environmental variables
nvs <- dim(PCs$rotation)[1]  # Number of variables in PCA rotation
vec <- c("BIO02", "BIO03", "BIO04", "BIO10", "BIO12", "BIO15")  # The selected variables
lv <- length(vec)  # Length of the vector (number of variables)

# Get the indices for the variables to be plotted in the PCA
vind <- rownames(PCs$rotation) %in% vec

# Scale factor for better visualization of the PCA arrows
scal <- 15

# Define the range for the x and y axes of the PCA plot
xrng <- range(PCs$x[, 1], PCs$rotation[, 1] / scal) * 1.1
yrng <- range(PCs$x[, 2], PCs$rotation[, 2] / scal) * 1.1

# Plot the PCA biplot
pdf(file = "PCA biplot.pdf")
plot((PCs$x[, 1:2]), xlim = xrng, ylim = yrng, pch = ".", cex = 4, col = rgb(r, g, b, max = 255), asp = 1)
arrows(rep(0, lv), rep(0, lv), PCs$rotation[vec, 1] / scal, PCs$rotation[vec, 2] / scal, length = 0.0625, lwd = 1)
jit <- 0.004
text(PCs$rotation[vec, 1] / scal + jit * sign(PCs$rotation[vec, 1]), 
     PCs$rotation[vec, 2] / scal + jit * sign(PCs$rotation[vec, 2]), labels = vec, cex = 1, family = "serif")
dev.off()

################# Calculate genetic components under future important variables #################

### ssp126_50 
fut_5026 <- read.csv("test_115000_ssp126_50_mean+28pop.csv")
Trns_grid_5026 <- cbind(fut_5026[, c("long", "lat")], predict(gfsnp_final, fut_5026[, imp.vars]))
write.csv(Trns_grid_5026, "Trns_grid_5026.csv")

## Calculate Genetic Offset (GO) values using Euclidean distance
euclidean <- function(a, b) sqrt(sum((a - b)^2))  # Euclidean distance function for genetic offset
df <- data.frame()
for (i in 1:nrow(fut_5026)) {
    enc_dist <- euclidean(Trns_grid[i, imp.vars], Trns_grid_5026[i, imp.vars])  # Calculate the genetic offset for each point
    df <- rbind(df, enc_dist)
    names(df) <- "GO_5026"
}
write.csv(df, "DH_GO_5026.csv")

# Combine the geographic coordinates with the Genetic Offset (GO) values
DH_GO_5026_final <- cbind(fut_5026[, c("long", "lat")], df)
write.csv(DH_GO_5026_final, "DH_GO_5026_final.csv")

## After output, remove the first 43 data points, check the latitude column to ensure at least two samples per latitude, 
## and then save the cleaned file as "test_DH_GO_5026_final.csv"

# Convert the points of genetic offset into a raster format for visualization
library(dplyr)
library(raster)
test_DH_GO_5026 <- read.csv("test_DH_GO_5026_final.csv")

# Remove the first column (cell number)
test <- test_DH_GO_5026[, -1]

# Convert the data to raster format using the geographic coordinates and genetic offset
DH5026_raster <- test %>% rasterFromXYZ(crs = "+proj=longlat +datum=WGS84 +no_defs")

# Plot and save the raster showing genetic offset under the SSP126-50 scenario
pdf(file = "DHGO5026.pdf")
plot(DH5026_raster)
dev.off()
writeRaster(DH5026_raster, "DH5026_raster.asc", format = "ascii")

### SSP126_90
fut_9026 <- read.csv("test_115000_ssp126_90_mean+28pop.csv")
Trns_grid_9026 <- cbind(fut_9026[,c("long","lat")], predict(gfsnp_final,fut_9026[,imp.vars]))
write.csv(Trns_grid_9026,"Trns_grid_9026.csv")

# Calculate Genetic Offset (GO) values using Euclidean distance
euclidean <- function(a, b) sqrt(sum((a - b)^2))  # Euclidean distance function for GO
df <- data.frame()

# Loop through each row of the future climate data and calculate the GO values
for (i in 1:nrow(fut_9026)){
  enc_dist <- euclidean(Trns_grid[i,imp.vars],Trns_grid_9026[i,imp.vars])  # Calculate the distance between current and future allele frequencies
  df <- rbind(df, enc_dist)  # Store the calculated GO values
  names(df) <- "GO_9026"
}
write.csv(df,"DH_GO_9026.csv")
DH_GO_9026_final <- cbind(fut_9026[,c("long","lat")], df)
write.csv(DH_GO_9026_final,"DH_GO_9026_final.csv") 

# Convert the calculated GO values to a raster format for spatial visualization
library(dplyr)  # Ensure dplyr is loaded
test_DH_GO_9026 <- read.csv("test_DH_GO_9026_final.csv")  # Read the final GO data
test <- test_DH_GO_9026[,-1]  

# Convert the data to a raster format with appropriate CRS (Coordinate Reference System)
DH9026_raster <- test %>% rasterFromXYZ(crs="+proj=longlat +datum=WGS84 +no_defs")

# Plot the raster and save the output as a PDF
pdf(file="DHGO9026.pdf")
plot(DH9026_raster)
dev.off()
writeRaster(DH9026_raster,"DH9026_raster.asc", format="ascii")

##### SSP585_50 
fut_5085 <- read.csv("test_115000_ssp585_50_mean+28pop.csv")
Trns_grid_5085 <- cbind(fut_5085[,c("long","lat")], predict(gfsnp_final,fut_5085[,imp.vars]))
write.csv(Trns_grid_5085,"Trns_grid_5085.csv")

# Calculate the Genetic Offset (GO) values for SSP585 (50th percentile)
euclidean <- function(a, b) sqrt(sum((a - b)^2))  # Euclidean distance function for GO
df <- data.frame()

# Loop through each row of the future climate data and calculate the GO values
for (i in 1:nrow(fut_5085)){
  enc_dist <- euclidean(Trns_grid[i,imp.vars],Trns_grid_5085[i,imp.vars])  # Calculate GO
  df <- rbind(df, enc_dist)  # Store the results
  names(df) <- "GO_5085"
}
write.csv(df,"DH_GO_5085.csv")
DH_GO_5085_final <- cbind(fut_5085[,c("long","lat")], df)
write.csv(DH_GO_5085_final,"DH_GO_5085_final.csv")  

# Convert the GO results to raster format for visualization
test_DH_GO_5085 <- read.csv("test_DH_GO_5085_final.csv")  
test <- test_DH_GO_5085[,-1] 

# Convert to raster format
DH5085_raster <- test %>% rasterFromXYZ(crs="+proj=longlat +datum=WGS84 +no_defs")

# Plot the raster and save it as a PDF
pdf(file="DHGO5085.pdf")
plot(DH5085_raster)
dev.off()
writeRaster(DH5085_raster,"DH5085_raster.asc", format="ascii")

### SSP585_90 
fut_9085 <- read.csv("test_115000_ssp585_90_mean+28pop.csv")
Trns_grid_9085 <- cbind(fut_9085[,c("long","lat")], predict(gfsnp_final,fut_9085[,imp.vars]))
write.csv(Trns_grid_9085,"Trns_grid_9085.csv")

# Calculate the Genetic Offset (GO) values for SSP585 (90th percentile)
euclidean <- function(a, b) sqrt(sum((a - b)^2))  # Euclidean distance for GO
df <- data.frame()

# Loop to calculate GO values for each point in the future data
for (i in 1:nrow(fut_9085)){
  enc_dist <- euclidean(Trns_grid[i,imp.vars],Trns_grid_9085[i,imp.vars])  # Calculate GO
  df <- rbind(df, enc_dist)  # Store the results
  names(df) <- "GO_9085"
}
write.csv(df,"DH_GO_9085.csv")
DH_GO_9085_final <- cbind(fut_9085[,c("long","lat")], df)
write.csv(DH_GO_9085_final,"DH_GO_9085_final.csv") 

# Convert the GO results for SSP585 (90th percentile) to raster format
test_DH_GO_9085 <- read.csv("test_DH_GO_9085_final.csv")  
test <- test_DH_GO_9085[,-1] 

# Create a raster object for visualization
DH9085_raster <- test %>% rasterFromXYZ(crs="+proj=longlat +datum=WGS84 +no_defs")

# Plot the raster and save the output as a PDF
pdf(file="DHGO9085.pdf")
plot(DH9085_raster)
dev.off()
writeRaster(DH9085_raster,"DH9085_raster.asc", format="ascii")

# Output the overall importance of environmental predictors used in the model
pdf(file="all_predictoroverallimportance.pdf")
plot(gfsnp_final, plot.type="Overall.Importance")
dev.off()

# Save the results of the model and the importance of environmental variables
write.csv(gfsnp_final$result, file="all_gf_result.csv")  
write.csv(gfsnp_final$overall.imp, file="all_gf_overall.imp.csv")   
write.csv(gfsnp_final$overall.imp2, file="all_gf_overall.imp2.csv")  

 
 
 
 ########################################################
 ####core_loci
 library(gradientForest)
gfData <- read.csv("core_loci2_AF.csv")
envGF <- gfData[,4:9] # Read climate variable data from the loaded file
snp <- gfData[,grep("RpChr",colnames(gfData))] # Extract allele frequency data from columns that contain "RpChr"
length(snp) # Check the number of loci extracted
#[1] 1198
maxLevel <- log2(0.368*nrow(envGF)/2) # Account for correlations, see ?gradientForest
# Build GF model
gfsnp_final <- gradientForest(cbind(envGF, snp), predictor.vars=colnames(envGF), response.vars=colnames(snp), ntree=500, maxLevel=maxLevel, trace=T, corr.threshold=0.50)
DH_current_bio <- read.csv("test_115000_current+43pop.csv")
most_important <- c("BIO02","BIO03","BIO04","BIO10","BIO12","BIO15") # Output the 6 most important BIO variables
imp.vars <- most_important
Trns_grid <- cbind(DH_current_bio[,c("long","lat")], predict(gfsnp_final,DH_current_bio[,imp.vars])) # Calculate allele frequencies under the current 6 important environmental variables and combine with sample coordinates (longitude and latitude)
write.csv(Trns_grid,"Trns_grid.csv")
PCs <- prcomp(Trns_grid[,imp.vars])
a1 <- PCs$x[,1]
a2 <- PCs$x[,2]
a3 <- PCs$x[,3]
r <- a1 + a2
g <- -a2
b <- a3 + a2 - a1
r <- (r - min(r)) / (max(r) - min(r)) * 255
g <- (g - min(g)) / (max(g) - min(g)) * 255
b <- (b - min(b)) / (max(b) - min(b)) * 255
pdf(file="current distribution of genetic variation.pdf")
plot(Trns_grid[,c("long","lat")],pch=".", cex=3, asp=1, col=rgb(r,g,b, max = 255)) 
dev.off() # Output the current distribution of genetic variation map

# Output the PCA plot for current important environmental variables
nvs <- dim(PCs$rotation)[1]
vec <- c("BIO02","BIO03","BIO04","BIO10","BIO12","BIO15")
lv <- length(vec)
vind <- rownames(PCs$rotation) %in% vec
scal <- 15
xrng <- range(PCs$x[, 1], PCs$rotation[, 1]/scal) *1.1
yrng <- range(PCs$x[, 2], PCs$rotation[, 2]/scal) *1.1
pdf(file="PCA biplot.pdf")
plot((PCs$x[, 1:2]), xlim = xrng, ylim = yrng, pch = ".", cex = 4, col = rgb(r, g, b, max = 255), asp = 1)
arrows(rep(0, lv), rep(0, lv), PCs$rotation[vec,1]/scal, PCs$rotation[vec, 2]/scal, length = 0.0625, lwd = 1)
jit <- 0.004
text(PCs$rotation[vec, 1]/scal + jit * sign(PCs$rotation[vec,1]), PCs$rotation[vec, 2]/scal + jit * sign(PCs$rotation[vec,2]), labels = vec, cex = 1, family = "serif")
dev.off() # Output PCA plot for the current environment variables
################# Calculate genetic components under future important variables #################

###ssp126_50
fut_5026 <- read.csv("test_115000_ssp126_50_mean+43pop.csv") 
Trns_grid_5026 <- cbind(fut_5026[,c("long","lat")], predict(gfsnp_final,fut_5026[,imp.vars]))
write.csv(Trns_grid_5026,"Trns_grid_5026.csv")
# Calculate GO values (Genetic Offset)
euclidean <- function(a, b) sqrt(sum((a - b)^2))  # Euclidean distance for genetic offset
df <- data.frame()
for (i in 1:nrow(fut_5026)){
    enc_dist <- euclidean(Trns_grid[i,imp.vars],Trns_grid_5026[i,imp.vars])
    df <- rbind(df, enc_dist)
    names(df) <- "GO_5026"
}
write.csv(df,"DH_GO_5026.csv")
DH_GO_5026_final <- cbind(fut_5026[,c("long","lat")], df)
write.csv(DH_GO_5026_final,"DH_GO_5026_final.csv") # After outputting the results, remove the first 43 sample points, check the latitude column, ensure that each latitude repeats at least 2 times, and then save it as "test_DH_GO_5026_final.csv"
# Convert points of average value to raster
library(dplyr)  # detach("package:dplyr") to remove the loaded dplyr package from the working environment
library(raster)
test_DH_GO_5026 <- read.csv("test_DH_GO_5026_final.csv")
test <- test_DH_GO_5026[,-1] # Remove the first column (cell number)
DH5026_raster <- test %>% rasterFromXYZ(crs="+proj=longlat +datum=WGS84 +no_defs")
pdf(file="DHGO5026.pdf")
plot(DH5026_raster)
dev.off()
writeRaster(DH5026_raster,"DH5026_raster.asc", format="ascii")

###ssp126_90
fut_9026 <- read.csv("test_115000_ssp126_90_mean+43pop.csv") 
Trns_grid_9026 <- cbind(fut_9026[,c("long","lat")], predict(gfsnp_final,fut_9026[,imp.vars]))
write.csv(Trns_grid_9026,"Trns_grid_9026.csv")

## Calculate GO values
euclidean <- function(a, b) sqrt(sum((a - b)^2))  # Euclidean distance for genetic offset
df <- data.frame()
for (i in 1:nrow(fut_9026)){
  enc_dist <- euclidean(Trns_grid[i,imp.vars],Trns_grid_9026[i,imp.vars])
  df <- rbind(df, enc_dist)
  names(df) <- "GO_9026"
}
write.csv(df,"DH_GO_9026.csv")
DH_GO_9026_final <- cbind(fut_9026[,c("long","lat")], df)
write.csv(DH_GO_9026_final,"DH_GO_9026_final.csv") ## After outputting the results, remove the first 43 sample points, check the latitude column, ensure that the same latitude is repeated at least twice, then save this file as test_DH_GO_9026_final.csv

# Convert points of average value to raster
library(dplyr)  ## detach("package:dplyr"), to remove the loaded dplyr package from the working environment
test_DH_GO_9026 <- read.csv("test_DH_GO_9026_final.csv")
test <- test_DH_GO_9026[,-1] ## Remove the first column (cell number)
DH9026_raster <- test %>% rasterFromXYZ(crs="+proj=longlat +datum=WGS84 +no_defs")
pdf(file="DHGO9026.pdf")
plot(DH9026_raster)
dev.off()
writeRaster(DH9026_raster,"DH9026_raster.asc", format="ascii")

##### ssp585_50
fut_5085 <- read.csv("test_115000_ssp585_50_mean+43pop.csv") 
Trns_grid_5085 <- cbind(fut_5085[,c("long","lat")], predict(gfsnp_final,fut_5085[,imp.vars]))
write.csv(Trns_grid_5085,"Trns_grid_5085.csv")

## Calculate GO values
euclidean <- function(a, b) sqrt(sum((a - b)^2))  # Euclidean distance for genetic offset
df <- data.frame()
for (i in 1:nrow(fut_5085)){
  enc_dist <- euclidean(Trns_grid[i,imp.vars],Trns_grid_5085[i,imp.vars])
  df <- rbind(df, enc_dist)
  names(df) <- "GO_5085"
}
write.csv(df,"DH_GO_5085.csv")
DH_GO_5085_final <- cbind(fut_5085[,c("long","lat")], df)
write.csv(DH_GO_5085_final,"DH_GO_5085_final.csv") ## After outputting the results, remove the first 43 sample points, check the latitude column, ensure that the same latitude is repeated at least twice, then save this file as test_DH_GO_5085_final.csv

# Convert points of average value to raster
library(dplyr)  ## detach("package:dplyr"), to remove the loaded dplyr package from the working environment
test_DH_GO_5085 <- read.csv("test_DH_GO_5085_final.csv")
test <- test_DH_GO_5085[,-1] ## Remove the first column (cell number)
DH5085_raster <- test %>% rasterFromXYZ(crs="+proj=longlat +datum=WGS84 +no_defs")
pdf(file="DHGO5085.pdf")
plot(DH5085_raster)
dev.off()
writeRaster(DH5085_raster,"DH5085_raster.asc", format="ascii")

### ssp585_90
fut_9085 <- read.csv("test_115000_ssp585_90_mean+43pop.csv") 
Trns_grid_9085 <- cbind(fut_9085[,c("long","lat")], predict(gfsnp_final,fut_9085[,imp.vars]))
write.csv(Trns_grid_9085,"Trns_grid_9085.csv")

## Calculate GO values
euclidean <- function(a, b) sqrt(sum((a - b)^2))  # Euclidean distance for genetic offset
df <- data.frame()
for (i in 1:nrow(fut_9085)){
  enc_dist <- euclidean(Trns_grid[i,imp.vars],Trns_grid_9085[i,imp.vars])
  df <- rbind(df, enc_dist)
  names(df) <- "GO_9085"
}
write.csv(df,"DH_GO_9085.csv")
DH_GO_9085_final <- cbind(fut_9085[,c("long","lat")], df)
write.csv(DH_GO_9085_final,"DH_GO_9085_final.csv") ## After outputting the results, remove the first 43 sample points, check the latitude column, ensure that the same latitude is repeated at least twice, then save this file as test_DH_GO_9085_final.csv

# Convert points of average value to raster
library(dplyr)  ## detach("package:dplyr"), to remove the loaded dplyr package from the working environment
test_DH_GO_9085 <- read.csv("test_DH_GO_9085_final.csv")
test <- test_DH_GO_9085[,-1] ## Remove the first column (cell number)
DH9085_raster <- test %>% rasterFromXYZ(crs="+proj=longlat +datum=WGS84 +no_defs")
pdf(file="DHGO9085.pdf")
plot(DH9085_raster)
dev.off()
writeRaster(DH9085_raster,"DH9085_raster.asc", format="ascii")

# Output variable importance results
pdf(file="all_predictoroverallimportance.pdf")
plot(gfsnp_final, plot.type="Overall.Importance")
dev.off()

# Output variable importance tables
write.csv(gfsnp_final$result, file="all_gf_result.csv") 
write.csv(gfsnp_final$overall.imp, file="all_gf_overall.imp.csv")   # Output importance of all climate data used to build the model
write.csv(gfsnp_final$overall.imp2, file="all_gf_overall.imp2.csv")  # Output R2 weighted importance of all climate data used to build the model







###################################################################################################################
########## core_east
library(gradientForest)
gfData <- read.csv("core_east_AF.csv")
envGF <- gfData[, 4:9]  # Read the climate variable data from the loaded file
snp <- gfData[, grep("RpChr", colnames(gfData))]  # Read allele frequency data from columns containing "RpChr" in the header
length(snp)  # Check the number of extracted loci
#[1] 1198
maxLevel <- log2(0.368 * nrow(envGF) / 2)  # Account for correlations, see ?gradientForest
# Build the GF model
gfsnp_final <- gradientForest(cbind(envGF, snp), predictor.vars = colnames(envGF), response.vars = colnames(snp), ntree = 500, maxLevel = maxLevel, trace = TRUE, corr.threshold = 0.50)
DH_current_bio <- read.csv("test_115000_current+15pop.csv")
most_important <- c("BIO02", "BIO03", "BIO04", "BIO10", "BIO12", "BIO15")  # Output the six most important BIO variables
imp.vars <- most_important
Trns_grid <- cbind(DH_current_bio[, c("long", "lat")], predict(gfsnp_final, DH_current_bio[, imp.vars]))  # Calculate allele frequencies under the current 6 most important environmental variables and combine results with sample coordinates
write.csv(Trns_grid, "Trns_grid.csv")
PCs <- prcomp(Trns_grid[, imp.vars])
a1 <- PCs$x[, 1]
a2 <- PCs$x[, 2]
a3 <- PCs$x[, 3]
r <- a1 + a2
g <- -a2
b <- a3 + a2 - a1
r <- (r - min(r)) / (max(r) - min(r)) * 255
g <- (g - min(g)) / (max(g) - min(g)) * 255
b <- (b - min(b)) / (max(b) - min(b)) * 255
pdf(file = "current distribution of genetic variation.pdf")
plot(Trns_grid[, c("long", "lat")], pch = ".", cex = 3, asp = 1, col = rgb(r, g, b, max = 255)) 
dev.off()  # Output a map of the geographical distribution of current genetic variation

## Output the PCA plot for the current important environmental variables
nvs <- dim(PCs$rotation)[1]
vec <- c("BIO02", "BIO03", "BIO04", "BIO10", "BIO12", "BIO15")
lv <- length(vec)
vind <- rownames(PCs$rotation) %in% vec
scal <- 15
xrng <- range(PCs$x[, 1], PCs$rotation[, 1] / scal) * 1.1
yrng <- range(PCs$x[, 2], PCs$rotation[, 2] / scal) * 1.1
pdf(file = "PCA biplot.pdf")
plot((PCs$x[, 1:2]), xlim = xrng, ylim = yrng, pch = ".", cex = 4, col = rgb(r, g, b, max = 255), asp = 1)
arrows(rep(0, lv), rep(0, lv), PCs$rotation[vec, 1] / scal, PCs$rotation[vec, 2] / scal, length = 0.0625, lwd = 1)
jit <- 0.004
text(PCs$rotation[vec, 1] / scal + jit * sign(PCs$rotation[vec, 1]), PCs$rotation[vec, 2] / scal + jit * sign(PCs$rotation[vec, 2]), labels = vec, cex = 1, family = "serif")
dev.off()  # Output the PCA plot for current important variables

################# Calculate genetic components under future important variables #################

### ssp126_50
fut_5026 <- read.csv("test_115000_ssp126_50_mean+15pop.csv") 
Trns_grid_5026 <- cbind(fut_5026[, c("long", "lat")], predict(gfsnp_final, fut_5026[, imp.vars]))
write.csv(Trns_grid_5026, "Trns_grid_5026.csv")

## Calculate GO values
euclidean <- function(a, b) sqrt(sum((a - b)^2))  # Euclidean distance for genetic offset
df <- data.frame()
for (i in 1:nrow(fut_5026)) {
    enc_dist <- euclidean(Trns_grid[i, imp.vars], Trns_grid_5026[i, imp.vars])
    df <- rbind(df, enc_dist)
    names(df) <- "GO_5026"
}
write.csv(df, "DH_GO_5026.csv")
DH_GO_5026_final <- cbind(fut_5026[, c("long", "lat")], df)
write.csv(DH_GO_5026_final, "DH_GO_5026_final.csv")  # After outputting the results, remove the first 43 sample points, check the latitude column, ensure the same latitude is repeated at least twice, then save this file as test_DH_GO_5026_final.csv

# Convert points of average value to raster
library(dplyr)  ## detach("package:dplyr"), to remove the loaded dplyr package from the working environment
library(raster)
test_DH_GO_5026 <- read.csv("test_DH_GO_5026_final.csv")
test <- test_DH_GO_5026[, -1]  ## Remove the first column (cell number)
DH5026_raster <- test %>% rasterFromXYZ(crs = "+proj=longlat +datum=WGS84 +no_defs")
pdf(file = "DHGO5026.pdf")
plot(DH5026_raster)
dev.off()
writeRaster(DH5026_raster, "DH5026_raster.asc", format = "ascii")

### ssp126_90
fut_9026 <- read.csv("test_115000_ssp126_90_mean+15pop.csv") 
Trns_grid_9026 <- cbind(fut_9026[, c("long", "lat")], predict(gfsnp_final, fut_9026[, imp.vars]))
write.csv(Trns_grid_9026, "Trns_grid_9026.csv")

## Calculate GO values
euclidean <- function(a, b) sqrt(sum((a - b)^2))  # Euclidean distance for genetic offset
df <- data.frame()
for (i in 1:nrow(fut_9026)) {
  enc_dist <- euclidean(Trns_grid[i, imp.vars], Trns_grid_9026[i, imp.vars])
  df <- rbind(df, enc_dist)
  names(df) <- "GO_9026"
}
write.csv(df, "DH_GO_9026.csv")
DH_GO_9026_final <- cbind(fut_9026[, c("long", "lat")], df)
write.csv(DH_GO_9026_final, "DH_GO_9026_final.csv")  # After outputting the results, remove the first 43 sample points, check the latitude column, ensure the same latitude is repeated at least twice, then save this file as test_DH_GO_9026_final.csv

# Convert points of average value to raster
library(dplyr)  ## detach("package:dplyr"), to remove the loaded dplyr package from the working environment
test_DH_GO_9026 <- read.csv("test_DH_GO_9026_final.csv")
test <- test_DH_GO_9026[,-1] ## Remove the first column (cell number)
DH9026_raster <- test %>% rasterFromXYZ(crs = "+proj=longlat +datum=WGS84 +no_defs")
pdf(file = "DHGO9026.pdf")
plot(DH9026_raster)
dev.off()
writeRaster(DH9026_raster, "DH9026_raster.asc", format = "ascii")

##### 5085
fut_5085 <- read.csv("test_115000_ssp585_50_mean+15pop.csv") 
Trns_grid_5085 <- cbind(fut_5085[, c("long", "lat")], predict(gfsnp_final, fut_5085[, imp.vars]))
write.csv(Trns_grid_5085, "Trns_grid_5085.csv")

## Calculate GO values
euclidean <- function(a, b) sqrt(sum((a - b)^2))  # Euclidean distance for genetic offset
df <- data.frame()
for (i in 1:nrow(fut_5085)) {
  enc_dist <- euclidean(Trns_grid[i, imp.vars], Trns_grid_5085[i, imp.vars])
  df <- rbind(df, enc_dist)
  names(df) <- "GO_5085"
}
write.csv(df, "DH_GO_5085.csv")
DH_GO_5085_final <- cbind(fut_5085[, c("long", "lat")], df)
write.csv(DH_GO_5085_final, "DH_GO_5085_final.csv") ## After outputting the results, remove the first 43 sample points, check the latitude column, ensure the same latitude is repeated at least twice, then save the file as test_DH_GO_5085_final.csv

# Convert points of average value to raster
library(dplyr)  ## detach("package:dplyr"), to remove the loaded dplyr package from the working environment
test_DH_GO_5085 <- read.csv("test_DH_GO_5085_final.csv")
test <- test_DH_GO_5085[,-1] ## Remove the first column (cell number)
DH5085_raster <- test %>% rasterFromXYZ(crs = "+proj=longlat +datum=WGS84 +no_defs")
pdf(file = "DHGO5085.pdf")
plot(DH5085_raster)
dev.off()
writeRaster(DH5085_raster, "DH5085_raster.asc", format = "ascii")

### ssp585_90
fut_9085 <- read.csv("test_115000_ssp585_90_mean+15pop.csv") 
Trns_grid_9085 <- cbind(fut_9085[, c("long", "lat")], predict(gfsnp_final, fut_9085[, imp.vars]))
write.csv(Trns_grid_9085, "Trns_grid_9085.csv")

## Calculate GO values
euclidean <- function(a, b) sqrt(sum((a - b)^2))  # Euclidean distance for genetic offset
df <- data.frame()
for (i in 1:nrow(fut_9085)) {
  enc_dist <- euclidean(Trns_grid[i, imp.vars], Trns_grid_9085[i, imp.vars])
  df <- rbind(df, enc_dist)
  names(df) <- "GO_9085"
}
write.csv(df, "DH_GO_9085.csv")
DH_GO_9085_final <- cbind(fut_9085[, c("long", "lat")], df)
write.csv(DH_GO_9085_final, "DH_GO_9085_final.csv") ## After outputting the results, remove the first 43 sample points, check the latitude column, ensure the same latitude is repeated at least twice, then save the file as test_DH_GO_9085_final.csv

# Convert points of average value to raster
library(dplyr)  ## detach("package:dplyr"), to remove the loaded dplyr package from the working environment
test_DH_GO_9085 <- read.csv("test_DH_GO_9085_final.csv")
test <- test_DH_GO_9085[,-1] ## Remove the first column (cell number)
DH9085_raster <- test %>% rasterFromXYZ(crs = "+proj=longlat +datum=WGS84 +no_defs")
pdf(file = "DHGO9085.pdf")
plot(DH9085_raster)
dev.off()
writeRaster(DH9085_raster, "DH9085_raster.asc", format = "ascii")

# Output variable importance results
pdf(file = "all_predictoroverallimportance.pdf")
plot(gfsnp_final, plot.type = "Overall.Importance")
dev.off()

# Output variable importance tables
write.csv(gfsnp_final$result, file = "all_gf_result.csv") 
write.csv(gfsnp_final$overall.imp, file = "all_gf_overall.imp.csv")   # Output importance of all climate data used to construct the model
write.csv(gfsnp_final$overall.imp2, file = "all_gf_overall.imp2.csv")  # Output R2-weighted importance of all climate data used to construct the model






########################################################
 #### core west
library(gradientForest)
gfData <- read.csv("core_west_AF.csv")
envGF <- gfData[,4:9] # Read the climate variable data from the loaded file
snp <- gfData[,grep("RpChr",colnames(gfData))] # Extract the allele frequency data based on column names containing "RpChr"
length(snp) # Check the number of loci extracted
#[1] 1198
maxLevel <- log2(0.368*nrow(envGF)/2) # Account for correlations, see ?gradientForest
# Build GF model
gfsnp_final <- gradientForest(cbind(envGF, snp), predictor.vars = colnames(envGF), response.vars = colnames(snp), ntree = 500, maxLevel = maxLevel, trace = T, corr.threshold = 0.50)
DH_current_bio <- read.csv("test_115000_current+28pop.csv")
most_important <- c("BIO02", "BIO03", "BIO04", "BIO10", "BIO12", "BIO15") ## Output the 6 most important BIO variables
imp.vars <- most_important
Trns_grid <- cbind(DH_current_bio[, c("long", "lat")], predict(gfsnp_final, DH_current_bio[, imp.vars])) # Calculate allele frequencies under the current 6 important environmental variables, and combine the genomic composition results with the sample coordinates into one file
write.csv(Trns_grid, "Trns_grid.csv")
PCs <- prcomp(Trns_grid[, imp.vars])
a1 <- PCs$x[, 1]
a2 <- PCs$x[, 2]
a3 <- PCs$x[, 3]
r <- a1 + a2
g <- -a2
b <- a3 + a2 - a1
r <- (r - min(r)) / (max(r) - min(r)) * 255
g <- (g - min(g)) / (max(g) - min(g)) * 255
b <- (b - min(b)) / (max(b) - min(b)) * 255
pdf(file = "current distribution of genetic variation.pdf")
plot(Trns_grid[, c("long", "lat")], pch = ".", cex = 3, asp = 1, col = rgb(r, g, b, max = 255)) 
dev.off() ## Output the geographic distribution map of contemporary genetic variation

## Output the PCA plot under the current important environmental variables
nvs <- dim(PCs$rotation)[1]
vec <- c("BIO02", "BIO03", "BIO04", "BIO10", "BIO12", "BIO15")
lv <- length(vec)
vind <- rownames(PCs$rotation) %in% vec
scal <- 15
xrng <- range(PCs$x[, 1], PCs$rotation[, 1] / scal) * 1.1
yrng <- range(PCs$x[, 2], PCs$rotation[, 2] / scal) * 1.1
pdf(file = "PCA biplot.pdf")
plot((PCs$x[, 1:2]), xlim = xrng, ylim = yrng, pch = ".", cex = 4, col = rgb(r, g, b, max = 255), asp = 1)
arrows(rep(0, lv), rep(0, lv), PCs$rotation[vec, 1] / scal, PCs$rotation[vec, 2] / scal, length = 0.0625, lwd = 1)
jit <- 0.004
text(PCs$rotation[vec, 1] / scal + jit * sign(PCs$rotation[vec, 1]), PCs$rotation[vec, 2] / scal + jit * sign(PCs$rotation[vec, 2]), labels = vec, cex = 1, family = "serif")
dev.off() ## Output the PCA plot for the current important environmental variables

################# Calculate genetic components under future important variables #################
### ssp126_50
fut_5026 <- read.csv("test_115000_ssp126_50_mean+28pop.csv") 
Trns_grid_5026 <- cbind(fut_5026[, c("long", "lat")], predict(gfsnp_final, fut_5026[, imp.vars]))
write.csv(Trns_grid_5026, "Trns_grid_5026.csv")
## Calculate GO values
euclidean <- function(a, b) sqrt(sum((a - b)^2))  # Euclidean distance for genetic offset
df <- data.frame()
for (i in 1:nrow(fut_5026)) {
  enc_dist <- euclidean(Trns_grid[i, imp.vars], Trns_grid_5026[i, imp.vars])
  df <- rbind(df, enc_dist)
  names(df) <- "GO_5026"
}
write.csv(df, "DH_GO_5026.csv")
DH_GO_5026_final <- cbind(fut_5026[, c("long", "lat")], df)
write.csv(DH_GO_5026_final, "DH_GO_5026_final.csv") ## After outputting the results, remove the first 43 sample points, check the latitude column, ensure the same latitude is repeated at least twice, then save the file as test_DH_GO_5026_final.csv

# Convert points of average value to raster
library(dplyr)  ## detach("package:dplyr"), to remove the loaded dplyr package from the working environment
library(raster)
test_DH_GO_5026 <- read.csv("test_DH_GO_5026_final.csv")
test <- test_DH_GO_5026[, -1] ## Remove the first column (cell number)
DH5026_raster <- test %>% rasterFromXYZ(crs = "+proj=longlat +datum=WGS84 +no_defs")
pdf(file = "DHGO5026.pdf")
plot(DH5026_raster)
dev.off()
writeRaster(DH5026_raster, "DH5026_raster.asc", format = "ascii")

### ssp126_90
fut_9026 <- read.csv("test_115000_ssp126_90_mean+28pop.csv") 
Trns_grid_9026 <- cbind(fut_9026[, c("long", "lat")], predict(gfsnp_final, fut_9026[, imp.vars]))
write.csv(Trns_grid_9026, "Trns_grid_9026.csv")

## Calculate GO values
euclidean <- function(a, b) sqrt(sum((a - b)^2))  # Euclidean distance for genetic offset
df <- data.frame()
for (i in 1:nrow(fut_9026)) {
  enc_dist <- euclidean(Trns_grid[i, imp.vars], Trns_grid_9026[i, imp.vars])
  df <- rbind(df, enc_dist)
  names(df) <- "GO_9026"
}
write.csv(df, "DH_GO_9026.csv")
DH_GO_9026_final <- cbind(fut_9026[, c("long", "lat")], df)
write.csv(DH_GO_9026_final, "DH_GO_9026_final.csv") ## After outputting the results, remove the first 43 sample points, check the latitude column, ensure the same latitude is repeated at least twice, then save the file as test_DH_GO_9026_final.csv

# Convert points of average value to raster
library(dplyr)  ## detach("package:dplyr"), to remove the loaded dplyr package from the working environment
test_DH_GO_9026 <- read.csv("test_DH_GO_9026_final.csv")
test <- test_DH_GO_9026[, -1] ## Remove the first column (cell number)
DH9026_raster <- test %>% rasterFromXYZ(crs = "+proj=longlat +datum=WGS84 +no_defs")
pdf(file = "DHGO9026.pdf")
plot(DH9026_raster)
dev.off()
writeRaster(DH9026_raster, "DH9026_raster.asc", format = "ascii")


##### 5085
fut_5085 <- read.csv("test_115000_ssp585_50_mean+28pop.csv") 
Trns_grid_5085 <- cbind(fut_5085[,c("long","lat")], predict(gfsnp_final, fut_5085[,imp.vars]))
write.csv(Trns_grid_5085,"Trns_grid_5085.csv")

## Calculate GO values
euclidean <- function(a, b) sqrt(sum((a - b)^2))  # Euclidean distance for genetic offset
df <- data.frame()
for (i in 1:nrow(fut_5085)) {
  enc_dist <- euclidean(Trns_grid[i, imp.vars], Trns_grid_5085[i, imp.vars])
  df <- rbind(df, enc_dist)
  names(df) <- "GO_5085"
}
write.csv(df, "DH_GO_5085.csv")
DH_GO_5085_final <- cbind(fut_5085[, c("long", "lat")], df)
write.csv(DH_GO_5085_final, "DH_GO_5085_final.csv")  ## After outputting the results, remove the first 43 sample points, check the latitude column, ensure the same latitude is repeated at least twice, then save the file as test_DH_GO_5085_final.csv

# Convert points of average value to raster
library(dplyr)  ## detach("package:dplyr"), remove the loaded dplyr package from the working environment
test_DH_GO_5085 <- read.csv("test_DH_GO_5085_final.csv")
test <- test_DH_GO_5085[,-1]  ## Remove the first column (cell number)
DH5085_raster <- test %>% rasterFromXYZ(crs = "+proj=longlat +datum=WGS84 +no_defs")
pdf(file = "DHGO5085.pdf")
plot(DH5085_raster)
dev.off()
writeRaster(DH5085_raster, "DH5085_raster.asc", format = "ascii")

### ssp585_90
fut_9085 <- read.csv("test_115000_ssp585_90_mean+28pop.csv") 
Trns_grid_9085 <- cbind(fut_9085[, c("long", "lat")], predict(gfsnp_final, fut_9085[, imp.vars]))
write.csv(Trns_grid_9085, "Trns_grid_9085.csv")

## Calculate GO values
euclidean <- function(a, b) sqrt(sum((a - b)^2))  # Euclidean distance for genetic offset
df <- data.frame()
for (i in 1:nrow(fut_9085)) {
  enc_dist <- euclidean(Trns_grid[i, imp.vars], Trns_grid_9085[i, imp.vars])
  df <- rbind(df, enc_dist)
  names(df) <- "GO_9085"
}
write.csv(df, "DH_GO_9085.csv")
DH_GO_9085_final <- cbind(fut_9085[, c("long", "lat")], df)
write.csv(DH_GO_9085_final, "DH_GO_9085_final.csv")  ## After outputting the results, remove the first 43 sample points, check the latitude column, ensure the same latitude is repeated at least twice, then save the file as test_DH_GO_9085_final.csv

# Convert points of average value to raster
library(dplyr)  ## detach("package:dplyr"), remove the loaded dplyr package from the working environment
test_DH_GO_9085 <- read.csv("test_DH_GO_9085_final.csv")
test <- test_DH_GO_9085[,-1]  ## Remove the first column (cell number)
DH9085_raster <- test %>% rasterFromXYZ(crs = "+proj=longlat +datum=WGS84 +no_defs")
pdf(file = "DHGO9085.pdf")
plot(DH9085_raster)
dev.off()
writeRaster(DH9085_raster, "DH9085_raster.asc", format = "ascii")

# Output variable importance results
pdf(file = "all_predictoroverallimportance.pdf")
plot(gfsnp_final, plot.type = "Overall.Importance")
dev.off()

# Output the importance table of variables
write.csv(gfsnp_final$result, file = "all_gf_result.csv") 
write.csv(gfsnp_final$overall.imp, file = "all_gf_overall.imp.csv")   # Output the importance of all the climate data used to construct the model
write.csv(gfsnp_final$overall.imp2, file = "all_gf_overall.imp2.csv")  # Output the R2-weighted importance of the climate data used to construct the model










######################################作图
## All populations_pan
library(raster)
DH5026_raster <- raster("DH5026_raster.asc")
DH9026_raster <- raster("DH9026_raster.asc")
DH5085_raster <- raster("DH5085_raster.asc")
DH9085_raster <- raster("DH9085_raster.asc")

pdf(file = "DH_all_GO.pdf")
par(mfrow = c(2, 2), mar = c(2, 2, 2, 1), oma = c(1, 1, 1, 2))
# First subplot
plot(DH5026_raster, zlim = c(0.10, 0.26), legend = FALSE)  # Plot DH5026_raster, set z-axis range, hide legend
points_data <- read.csv("pan_allgo.csv")  # Read GO values from CSV file
library(sp)  # Load spatial data processing library
# Create spatial points data frame
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object
# Set coordinate reference system (assuming longitude and latitude use WGS84)
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
plot(DH9026_raster, zlim = c(0.10, 0.26), legend = FALSE)  # Plot DH9026_raster
points_data <- read.csv("pan_allgo.csv")  # Read GO values from CSV file
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
plot(DH5085_raster, zlim = c(0.10, 0.26), legend = FALSE)  # Plot DH5085_raster
points_data <- read.csv("pan_allgo.csv")  # Read GO values from CSV file
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
plot(DH9085_raster, zlim = c(0.10, 0.26))  # Plot DH9085_raster and display legend
points_data <- read.csv("pan_allgo.csv")  # Read GO values from CSV file
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



## All populations_core
library(raster)
DH5026_raster <- raster("DH5026_raster.asc")
DH9026_raster <- raster("DH9026_raster.asc")
DH5085_raster <- raster("DH5085_raster.asc")
DH9085_raster <- raster("DH9085_raster.asc")

pdf(file = "DH_all_GO.pdf")
par(mfrow = c(2, 2), mar = c(2, 2, 2, 1), oma = c(1, 1, 1, 2))
# First subplot
plot(DH5026_raster, zlim = c(0.10, 0.36), legend = FALSE)  # Plot DH5026_raster, set z-axis range, hide legend
points_data <- read.csv("CoreGO_ALL.csv")  # Read GO values from CSV file
library(sp)  # Load spatial data processing library
# Create spatial points data frame
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object
# Set coordinate reference system (assuming longitude and latitude use WGS84)
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
plot(DH9026_raster, zlim = c(0.10, 0.36), legend = FALSE)  # Plot DH9026_raster
points_data <- read.csv("CoreGO_ALL.csv")  # Read GO values from CSV file
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
plot(DH5085_raster, zlim = c(0.10, 0.36), legend = FALSE)  # Plot DH5085_raster
points_data <- read.csv("CoreGO_ALL.csv")  # Read GO values from CSV file
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
plot(DH9085_raster, zlim = c(0.10, 0.36))  # Plot DH9085_raster and show legend
points_data <- read.csv("CoreGO_ALL.csv")  # Read GO values from CSV file
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




## East pan
library(raster)
DH5026_raster <- raster("DH5026_paneast.tif")
DH9026_raster <- raster("DH9026_paneast.tif")
DH5085_raster <- raster("DH5085_paneast.tif")
DH9085_raster <- raster("DH9085_paneast.tif")

pdf(file = "DH_all_GO2.pdf")
par(mfrow = c(2, 2), mar = c(2, 2, 2, 1), oma = c(1, 1, 1, 2))
# First subplot
plot(DH5026_raster, zlim = c(0.10, 0.26), legend = FALSE)  # Plot DH5026_raster, set z-axis range, hide legend
points_data <- read.csv("paneast_go.csv")  # Read GO value data from the CSV file
library(sp)  # Load spatial data processing library
# Create spatial points data frame
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object
# Set coordinate reference system (assuming longitude and latitude use WGS84)
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
# Calculate maximum and minimum GO values
max_GO <- max(points_data$GO_5026)  # Get the maximum GO value
min_GO <- min(points_data$GO_5026)  # Get the minimum GO value
# Calculate symbol radius
points_data$radius <- 0.05 + 0.3 * (points_data$GO_5026 - min_GO) / (max_GO - min_GO)
# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(a)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label
# Second subplot
plot(DH9026_raster, zlim = c(0.10, 0.26), legend = FALSE)  # Plot DH9026_raster
points_data <- read.csv("paneast_go.csv")  # Read GO value data from the CSV file
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
plot(DH5085_raster, zlim = c(0.10, 0.26), legend = FALSE)  # Plot DH5085_raster
points_data <- read.csv("paneast_go.csv")  # Read GO value data from the CSV file
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
plot(DH9085_raster, zlim = c(0.10, 0.26))  # Plot DH9085_raster and show legend
points_data <- read.csv("paneast_go.csv")  # Read GO value data from the CSV file
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



## West pan
library(raster)
DH5026_raster <- raster("DH5026_panwest.tif")
DH9026_raster <- raster("DH9026_panwest.tif")
DH5085_raster <- raster("DH5085_panwest.tif")
DH9085_raster <- raster("DH9085_panwest.tif")

pdf(file = "DH_all_GO.pdf")
par(mfrow = c(2, 2), mar = c(2, 2, 2, 1), oma = c(1, 1, 1, 2))
# First subplot
plot(DH5026_raster, zlim = c(0.10, 0.26), legend = FALSE)  # Plot DH5026_raster, set z-axis range, hide legend
points_data <- read.csv("panwest_go.csv")  # Read GO values data from the CSV file
library(sp)  # Load spatial data processing library
# Create spatial points data frame
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object
# Set coordinate reference system (assuming longitude and latitude use WGS84)
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
# Calculate maximum and minimum GO values
max_GO <- max(points_data$GO_5026)  # Get the maximum GO value
min_GO <- min(points_data$GO_5026)  # Get the minimum GO value
# Calculate symbol radius
points_data$radius <- 0.05 + 0.3 * (points_data$GO_5026 - min_GO) / (max_GO - min_GO)
# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(a)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label
# Second subplot
plot(DH9026_raster, zlim = c(0.10, 0.26), legend = FALSE)  # Plot DH9026_raster
points_data <- read.csv("panwest_go.csv")  # Read GO values data from the CSV file
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
plot(DH5085_raster, zlim = c(0.10, 0.26), legend = FALSE)  # Plot DH5085_raster
points_data <- read.csv("panwest_go.csv")  # Read GO values data from the CSV file
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
plot(DH9085_raster, zlim = c(0.10, 0.26))  # Plot DH9085_raster and show legend
points_data <- read.csv("panwest_go.csv")  # Read GO values data from the CSV file
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



## Core
## West
library(raster)
DH5026_raster <- raster("DH5026_corewest.tif")
DH9026_raster <- raster("DH9026_corewest.tif")
DH5085_raster <- raster("DH5085_corewest.tif")
DH9085_raster <- raster("DH9085_corewest.tif")

pdf(file = "DH_all_GO.pdf")
par(mfrow = c(2, 2), mar = c(2, 2, 2, 1), oma = c(1, 1, 1, 2))
# First subplot
plot(DH5026_raster, zlim = c(0.10, 0.36), legend = FALSE)  # Plot DH5026_raster, set z-axis range, hide legend
points_data <- read.csv("corewest_go.csv")  # Read GO values data from CSV file
library(sp)  # Load spatial data processing library
# Create spatial points data frame
coordinates(points_data) <- c("long", "lat")  # Define longitude and latitude columns
points_sp <- points_data  # Convert to SpatialPoints object
# Set coordinate reference system (assuming longitude and latitude use WGS84)
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
# Calculate maximum and minimum GO values
max_GO <- max(points_data$GO_5026)  # Get the maximum GO value
min_GO <- min(points_data$GO_5026)  # Get the minimum GO value
# Calculate symbol radius
points_data$radius <- 0.05 + 0.3 * (points_data$GO_5026 - min_GO) / (max_GO - min_GO)
# Add points to the plot
symbols(points_data$long, points_data$lat, circles = points_data$radius, inches = FALSE, add = TRUE, bg = "blue")
mtext("(a)", side = 3, line = -2, at = par("usr")[1], adj = 0, cex = 1.2)  # Add subplot label
# Second subplot
plot(DH9026_raster, zlim = c(0.10, 0.36), legend = FALSE)  # Plot DH9026_raster
points_data <- read.csv("corewest_go.csv")  # Read GO values data from CSV file
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
plot(DH5085_raster, zlim = c(0.10, 0.36), legend = FALSE)  # Plot DH5085_raster
points_data <- read.csv("corewest_go.csv")  # Read GO values data from CSV file
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
plot(DH9085_raster, zlim = c(0.10, 0.36))  # Plot DH9085_raster and display the legend
points_data <- read.csv("corewest_go.csv")  # Read GO values data from CSV file
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



# Core
## East
library(raster)
DH5026_raster <- raster("DH5026_coreeast.tif")
DH9026_raster <- raster("DH9026_coreeast.tif")
DH5085_raster <- raster("DH5085_coreeast.tif")
DH9085_raster <- raster("DH9085_coreeast.tif")

pdf(file = "DH_all_GO.pdf")
par(mfrow = c(2, 2), mar = c(2, 2, 2, 1), oma = c(1, 1, 1, 2))
# First subplot
plot(DH5026_raster, zlim = c(0.10, 0.36), legend = FALSE)  # Plot DH5026_raster, set z-axis range, hide legend
points_data <- read.csv("coreeast_go.csv")  # Read GO values data from the CSV file
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
plot(DH9026_raster, zlim = c(0.10, 0.36), legend = FALSE)  # Plot DH9026_raster
points_data <- read.csv("coreeast_go.csv")  # Read GO values data from the CSV file
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
plot(DH5085_raster, zlim = c(0.10, 0.36), legend = FALSE)  # Plot DH5085_raster
points_data <- read.csv("coreeast_go.csv")  # Read GO values data from the CSV file
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
plot(DH9085_raster, zlim = c(0.10, 0.36))  # Plot DH9085_raster and display legend
points_data <- read.csv("coreeast_go.csv")  # Read GO values data from the CSV file
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
