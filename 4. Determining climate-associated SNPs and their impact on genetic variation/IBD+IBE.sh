## Preparing files for adaptive and neutral loci
# Convert the input PLINK binary file to VCF format, output as 213snp3.vcf
plink --bfile 213snp3 --export vcf --out 213snp3 --allow-extra-chr
plink --bfile 213snp3 --extract core_snp.txt --make-bed --out core_adapt
plink --bfile core_adapt --recode structure --out core_adapt --allow-extra-chr  ## Outputs a STRUCTURE (.str) file


plink --bfile 213snp3 --extract pan_snp.txt --make-bed --out pan_adapt --allow-extra-chr
plink --bfile pan_adapt --recode structure --out pan_adapt --allow-extra-chr 

# Extract random neutral loci and create a new binary file
plink --bfile 213snp3 --extract random_samples_200000.txt --make-bed --out nutural --allow-extra-chr
# Recode the neutral dataset to STRUCTURE format
plink --bfile nutural --recode structure --out nutural --allow-extra-chr  ## Outputs a STRUCTURE (.str) file

# Convert the STRUCTURE format files by replacing spaces with tabs
tr ' ' '\t' < core_adapt.recode.strct_in > core_adapt.str
tr ' ' '\t' < pan_adapt.recode.strct_in > pan_adapt.str
tr ' ' '\t' < nutural.recode.strct_in > nutural.str  ## Process .str files by replacing spaces with tabs




library(hierfstat)  # For working with genetic data in the hierfstat format
library(reshape2)   # For reshaping data frames
library(adegenet)   # For reading .str data files
###################Neutral loci IBD (Identity by Descent) + IBE (Identity by State)###############################
str <- read.structure("nutural.str", n.ind = 201, n.loc = 200000, col.lab = 1, col.pop = 2, ask = FALSE, quiet = FALSE)
# Convert data from a STRUCTURE .str file to a genind object...
data <- genind2df(str)  # Convert STRUCTURE data to a genind data frame
write.csv(data, "neutral_genind2df.csv")  # Save the genind formatted data frame to a CSV file
test <- genind2hierfstat(str)  # Convert the genind object to a hierfstat format for genetic statistics
fst <- genet.dist(test, method = "Ds")  # Compute genetic distance using Nei's standard genetic distance method
# Function to reshape the distance matrix to a table and remove redundant data
fun2 <- function(inDist) {
  df <- melt(as.matrix(inDist), varnames = c("row", "col"))  # Melt the distance matrix into a long format table
  df[as.numeric(df$row) > as.numeric(df$col), ]  # Remove duplicated distances (only keep one triangle of the matrix)
}
col_GD <- fun2(fst)  # Apply the function to the genetic distance matrix
write.csv(col_GD, "col_GD_neutral.csv") 
Fst <- fst / (1 - fst)  # Compute Fst/(1-Fst) to adjust genetic distance values
Fst_matrix <- as.data.frame(as.matrix(Fst))  
write.csv(Fst_matrix, "Fst_matrix_neutral.csv")  
 

library(vegan)  
library(geosphere)  
env = read.csv("WGR_env_6.csv", head = T) 
ENV = as.data.frame(env[, c(4:9)]) 
ENV <- scale(ENV, center = TRUE, scale = TRUE)  # Scale the environmental data
rownames(ENV) = env$POP  # Set row names to population names
env_dist = vegdist(ENV, method = "euclidean", binary = FALSE, diag = FALSE, upper = FALSE, na.rm = FALSE)  # Compute Euclidean distance between samples
write.csv(as.matrix(env_dist), file = "env_scale_panadapt.csv", quote = F) 
fst = read.csv("Fst_matrix_neutral.csv", head = T, row.names = 1)  
geo_dist = read.csv("geo_dist_panadapt.csv", head = T, row.names = 1)  
# Perform Mantel test to correlate geographical distance with environmental distance
env_dist = as.matrix(env_dist)  # Convert environment distance to matrix
mantel(geo_dist, env_dist, method = "pearson", permutations = 999)  # Mantel test using Pearson's correlation

# Perform Mantel test for Fst vs geographical distance (IBD)
mantel(fst, geo_dist, method = "pearson", permutations = 999)  
# Perform Mantel test for Fst vs environmental distance (IBE)
mantel(fst, env_dist, method = "pearson", permutations = 999)  
# Partial Mantel test for IBE controlling for geographical distance
mantel.partial(fst, env_dist, geo_dist, method = "pearson", permutations = 999)  
# Partial Mantel test for IBD controlling for environmental distance
mantel.partial(fst, geo_dist, env_dist, method = "pearson", permutations = 999)  

library(ggplot2)
library(cowplot)
# Read the distance matrices for plotting
all_geo_dit = read.csv("geo_dist_panadapt.csv", head = FALSE)  
all_env_dist = read.csv("env_scale_panadapt.csv", head = FALSE)  
all_Fst = read.csv("Fst_matrix_neutral.csv", head = FALSE)  

# Define a function to transform the data into a suitable format for plotting
trans <- function(raw_data) {
    raw_data = raw_data[-1, -1]  # Remove the first row and column
    out_data = data.frame(raw_data[, 1])  # Initialize a new data frame with the first column
    colnames(out_data) = "value"  # Set column name to "value"
    
    # Loop to extract values from other columns and append them to the output data frame
    for (i in 2:43) {
        temp = data.frame(raw_data[i:43, i])  # Extract the i-th column
        colnames(temp) = "value"  # Set column name to "value"
        out_data = rbind(out_data, temp)  # Append to the output data frame
    }
    out_data = na.omit(out_data)  # Remove NA values from the data
    return(out_data)  # Return the transformed data
}

plot_data = as.data.frame(matrix(nrow = 946, ncol = 0))  
# Fill the data frame with transformed values from the distance matrices
plot_data$all_geo_dit = trans(all_geo_dit)$value  # Add geographical distance data
plot_data$all_env_dist = trans(all_env_dist)$value  # Add environmental distance data
plot_data$all_Fst = trans(all_Fst)$value  # Add Fst data
colnames(plot_data) = c("geo_dist", "env_dist", "all_fst")
write.csv(plot_data, file = "plot_data_neutral.csv", quote = FALSE, row.names = FALSE)  # Save plot data to CSV




###################################################################################################################
#########pan 
library(hierfstat) 
library(reshape2) 
library(adegenet)  
str <- read.structure("pan_adapt.str", n.ind = 201, n.loc = 16709, col.lab = 1, col.pop = 2, ask = FALSE, quiet = FALSE)  # Read STRUCTURE .str file
# Converting data from a STRUCTURE .str file to a genind object (genetic data format)
data <- genind2df(str)  # Convert to a data frame
write.csv(data, "pan_adapt_genind2df.csv")  
test <- genind2hierfstat(str)  # Convert to hierfstat format
fst <- genet.dist(test, method = "Ds")  # Compute genetic distance using Nei's standard genetic distance method

# Function to transform the distance matrix into a data frame and remove duplicate data
fun2 <- function(inDist) {
  df <- melt(as.matrix(inDist), varnames = c("row", "col"))  # Use melt from reshape2 to reshape the distance matrix
  df[as.numeric(df$row) > as.numeric(df$col), ]  # Keep only the upper triangle of the matrix to remove duplicates
}  

col_GD <- fun2(fst)  # Apply the function to the Fst matrix
write.csv(col_GD, "col_GD_panadapt.csv")  
Fst <- fst / (1 - fst)  # Compute Fst / (1 - Fst)
Fst_matrix <- as.data.frame(as.matrix(Fst))  # Convert genetic distance to Fst/(1-Fst) matrix and then to a data frame
write.csv(Fst_matrix, "Fst_matrix_panadapt.csv") 
library(vegan)  
library(geosphere)  
env = read.csv("WGR_env_6.csv", head = TRUE)  
ENV = as.data.frame(env[, c(4:9)])  # Extract columns 4 to 9 for environmental variables
ENV <- scale(ENV, center = TRUE, scale = TRUE)  # Standardize the environmental data
rownames(ENV) = env$POP  # Set row names to population names from the "POP" column
env_dist = vegdist(ENV, method = "euclidean", binary = FALSE, diag = FALSE, upper = FALSE, na.rm = FALSE)  # Compute Euclidean distance between environmental variables
write.csv(as.matrix(env_dist), file = "env_scale_panadapt.csv", quote = FALSE)  
fst = read.csv("Fst_matrix_panadapt.csv", head = TRUE, row.names = 1) 
# Calculate geographic distance
dist = as.data.frame(env[, c(2:3)])  # Extract the geographic coordinates columns as a data frame
dist <- scale(dist, center = TRUE, scale = TRUE)  # Standardize geographic distance data
rownames(dist) = env$POP  # Set row names to population names
geo_dist = vegdist(dist, method = "euclidean", binary = FALSE, diag = FALSE, upper = FALSE, na.rm = FALSE)  # Compute Euclidean distance between geographic coordinates
write.csv(as.matrix(geo_dist), file = "geo_dist_panadapt.csv", quote = FALSE) 
#### Geo-env distance Pearson correlation
env_dist = as.matrix(env_dist)  # Convert environmental distance to matrix
mantel(geo_dist, env_dist, method = "pearson", permutations = 999)  # Perform Mantel test (Geo-env correlation)

#### IBD (Isolation by Distance)
mantel(fst, geo_dist, method = "pearson", permutations = 999)  # Perform Mantel test for Fst vs geo_dist (IBD)
#### IBE (Isolation by Environment)
mantel(fst, env_dist, method = "pearson", permutations = 999)  # Perform Mantel test for Fst vs env_dist (IBE)
#### Partial Mantel Test for IBE (controlling for geographical distance)
mantel.partial(fst, env_dist, geo_dist, method = "pearson", permutations = 999)  # Perform partial Mantel test (IBE controlling for geo_dist)
#### Partial Mantel Test for IBD (controlling for environmental distance)
mantel.partial(fst, geo_dist, env_dist, method = "pearson", permutations = 999)  # Perform partial Mantel test (IBD controlling for env_dist)



library(ggplot2)  
library(cowplot)
all_geo_dit = read.csv("geo_dist_panadapt.csv", head = FALSE)  
all_env_dist = read.csv("env_scale_panadapt.csv", head = FALSE)  
all_Fst = read.csv("Fst_matrix_panadapt.csv", head = FALSE)  
# Define function to transform the raw distance data into a format suitable for plotting
trans <- function(raw_data) {
    raw_data = raw_data[-1, -1]  # Remove the first row and the first column
    out_data = data.frame(raw_data[, 1])  # Initialize the output data frame with the first column of data
    colnames(out_data) = "value"  # Set the column name to "value"
    
    # Loop to extract data from other columns
    for (i in 2:43) {
        temp = data.frame(raw_data[i:43, i])  # Extract the i-th column of data
        colnames(temp) = "value"  # Set the column name to "value"
        out_data = rbind(out_data, temp)  # Combine this data into the output data frame
    }
    
    out_data = na.omit(out_data)  # Remove missing values (NA)
    return(out_data)  # Return the transformed data
}

# Create an empty data frame to store the plot data
plot_data = as.data.frame(matrix(nrow = 946, ncol = 0))  # Initialize an empty data frame with 946 rows and 0 columns
# Populate the plot data with the transformed data
plot_data$all_geo_dit = trans(all_geo_dit)$value  # Add geographic distance data
plot_data$all_env_dist = trans(all_env_dist)$value  # Add environmental distance data
plot_data$all_Fst = trans(all_Fst)$value  # Add Fst data
# Set the column names for the plot data
colnames(plot_data) = c("geo_dist", "env_dist", "all_fst")
write.csv(plot_data, file = "plot_data_panadapt.csv", quote = FALSE, row.names = FALSE)



#########core_adapt
library(hierfstat)
library(reshape2)
library(adegenet)
str <- read.structure("core_adapt.str", n.ind = 201, n.loc = 1198, col.lab = 1, col.pop = 2, ask = FALSE, quiet = FALSE)
data <- genind2df(str)
write.csv(data,"core_adapt_genind2df.csv") 
# Convert to hierfstat format
test <- genind2hierfstat(str)
# Calculate genetic distance using Nei's standard method
fst <- genet.dist(test, method="Ds")
# Function to reshape the distance matrix and remove duplicates
fun2 <- function(inDist) {
  df <- melt(as.matrix(inDist), varnames = c("row", "col"))
  df[as.numeric(df$row) > as.numeric(df$col), ]
}  
col_GD <- fun2(fst)
write.csv(col_GD,"col_GD_coreadapt.csv") 
# Calculate Fst/(1-Fst)
Fst <- fst / (1 - fst)     
Fst_matrix <- as.data.frame(as.matrix(Fst))  # Convert to Fst/(1-Fst) matrix
write.csv(Fst_matrix,"Fst_matrix_coreadapt.csv") 
library(vegan)  
library(geosphere)
env = read.csv("WGR_env_6.csv", head = TRUE)  
ENV = as.data.frame(env[, c(4:9)])  
ENV <- scale(ENV, center = TRUE, scale = TRUE)  # Standardize the data
rownames(ENV) = env$POP  # Set row names from the POP column
env_dist = vegdist(ENV, method = "euclidean", binary = FALSE, diag = FALSE, upper = FALSE, na.rm = FALSE)  # Compute Euclidean distance for environmental variables
write.csv(as.matrix(env_dist), file = "env_scale_coreadapt.csv", quote = FALSE) 

# Read Fst matrix data
fst = read.csv("Fst_matrix_coreadapt.csv", head = TRUE, row.names = 1)  

# Calculate geographic distance
dist = as.data.frame(env[, c(2:3)])  # Extract geographic coordinates
dist <- scale(dist, center = TRUE, scale = TRUE)  # Standardize geographic data
rownames(dist) = env$POP  # Set row names from the POP column
geo_dist = vegdist(dist, method = "euclidean", binary = FALSE, diag = FALSE, upper = FALSE, na.rm = FALSE)  # Compute Euclidean distance for geographic variables
write.csv(as.matrix(geo_dist), file = "geo_dist_coreadapt.csv", quote = FALSE)  

# Pearson correlation for geo-env distance
env_dist = as.matrix(env_dist)
mantel(geo_dist, env_dist, method="pearson", permutations=999)  # Mantel test between geographic and environmental distances
# Pearson correlation for IBD (Isolation by Distance)
mantel(fst, geo_dist, method="pearson", permutations=999)  # Mantel test between Fst and geographic distance
# Pearson correlation for IBE (Isolation by Environment)
mantel(fst, env_dist, method="pearson", permutations=999)  # Mantel test between Fst and environmental distance
# Partial Mantel test for IBE
mantel.partial(fst, env_dist, geo_dist, method="pearson", permutations=999)  # Partial Mantel test between Fst and environmental distance controlling for geographic distance
# Partial Mantel test for IBD
mantel.partial(fst, geo_dist, env_dist, method="pearson", permutations=999)  # Partial Mantel test between Fst and geographic distance controlling for environmental distance


library(ggplot2)  
library(cowplot)  

all_geo_dit = read.csv("geo_dist_coreadapt.csv", head = FALSE)  
all_env_dist = read.csv("env_scale_coreadapt.csv", head = FALSE) 
all_Fst = read.csv("Fst_matrix_coreadapt.csv", head = FALSE) 

# Define function to transform data
trans <- function(raw_data) {
    raw_data = raw_data[-1, -1]  # Remove the first row and first column
    out_data = data.frame(raw_data[, 1])  # Initialize output data frame with the first column
    colnames(out_data) = "value"  # Set column name to "value"
    
    # Loop to extract data from other columns
    for (i in 2:43) {
        temp = data.frame(raw_data[i:43, i])  # Extract the i-th column of data
        colnames(temp) = "value"  # Set column name to "value"
        out_data = rbind(out_data, temp)  # Combine into output data frame
    }
    
    out_data = na.omit(out_data)  # Remove missing values
    return(out_data)  # Return transformed data
}

plot_data = as.data.frame(matrix(nrow = 946, ncol = 0))  # Initialize empty data frame

# Fill the plot data with the transformed data
plot_data$all_geo_dit = trans(all_geo_dit)$value  # Add geographic distance data
plot_data$all_env_dist = trans(all_env_dist)$value  # Add environmental distance data
plot_data$all_Fst = trans(all_Fst)$value  # Add Fst data
# Set column names for the plot data
colnames(plot_data) = c("geo_dist", "env_dist", "all_fst")
write.csv(plot_data, file = "plot_data_coreadapt.csv", quote = FALSE, row.names = FALSE)  

#########################################################################################################################################
## Adaptation (pan, core)

# Plot 1: Geographical distance vs Fst (adaptation)
p1 = ggplot(plot_data) +
  geom_point(aes(x = geo_dist, y = fst_adapt), 
             size = 3, alpha = 0.7, color = "#dc6f82", shape = 21, fill = "#fbf3d5", stroke = 0.8) +  # Scatter plot
  geom_smooth(aes(x = geo_dist, y = fst_adapt), 
              alpha = 0.7, formula = y ~ x, method = lm, se = TRUE, level = 0.95, color = "#d14365", 
              fill = "#f5e0ca", size = 1) +  # Linear regression line
  labs(x = "Geographical Distance (100km)", 
       y = expression(italic(F)[italic(ST)] / (1 - italic(F)[italic(ST)]))) +  # Axis labels
  panel_border(color = "black", size = 0.6) +  # Panel border
  theme_bw() +  # Black and white theme
  theme(text = element_text(family = "sans"),
        axis.ticks.length = unit(0.25, "lines"), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        plot.title = element_text(size = 15, hjust = 0),
        axis.title = element_text(size = 15), 
        panel.background = element_rect(fill = "white"), 
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "mm"))

# Plot 2: Environmental distance vs Fst (adaptation)
p2 = ggplot(plot_data) +
  geom_point(aes(x = env_dist, y = fst_adapt), 
             size = 3, alpha = 0.7, color = "#dc6f82", shape = 21, fill = "#fbf3d5", stroke = 0.8) +  # Scatter plot
  geom_smooth(aes(x = env_dist, y = fst_adapt), 
              alpha = 0.7, formula = y ~ x, method = lm, se = TRUE, level = 0.95, color = "#d14365", 
              fill = "#f5e0ca", size = 1) +  # Linear regression line
  labs(x = "Environment Distance", 
       y = expression(italic(F)[italic(ST)] / (1 - italic(F)[italic(ST)]))) +  # Axis labels
  panel_border(color = "black", size = 0.6) +  # Panel border
  theme_bw() +  # Black and white theme
  theme(text = element_text(family = "sans"),
        axis.ticks.length = unit(0.25, "lines"), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        plot.title = element_text(size = 15, hjust = 0),
        axis.title = element_text(size = 15), 
        panel.background = element_rect(fill = "white"), 
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "mm"))

## Neutral Fst
# Plot 1: Geographical distance vs Neutral Fst
p1 = ggplot(plot_data) +
  geom_point(aes(x = geo_dist, y = fst_neutral), 
             size = 3, alpha = 0.7, color = "#609ac6", shape = 21, fill = "#d6dde2", stroke = 0.8) +  # Scatter plot
  geom_smooth(aes(x = geo_dist, y = fst_neutral), 
              alpha = 0.7, formula = y ~ x, method = lm, se = TRUE, level = 0.95, color = "#4e4f52", 
              fill = "#d6d6d6", size = 1, fullrange = FALSE) +  # Linear regression line
  labs(x = "Geographical Distance (100km)", 
       y = expression(italic(F)[italic(ST)] / (1 - italic(F)[italic(ST)]))) +  # Axis labels
  panel_border(color = "black", size = 0.6, linetype = 1, remove = FALSE) +  # Panel border
  theme_bw() + 
  theme(text = element_text(family = "sans"),
        axis.ticks.length = unit(0.25, "lines"), 
        axis.ticks = element_line(colour = "black", unit(0.6, "line")),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"), 
        plot.title = element_text(size = 15L, hjust = 0),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "mm"))

# Plot 2: Environmental distance vs Neutral Fst
p2 = ggplot(plot_data) +
  geom_point(aes(x = env_dist, y = fst_neutral), 
             size = 3, alpha = 0.7, color = "#609ac6", shape = 21, fill = "#d6dde2", stroke = 0.8) +  # Scatter plot
  geom_smooth(aes(x = env_dist, y = fst_neutral), 
              alpha = 0.7, formula = y ~ x, method = lm, se = TRUE, level = 0.95, color = "#4e4f52", 
              fill = "#d6d6d6", size = 1, fullrange = FALSE) +  # Linear regression line
  labs(x = "Environment Distance", 
       y = expression(italic(F)[italic(ST)] / (1 - italic(F)[italic(ST)]))) +  # Axis labels
  panel_border(color = "black", size = 0.6, linetype = 1, remove = FALSE) +  # Panel border
  theme_bw() +  
  theme(text = element_text(family = "sans"),
        axis.ticks.length = unit(0.25, "lines"), 
        axis.ticks = element_line(colour = "black", unit(0.6, "line")),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"), 
        plot.title = element_text(size = 15L, hjust = 0),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "mm"))
