###pan_loci DAPC+PCA
library(ggplot2)
library(adegenet)
# Read the structure file
Rheum <- read.structure("pan_adapt.str", n.ind = 201, n.loc = 16709, col.lab = 1, col.pop = 2, ask = FALSE, quiet = FALSE)
# Perform clustering analysis
results <- find.clusters(Rheum, max.n.clust = 60, n.pca = 100)  # Set the maximum number of clusters and use 100 PCs
DAPC1 <- dapc(Rheum, results$grp)  # Perform DAPC analysis
# Save the cluster group information to a CSV file
grp <- results$grp
write.csv(grp, "grp.csv")
# Extract cluster assignment and get unique cluster labels
clusters <- results$grp
unique_clusters <- unique(clusters)
# Define colors for each cluster
myCols <- c("olivedrab3", "lightslateblue", "coral1")  # Colors for cluster 1, 2, and 3
# Assign colors to points based on their cluster
point_colors <- myCols[as.numeric(as.factor(clusters))]
# Create a color plot for DAPC
colorplot(DAPC1$ind, DAPC1$ind, transp = TRUE, cex = 3, xlab = "DAPC 1", ylab = "DAPC 2")
abline(v = 0, h = 0, col = "grey", lty = 2)  # Add center lines at 0,0
# Allow points to extend outside the plot area
par(xpd = TRUE)
# Plot the points with the assigned colors for each cluster
points(DAPC1$ind.coord[, 1], DAPC1$ind.coord[, 2], 
       pch = 21, cex = 3, col = "black", bg = point_colors, lwd = 1.5)
# Add a legend to the plot
legend("topright", 
       legend = c("East", "West1", "West2"),  # Labels for the clusters
       pch = 21, 
       pt.bg = myCols,  # Set the background color of the legend points
       title = "Clusters", 
       cex = 0.8, 
       bty = "n")  # No box around the legend



plink --threads 16 --bfile pan_adapt --pca 3 --out pan_adapt_pca --allow-extra-chr
##This step produces two files, one ending in .eigenval, which records the eigenvalues, which are used to calculate the weight of each PC. The other is a file ending in .eigenvec that records the eigenvectors
##Where .eigenval represents the weighting of each PCA and the other records the eigenvectors, which are used to plot the axes
104.815
42.6822
31.8448

 PC1=58.44
 PC2=23.80

library(ggplot2)
#install.packages("scatterplot3d")
library(scatterplot3d)
library(ggrepel)
# Load necessary libraries
library(ggplot2)      # For data visualization
library(scatterplot3d) # For 3D plotting (not used here, but can be helpful for 3D plots)
library(ggrepel)      # For text labels that don't overlap
# Set colors for the different groups
color <- c("olivedrab3", "lightslateblue", "coral1")
# Load the PCA results data
df = read.csv("pan_pca.csv", header = TRUE)
# Convert the 'group' column to a factor (for categorical color mapping)
df$group = factor(df$group)
# Assign colors to groups based on the 'group' column
colors <- color[as.numeric(df$group)]
# Create a basic scatter plot using ggplot2
p1 <- ggplot(data = df, mapping = aes(PC1, PC2, color = colors))
# Add points to the plot, set point size and color
p2 <- p1 + geom_point(color = colors, size = 3.5) +
  geom_text_repel(label = "", color = colors)  # Add labels (currently empty, but can be customized)
# Customize plot labels and theme
p2 + xlab("PC1(58.44%)") + ylab("PC2(23.80%)") +
  theme(text = element_text(size = 16))  # Increase font size for better readability



#######
##CORE_loci DAPC+PCA

library(ggplot2)
library(adegenet)
# Read the structure file containing genetic data
Rheum <-  read.structure("core_adapt.str", n.ind = 201, n.loc = 1198, col.lab = 1, col.pop = 2, ask = FALSE, quiet = FALSE)
# Perform clustering analysis based on PCA results
results <- find.clusters(Rheum, max.n.clust = 60, n.pca = 100)  # Maximum 60 clusters, PCA with 100 components
# Perform Discriminant Analysis of Principal Components (DAPC)
DAPC1 <- dapc(Rheum, results$grp)  # Perform DAPC using the clusters found
# Save the cluster group information into a CSV file
grp <- results$grp
write.csv(grp, "grp.csv")  # Save cluster information into a file
# Extract unique clusters for further use
clusters <- results$grp
unique_clusters <- unique(clusters)
# Define colors for each cluster
myCols <- c("lightslateblue", "coral1", "olivedrab3")
# Assign colors based on the cluster labels
point_colors <- myCols[as.numeric(as.factor(clusters))]
# Create a color plot for the DAPC results
colorplot(DAPC1$ind, DAPC1$ind, transp = TRUE, cex = 3, xlab = "DAPC 1", ylab = "DAPC 2")
abline(v = 0, h = 0, col = "grey", lty = 2)  # Add center lines at 0,0
# Allow points to exceed the plot region and add the points with assigned colors
par(xpd = TRUE)  # Allow plotting outside the default area
points(DAPC1$ind.coord[, 1], DAPC1$ind.coord[, 2], 
       pch = 21, cex = 3, col = "black", bg = point_colors, lwd = 1.5)
# Add a legend to the plot
legend("topright", 
       legend = c("East", "West1", "West2"),  # Labels for clusters
       pch = 21, 
       pt.bg = myCols,  # Background colors for legend items
       title = "Clusters", 
       cex = 0.8, 
       bty = "n")  # No border around the legend

	   
plink --threads 16 --bfile core_adapt --pca 3 --out core_adapt_pca --allow-extra-chr
256.298
21.1057
8.88698
 PC1=89.52
 PC2=7.37

# Load ggplot2 and ggrepel for creating the PCA plot
library(ggplot2)
library(ggrepel)

# Load PCA data from CSV file
df = read.csv("core_pca.csv", header = T)

# Define the color mapping based on the groups
color <- c("lightslateblue", "coral1", "olivedrab3")
df$group = factor(df$group)
colors <- color[as.numeric(df$group)]

# Plot the PCA results
p1 <- ggplot(data = df, mapping = aes(PC1, PC2, color = colors))
p2 <- p1 + geom_point(color = colors, size = 3.5) + geom_text_repel(label = "", color = colors)
p2 + xlab("PC1(89.52%)") + ylab("PC2(7.37%)") + theme(text = element_text(size = 16))

 
 
 
 ##all
 ##all_snp pca
plink --threads 16 --bfile 213snp7 --pca 3 --out 201snp_pca --allow-extra-chr
 PC1=46.46
 PC2=34.33

## Principal Component Analysis (PCA) and Visualization in R
# Load required libraries
library(ggplot2)  # For plotting
library(scatterplot3d)  # For 3D scatter plot (not used directly here, but loaded)
library(ggrepel)  # For label positioning in plots
# Read the PCA data
df <- read.csv("201snp_pca.csv", header = TRUE)
# Define colors for different groups
color <- c("coral1", "lightslateblue", "olivedrab3")
# Convert the 'group' column to a factor for proper grouping
df$group <- factor(df$group)

# Create a PCA plot using ggplot2
p1 <- ggplot(data = df, mapping = aes(x = PC1, y = PC2, color = group, label = name)) +
  geom_point(size = 3.5) +  # Plot points with size 3.5
  scale_color_manual(values = color,  # Manually set colors for groups
                     labels = c("1" = "East", "2" = "West1", "3" = "West2")) +  # Label groups
  xlab("PC1(46.46%)") +  # X-axis label with variance explained
  ylab("PC2(34.33%)") +  # Y-axis label with variance explained
  theme(text = element_text(size = 16))  # Increase text size for readability
print(p1)


## Discriminant Analysis of Principal Components (DAPC)
# Load necessary libraries
library(adegenet)  # For genetic data analysis
# Read the structure file containing genetic data
Rheum <- read.structure("201snpfinal_2.str", n.ind = 201, n.loc = 311084, col.lab = 1, col.pop = 2, ask = FALSE, quiet = FALSE)
# Perform clustering based on PCA
results <- find.clusters(Rheum, max.n.clust = 60, n.pca = 100)
# Perform DAPC (Discriminant Analysis of Principal Components)
DAPC1 <- dapc(Rheum, results$grp)
# Save the clustering results
grp <- results$grp
write.csv(grp, "grp.csv")
# Extract cluster information
clusters <- results$grp
# Define custom colors for each cluster
unique_clusters <- unique(clusters)
myCols <- c("lightslateblue", "coral1", "olivedrab3")
# Assign colors to points based on their cluster
point_colors <- myCols[as.factor(clusters)]  # Mapping cluster numbers to colors
# Create the color plot for DAPC results
colorplot(DAPC1$ind, DAPC1$ind, transp = TRUE, cex = 3, xlab = "DAPC 1", ylab = "DAPC 2")
abline(v = 0, h = 0, col = "grey", lty = 2)  # Add grey dashed lines at (0,0)
# Adjust plotting parameters to allow points outside the plot area
par(xpd = TRUE)
# Add points to the plot with color coding
points(DAPC1$ind.coord[, 1], DAPC1$ind.coord[, 2], 
       pch = 21, cex = 3, col = "black", bg = point_colors, lwd = 1.5)

