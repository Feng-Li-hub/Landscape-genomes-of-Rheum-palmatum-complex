————————————————————————  SNP raw data file generation —————————————————————————
# Using plink software to convert the data into .raw format
plink --file 213snp3 --recodeA --out mydata --allow-extra-chr
# Generated files: mydata.log, mydata.nosex, mydata.raw
# The following files are generated after the conversion:
# mydata.log: A log file containing the details of the conversion process.
# mydata.nosex: A file indicating individuals with missing sex information.
# mydata.raw: The raw genotype data in tab-delimited format, which will be used for further analysis.
# Next, we will work with the mydata.raw file to create a raw.csv file and remove columns 2 to 6.

#### wc -l 213snp3.map 213snp3.ped
# The 'wc -l' command counts the number of lines in the specified files.
#   1432892 213snp3.map
#      201 213snp3.ped
#   1433093 total
# This output indicates the number of lines in the map and ped files used for SNP analysis.





———————————————————————— lfmm_lea ———————————————————————— 
library(vegan)  # Load vegan package for ecological and environmental data analysis
library(LEA)    # Load LEA package for latent factor model analysis

# Reading environmental data from CSV
env <- read.csv("env_6.csv")
pred <- env[,5:10]
# Perform PCA using redundancy analysis (RDA)
pred.pca <- rda(pred, scale=T)  # Scale the predictor data before performing PCA
summary(pred.pca)$cont  # Show summary of PCA results (contingency table)

# Plot the screeplot of the PCA to visualize the variance explained by each principal component
screeplot(pred.pca, main = "Screeplot of Rheum Predictor Variables with Broken Stick", 
          bstick=TRUE, type="barplot")

# Calculate and save the scores for the first six principal components
correction <- round(scores(pred.pca, choices=1:6, display="species", scaling=0), digits=3)
write.csv(correction, "correction.csv")  # Save the scores to a CSV file

# Get the scores for the first three principal components and save to a CSV
pred.PC1_3 <- scores(pred.pca, choices=c(1, 2, 3), display="sites", scaling=0)
write.csv(pred.PC1_3, "pred.PC1_3.csv")  # Save PC1, PC2, PC3 scores

# File processing: convert commas to tabs for LFMM input file
sed 's/,/\t/g' 201.lfmm > 201_tab.lfmm  ## Convert CSV to tab-separated format

# Run LFMM (Latent Factor Mixed Modeling) on the prepared data
project = lfmm("201_tab.lfmm", "clima2.env", K = 3, repetitions = 5, CPU = 8, 
               iterations = 1000, burnin = 500, project = "new")  ## Perform LFMM analysis with K=3 latent factors


setwd("C:\\Users\\feng\\Desktop\\景观\\LEA_lfmm")
# Load the LFMM project from the specified directory
project <- load.lfmmProject("201_tab_clima2.lfmmProject")

# Extract z-scores for the first principal component (PC1) and calculate median values
z.pc1 <- z.scores(project, K = 3, d = 1)
z.pc1 <- apply(z.pc1, 1, median)  ## Take median of repeated z-scores for PC1
lambda.pc1 = median(z.pc1^2) / qchisq(0.5, df = 1)  ## Calculate lambda for PC1
lambda.pc1
#[1] 9.556868

# Compute the adjusted p-values for PC1 based on chi-squared distribution
p.pc1.adj = pchisq(z.pc1^2 / lambda.pc1, df = 1, lower = FALSE)

# Repeat the process for the second principal component (PC2)
z.pc2 <- z.scores(project, K = 3, d = 2)  ## Extract z-scores for PC2
z.pc2 <- apply(z.pc2, 1, median)  ## Take median of repeated z-scores for PC2
lambda.pc2 = median(z.pc2^2) / qchisq(0.5, df = 1)  ## Calculate lambda for PC2
lambda.pc2
#[1] 6.016021

# Compute the adjusted p-values for PC2
p.pc2.adj = pchisq(z.pc2^2 / lambda.pc2, df = 1, lower = FALSE)

# Repeat the process for the third principal component (PC3)
z.pc3 <- z.scores(project, K = 3, d = 3)  ## Extract z-scores for PC3
z.pc3 <- apply(z.pc3, 1, median)  ## Take median of repeated z-scores for PC3
lambda.pc3 = median(z.pc3^2) / qchisq(0.5, df = 1)  ## Calculate lambda for PC3
lambda.pc3
#[1] 2.031029

# Compute the adjusted p-values for PC3
p.pc3.adj = pchisq(z.pc3^2 / lambda.pc3, df = 1, lower = FALSE)

# Plot histograms for the adjusted p-values of PC1, PC2, and PC3
hist(p.pc1.adj, col = "blue", main = "pc1", xlab='')  ## Histogram for adjusted p-values of PC1
hist(p.pc2.adj, col = "blue", main = "pc2", xlab='')  ## Histogram for adjusted p-values of PC2
hist(p.pc3.adj, col = "blue", main = "pc3", xlab='')  ## Histogram for adjusted p-values of PC3

# Adjust the p-values using the q-value method to control for false discovery rate (FDR)
library(qvalue)
q.pc1 <- qvalue(p.pc1.adj)$qvalues  ## Adjusted q-values for PC1
q.pc2 <- qvalue(p.pc2.adj)$qvalues  ## Adjusted q-values for PC2
q.pc3 <- qvalue(p.pc3.adj)$qvalues  ## Adjusted q-values for PC3

# Filter significant results where q-value < 0.01 and |z| > 2
sum(q.pc1 < 0.01 & abs(z.pc1) > 2)  ## Count significant results for PC1
sum(q.pc2 < 0.01 & abs(z.pc2) > 2)  ## Count significant results for PC2
sum(q.pc3 < 0.01 & abs(z.pc3) > 2)  ## Count significant results for PC3


# Combine the z-scores and p-values from the three principal components (PCs)
lfmm.results <- cbind(z.pc1, q.pc1, z.pc2, q.pc2, z.pc3, q.pc3)
write.csv(lfmm.results, "lfmm.results.csv")

# Count how many rows meet the following conditions:
# 1) q.pc1 < 0.01 and |z.pc1| > 2
# 2) q.pc2 < 0.01 and |z.pc2| > 2
# 3) q.pc3 < 0.01 and |z.pc3| > 2
sum(q.pc1<0.01 & abs(z.pc1)>2 | q.pc2<0.01 & abs(z.pc2)>2 | q.pc3<0.01 & abs(z.pc3)>2)


load("C:/Users/feng/Desktop/景观/lfmm+RDA/lfmmpc.rdata")
gen <- gen[,-1]

# Add column names to the results for each PC (to keep track of the loci names)
lfmm_pc1.results <- cbind(colnames(gen[1,]), z.pc1, q.pc1)
lfmm_pc2.results <- cbind(colnames(gen[1,]), z.pc2, q.pc2)
lfmm_pc3.results <- cbind(colnames(gen[1,]), z.pc3, q.pc3)

# Save the results for each PC as separate CSV files
write.csv(lfmm_pc1.results, "lfmm_pc1.results.csv")
write.csv(lfmm_pc2.results, "lfmm_pc2.results.csv")
write.csv(lfmm_pc3.results, "lfmm_pc3.results.csv")


lfmm_pc1.results <- read.csv("lfmm_pc1.results.csv", header = TRUE)
lfmm_pc2.results <- read.csv("lfmm_pc2.results.csv", header = TRUE)
lfmm_pc3.results <- read.csv("lfmm_pc3.results.csv", header = TRUE)


library(dplyr)
# Filter data for PC1 based on the conditions (q.pc1 < 0.01 and |z.pc1| > 2)
filtered_data_1 <- lfmm_pc1.results %>%
  filter( lfmm_pc1.results$q.pc1 < 0.01, abs(lfmm_pc1.results$z.pc1) > 2)
write.csv(filtered_data_1, "filtered_lfmm_pc1.results.csv", row.names = FALSE)

# Filter data for PC2 based on the same conditions
filtered_data_2 <- lfmm_pc2.results %>%
  filter( lfmm_pc2.results$q.pc2 < 0.01, abs(lfmm_pc2.results$z.pc2) > 2)
write.csv(filtered_data_2, "filtered_lfmm_pc2.results.csv", row.names = FALSE)

# Filter data for PC3 based on the same conditions
filtered_data_3 <- lfmm_pc3.results %>%
  filter( lfmm_pc3.results$q.pc3 < 0.01, abs(lfmm_pc3.results$z.pc3) > 2)
write.csv(filtered_data_3, "filtered_lfmm_pc3.results.csv", row.names = FALSE)
# Combine the filtered results from all three PCs into one data frame
combined_data <- bind_rows(filtered_data_1, filtered_data_2, filtered_data_3)
# Extract unique loci from the combined data
combined_data2 <- unique(combined_data$X.1)
write.csv(combined_data2, "lfmm_loci.csv")


————————————RDA——————————————————————

# Load genotype data and replace missing values
library(data.table)
gen <- fread("201_no_NA.csv") # If the first column is missing headers, use row.names=1, otherwise read allele frequency data for all populations
dim(gen) # Check the dimensions of the data
sum(is.na(gen)) # Check how many genotypes are missing

# Load environmental data
env <- read.csv("env_6.csv")
str(env) # Check the structure of the environmental data
env$pop <- as.character(env$pop) # Convert the population column to character type
env$Group <- as.factor(env$Group) # Convert the Group column to factor type
identical(rownames(gen), env[,1]) # Ensure that row names in gen match the first column in env

# Load necessary libraries for analysis
library(vegan)  # For ecological statistics and multivariate analysis
library(lfmm) 
library(psych)

# Perform correlation analysis on environmental data and save the results
pdf("pca2.pdf") 
pairs.panels(env[,5:10], scale=T)    # Columns 5-10 are environmental variables based on GF predictions, remove variables with correlation > 0.8, keeping 6 environmental variables
dev.off()

# Define the environmental predictors
pred <- env[,5:10]

# Perform Redundancy Analysis (RDA)
Rheum.rda <- rda(gen[, -1] ~ ., data=pred, scale=T)
Rheum.rda # Display the RDA model results

# Model summary
# Inertia represents the total variation explained by the model
# Constrained inertia represents the variation explained by the environmental variables
# Unconstrained inertia represents the unexplained variation

# Eigenvalues for constrained and unconstrained axes
# Eigenvalues represent the variance explained by each axis

r2 <- RsquareAdj(Rheum.rda)  # Calculate R-squared for RDA
r2 # Display R-squared values

rda_adj <- r2$adj.r.squared  # Adjusted R-squared
rda_adj # Display adjusted R-squared value

# Summary of eigenvalues for constrained axes
summary(eigenvals(Rheum.rda, model = "constrained")) # Eigenvalue summary for constrained axes

# Create a scree plot to visualize the eigenvalues of each constrained axis
screeplot(Rheum.rda)

# Perform permutation test for the RDA model to assess the significance
signif.full <- anova.cca(Rheum.rda, parallel=getOption("mc.cores"))
signif.full # Display permutation test results

# Calculate Variance Inflation Factor (VIF) for each environmental predictor to check for multicollinearity
vif.cca(Rheum.rda)

# Identify selected loci based on RDA scores
load.rda <- scores(Rheum.rda, choices=c(1:3), display="species")  # Extract SNP scores from the first three constrained axes

# Plot histograms for loadings on RDA axes 1, 2, and 3
hist(load.rda[,1], main="Loadings on RDA1", col="orange")
hist(load.rda[,2], main="Loadings on RDA2", col="orange")
hist(load.rda[,3], main="Loadings on RDA3", col="orange") 

# Define a function to identify outliers (loci with extreme values in loadings)
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # Find outliers based on z standard deviations from the mean
  x[x < lims[1] | x > lims[2]]               # Return loci names that are outliers
}

# Identify outliers (selected loci) in each of the first three RDA axes
cand1 <- outliers(load.rda[,1],2.7) # Identify selected loci in RDA1 (0 loci)
cand2 <- outliers(load.rda[,2],2.7) # Identify selected loci in RDA2 (4929 loci)
cand3 <- outliers(load.rda[,3],2.7) # Identify selected loci in RDA3 (1990 loci)
ncand <- length(cand1) + length(cand2) + length(cand3) # Total number of selected loci across all three RDA axes

# Combine identified selected loci into one data frame
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)  # Convert SNP names to character type
str(cand) # Check the structure of the selected loci data

# For each selected SNP, calculate the correlation with the 6 environmental predictors
foo <- matrix(nrow=(ncand), ncol=6)  # Create a matrix with 6 columns for the 6 predictors
colnames(foo) <- c("BIO02","BIO03","BIO04","BIO10","BIO12","BIO15")
for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]  # Get the SNP name
  snp.gen <- gen[,nam]  # Extract SNP genotype data
  foo[i,] <- apply(pred,2,function(x) cor(x,snp.gen))  # Calculate correlations with each environmental predictor
}

# Add correlation results to the selected loci data
cand <- cbind.data.frame(cand,foo)  
head(cand) # Display the first few rows of the updated data

# Write the results to a CSV file, including the correlations with environmental predictors
length(cand$snp[duplicated(cand$snp)])  # Check how many duplicated SNPs there are
cand <- cand[!duplicated(cand$snp),] # Remove duplicate SNPs
write.csv(cand, "adaptive loci.csv")  # Save the final selected loci to a CSV file

# Identify the strongest predictor for each SNP and save the results
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,10] <- names(which.max(abs(bar[4:9]))) # Assign the strongest predictor to the 10th column
  cand[i,11] <- max(abs(bar[4:9]))              # Assign the strongest correlation to the 11th column
}

# Rename the columns for clarity
colnames(cand)[10] <- "predictor"
colnames(cand)[11] <- "correlation"

# Count how many selected loci correspond to each predictor
table(cand$predictor)  ## Count the number of selected loci for each predictor
write.csv(cand, "adaptive loci_most strong predictor.csv")  






————————————————————venn————————————————————
##########################################################
library(VennDiagram)   # For creating Venn diagrams
library(readr)         # For reading CSV files

# Read the LFMM and RDA data
lfmm_data <- read.csv("lfmm_loci.csv", header = FALSE)  # Read the LFMM data (without headers)
rda_data <- read.csv("RDA_2.7.csv", header = FALSE)    # Read the RDA data (without headers)

# Extract points (loci names) from both LFMM and RDA datasets
lfmm_points <- lfmm_data[[1]]  
rda_points <- rda_data[[1]]   

# Create the Venn diagram to show the overlap between LFMM and RDA loci
venn.plot <- venn.diagram(
  x = list(LFMM = lfmm_points, RDA = rda_points),  
  category.names = c("LFMM", "RDA"),               
  filename = NULL,                                 
  fill = c("skyblue", "pink"),                     
  alpha = 0.5,                                    
  cex = 2,                                         
  cat.cex = 2,                                    
  cat.pos = c(-20, 20),                            
  main = "Venn Diagram of LFMM and RDA Points"     
)

# Save the Venn diagram as a PDF file
pdf("venn_diagram.pdf")   
grid.draw(venn.plot)      
dev.off()                 

# Find the intersection of loci between LFMM and RDA
intersection_points <- intersect(lfmm_points, rda_points)  # Get the common loci between LFMM and RDA
write.csv(intersection_points, "intersection_points.csv", row.names = FALSE)  


























