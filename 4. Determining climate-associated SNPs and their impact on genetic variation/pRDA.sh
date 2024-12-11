library(dplyr)
# Read the 213snp_LD.prune.in file
ld_prune_data <- read.table("213snp_LD.prune.in", header = FALSE, stringsAsFactors = FALSE)
# Read the pan_loci2.csv file
pan_loci_data <- read.csv("pan_loci2.csv", header = FALSE, stringsAsFactors = FALSE)
# Remove loci from ld_prune_data that are present in pan_loci_data
remaining_ld_prune_data <- ld_prune_data[!ld_prune_data$V1 %in% pan_loci_data$V1, ]
# Read the remaining LD prune data (in case you want to overwrite or use another file)
remaining_ld_prune_data  <- read.csv("remaining_ld_prune_data.csv")
# Randomly sample 200,000 loci from the remaining LD prune data
random_samples <- remaining_ld_prune_data %>% sample_n(200000, replace = FALSE)
write.csv(random_samples, "random_samples_200000.csv")
# Load necessary libraries for PCA analysis
library(vegan)
library(dplyr)
library(data.table)
# Read the allele frequency data
allele_freq <- fread("F:\\王聪颖\\景观_2\\AF\\allele_freq_AF2.csv")
# Read the random samples
netural <- fread("random_samples_200000.csv", header = FALSE)  # Assuming no column names, only locus names
# Extract locus names from the data
snp_names <- netural$V1
# Filter out the ID column and select the columns related to the loci names
columns_to_extract <- c(snp_names)  # Includes the first column and the loci columns of interest
filtered_data_core <- allele_freq %>% select(all_of(columns_to_extract))
# Save the filtered allele frequency data to a new file
fwrite(filtered_data_core, "netural_loci2_AF.csv", row.names = FALSE)

##########
### Perform PCA analysis on the neutral genetic data to get the PCs.
# Read the allele frequency data for neutral loci
Freq_neutral <- fread("netural_loci2_AF.csv", head = TRUE)
# Perform PCA on the neutral genetic loci (excluding the first column)
pca <- rda(Freq_neutral[,-1], scale = TRUE)
# Plot the scree plot to show the eigenvalues for the first 10 PCs
screeplot(pca, type = "barplot", npcs = 10, main = "PCA Eigenvalues")
# Extract the first 3 principal components (PCs)
PCs <- scores(pca, choices = c(1:3), display = "sites", scaling = 0)
# Create a data frame with population names and their corresponding PC scores
PopStruct <- data.frame(POP = Freq_neutral[, 1], PCs)
# Rename columns to PC1, PC2, and PC3
colnames(PopStruct) <- c("POP", "PC1", "PC2", "PC3")
# Save the population structure with the first 3 PCs to a CSV file
write.csv(PopStruct, file = "PopStruct.csv", quote = FALSE)



##################################################################################
#####pan
# Load geographic data (longitude and latitude)
geog = read.csv("long_lat.csv", head = TRUE)
# Load environmental data (variables such as BIO02, BIO03, etc.)
env = read.csv("env_6.csv", head = TRUE)
ENV = as.data.frame(env[, c(2:7)])
# Standardize the environmental variables (mean = 0, sd = 1)
ENV <- scale(ENV, center = TRUE, scale = TRUE)
write.csv(ENV, file = "env.csv", quote = FALSE)
# Load allele frequency data
AllFreq <- fread("pan_loci2_AF.csv")
### Combine all variables into one data frame
Variables <- data.frame(geog, PopStruct[, -1], ENV)
write.csv(Variables, file = "Variables.csv", quote = FALSE)  # Output all variables again

### Variance partitioning: Decompose the factors driving genetic variation
# Full model 
pRDAfull <- rda(AllFreq[, -1] ~ PC1 + PC2 + PC3 + longitude + latitude + BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15, Variables)

# Adjusted R-squared value for the full model
RsquareAdj(pRDAfull)

# Output R-squared and adjusted R-squared for the full model
$r.squared  # [1] 0.7833547
$adj.r.squared  # [1] 0.7064806

# Perform ANOVA on the full model to test the significance
anova(pRDAfull)

# Permutation test for the RDA under reduced model
# Number of permutations = 999
# The full model has a significant result (p-value < 0.001)
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = AllFreq[, -1] ~ PC1 + PC2 + PC3 + longitude + latitude + BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15, data = Variables)
         Df Variance     F Pr(>F)    
Model    11  1522.31 10.19  0.001 ***
Residual 31   421.01                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Climatic model
pRDAclim <- rda(AllFreq[, -1] ~ BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15 + Condition(longitude + latitude + PC1 + PC2 + PC3), Variables)

# Adjusted R-squared value for the climatic model
RsquareAdj(pRDAclim)

# Output R-squared and adjusted R-squared for the climatic model
$r.squared  # [1] 0.09782426
$adj.r.squared  # [1] 0.06344601

# Perform ANOVA on the climatic model to test the significance
anova(pRDAclim)

# Permutation test for the RDA under the reduced model
# The model is significant with a p-value < 0.001
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = AllFreq[, -1] ~ BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15 + Condition(longitude + latitude + PC1 + PC2 + PC3), data = Variables)
         Df Variance     F Pr(>F)    
Model     6   190.10 2.333  0.001 ***
Residual 31   421.01                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Structural model
pRDAstruct <- rda(AllFreq[, -1] ~ PC1 + PC2 + PC3 + Condition(longitude + latitude + BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15), Variables)

# Adjusted R-squared value for the structural model
RsquareAdj(pRDAstruct)

# Output R-squared and adjusted R-squared for the structural model
$r.squared  # [1] 0.1616813
$adj.r.squared  # [1] 0.1738251

# Perform ANOVA on the structural model to test the significance
anova(pRDAstruct)

# Permutation test for the RDA under the reduced model
# The model is significant with a p-value < 0.001
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = AllFreq[, -1] ~ PC1 + PC2 + PC3 + Condition(longitude + latitude + BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15), data = Variables)
         Df Variance      F Pr(>F)    
Model     3   314.20 7.7117  0.001 ***
Residual 31   421.01                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Pure geography model
pRDAgeog <- rda(AllFreq[, -1] ~ longitude + latitude + Condition(BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15 + PC1 + PC2 + PC3), Variables)

# Adjusted R-squared value for the geography-only model
RsquareAdj(pRDAgeog)

# Output R-squared and adjusted R-squared for the geography-only model
$r.squared  # [1] 0.02318248
$adj.r.squared  # [1] 0.01171592

# Perform ANOVA on the geography-only model to test the significance
anova(pRDAgeog)

# Permutation test for the RDA under the reduced model
# The model is significant with a p-value = 0.021
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = AllFreq[, -1] ~ longitude + latitude + Condition(BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15 + PC1 + PC2 + PC3), data = Variables)
         Df Variance      F Pr(>F)  
Model     2    45.05 1.6586  0.021 *
Residual 31   421.01                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1





########################################################################################################
##########pan west
# Load geographic data (longitude and latitude for 28 populations)
geog = read.csv("28pop_longlat.csv", head = TRUE)
# Load allele frequency data for all populations
AllFreq = fread("pan_west_AF.csv")
# Load neutral loci frequency data
Freq_neutral = fread("F:\\王聪颖\\景观_2\\中性位点Prda\\netural_west_AF.csv", head = TRUE)
# Load environmental variables data
env = read.csv("Env_west.csv", head = TRUE)
# Select relevant columns (environmental variables) and scale them (standardization: mean = 0, sd = 1)
ENV = as.data.frame(env[, c(2:7)])
ENV <- scale(ENV, center = TRUE, scale = TRUE)
write.csv(ENV, file = "env.csv", quote = FALSE)


library(vegan)
### Perform PCA on the neutral genetic loci
# RDA is used here to perform PCA (Principal Component Analysis) for scaling the neutral loci data
pca <- rda(Freq_neutral[, -1], scale = TRUE)
# Create a scree plot to visualize eigenvalues from the PCA
screeplot(pca, type = "barplot", npcs = 10, main = "PCA Eigenvalues")
### Extract the first three principal components from PCA
PCs <- scores(pca, choices = c(1:3), display = "sites", scaling = 0
# Create a data frame with population names and their corresponding first three PCs
PopStruct <- data.frame(POP = Freq_neutral[, 1], PCs)
# Rename columns for better clarity (POP, PC1, PC2, PC3)
colnames(PopStruct) <- c("POP", "PC1", "PC2", "PC3")
write.csv(PopStruct, file = "PopStruct.csv", quote = FALSE)
### Combine all variables into one data frame (geographic, PCA, and environmental data)
Variables <- data.frame(geog, PopStruct[, -1], ENV)
write.csv(Variables, file = "Variables.csv", quote = FALSE)

### Variance Partitioning: Decompose the factors influencing genetic variation
# Perform redundancy analysis (RDA) to analyze the impact of environmental variables, geographic coordinates, and PCs on allele frequencies
pRDAfull <- rda(AllFreq[, -1] ~ PC1 + PC2 + PC3 + longitude + latitude + BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15, Variables)
# Calculate the adjusted R-squared value for the full model
RsquareAdj(pRDAfull)
# Output the R-squared and adjusted R-squared values for the full model
$r.squared  # [1] 0.8805174
$adj.r.squared  # [1] 0.7983732
# Perform ANOVA for the full model to test the significance of each factor
anova(pRDAfull)
# Permutation test for the RDA under the reduced model (999 permutations)
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = AllFreq[, -1] ~ PC1 + PC2 + PC3 + longitude + latitude + BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15, data = Variables)
         Df Variance      F Pr(>F)    
Model    11  1694.14 10.719  0.001 ***
Residual 16   229.89                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Climatic Model: Include only the environmental variables, conditioned on geographic and PCA variables
pRDAclim <- rda(AllFreq[, -1] ~ BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15 + Condition(longitude + latitude + PC1 + PC2 + PC3), Variables)
# Calculate the adjusted R-squared for the climatic model
RsquareAdj(pRDAclim)

# Output the R-squared and adjusted R-squared for the climatic model
$r.squared  # [1] 0.06487269
$adj.r.squared  # [1] 0.02462735

# Perform ANOVA on the climatic model to test the significance of the climatic variables
anova(pRDAclim)

# Permutation test for the RDA under the reduced model
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = AllFreq[, -1] ~ BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15 + Condition(longitude + latitude + PC1 + PC2 + PC3), data = Variables)
         Df Variance      F Pr(>F)   
Model     6   124.82 1.4479  0.002 **
Residual 16   229.89                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Structural Model: Include only the PCs, conditioned on geographic and environmental variables
pRDAstruct <- rda(AllFreq[, -1] ~ PC1 + PC2 + PC3 + Condition(longitude + latitude + BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15), Variables)
# Calculate the adjusted R-squared for the structural model
RsquareAdj(pRDAstruct)

# Output the R-squared and adjusted R-squared for the structural model
$r.squared  # [1] 0.1373694
$adj.r.squared  # [1] 0.1633733

# Perform ANOVA on the structural model to test the significance of the PCs
anova(pRDAstruct)

# Permutation test for the RDA under the reduced model
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = AllFreq[, -1] ~ PC1 + PC2 + PC3 + Condition(longitude + latitude + BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15), data = Variables)
         Df Variance      F Pr(>F)    
Model     3   264.30 6.1317  0.001 ***
Residual 16   229.89                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Pure Geography Model: Only include geographic coordinates, conditioned on environmental and PCA variables
pRDAgeog <- rda(AllFreq[, -1] ~ longitude + latitude + Condition(BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15 + PC1 + PC2 + PC3), Variables)

# Calculate the adjusted R-squared for the geography-only model
RsquareAdj(pRDAgeog)

# Output the R-squared and adjusted R-squared for the geography-only model
$r.squared  # [1] 0.02702213
$adj.r.squared  # [1] 0.01813022

# Perform ANOVA on the geography-only model to test the significance of the geographic factors
anova(pRDAgeog)

# Permutation test for the RDA under the reduced model
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = AllFreq[, -1] ~ longitude + latitude + Condition(BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15 + PC1 + PC2 + PC3), data = Variables)
         Df Variance      F Pr(>F)   
Model     2   51.991 1.8093  0.002 **
Residual 16  229.887                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



########################################################################################################
##########pan east
# Load geographic data (longitude and latitude for 15 populations)
geog = read.csv("15pop_longlat.csv", head = TRUE)
# Load allele frequency data for all populations
AllFreq = fread("pan_east_AF.csv")
# Load neutral loci frequency data
Freq_neutral = fread("F:\\王聪颖\\景观_2\\中性位点Prda\\netural_east_AF.csv", head = TRUE)
# Load environmental variables data
env = read.csv("east_env.csv", head = TRUE)
# Select relevant columns (environmental variables) and scale them (standardization: mean = 0, sd = 1)
ENV = as.data.frame(env[, c(2:7)])
ENV <- scale(ENV, center = TRUE, scale = TRUE)
write.csv(ENV, file = "env.csv", quote = FALSE)


library(vegan)
### Perform PCA on the neutral genetic loci
# RDA is used here to perform PCA (Principal Component Analysis) for scaling the neutral loci data
pca <- rda(Freq_neutral[, -1], scale = TRUE)
# Create a scree plot to visualize eigenvalues from the PCA
screeplot(pca, type = "barplot", npcs = 10, main = "PCA Eigenvalues")
### Extract the first three principal components from PCA
PCs <- scores(pca, choices = c(1:3), display = "sites", scaling = 0)
# Create a data frame with population names and their corresponding first three PCs
PopStruct <- data.frame(POP = Freq_neutral[, 1], PCs)
# Rename columns for better clarity (POP, PC1, PC2, PC3)
colnames(PopStruct) <- c("POP", "PC1", "PC2", "PC3")
write.csv(PopStruct, file = "PopStruct.csv", quote = FALSE)
### Combine all variables into one data frame (geographic, PCA, and environmental data)
Variables <- data.frame(geog, PopStruct[, -1], ENV)
write.csv(Variables, file = "Variables.csv", quote = FALSE)

### Variance Partitioning: Decompose the factors influencing genetic variation
# Perform redundancy analysis (RDA) to analyze the impact of environmental variables, geographic coordinates, and PCs on allele frequencies
pRDAfull <- rda(AllFreq[, -1] ~ PC1 + PC2 + PC3 + longitude + latitude + BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15, Variables)

# Calculate the adjusted R-squared value for the full model
RsquareAdj(pRDAfull)

# Output the R-squared and adjusted R-squared for the full model
$r.squared  # [1] 0.9327759
$adj.r.squared  # [1] 0.6862876

# Perform ANOVA for the full model to test the significance of each factor
anova(pRDAfull)

# Permutation test for the RDA under the reduced model (999 permutations)
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = AllFreq[, -1] ~ PC1 + PC2 + PC3 + longitude + latitude + BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15, data = Variables)
         Df Variance      F Pr(>F)    
Model    11    852.0 3.7843  0.001 ***
Residual  3     61.4                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Climatic Model: Include only the environmental variables, conditioned on geographic and PCA variables
pRDAclim <- rda(AllFreq[, -1] ~ BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15 + Condition(longitude + latitude + PC1 + PC2 + PC3), Variables)

# Calculate the adjusted R-squared for the climatic model
RsquareAdj(pRDAclim)

# Output the R-squared and adjusted R-squared for the climatic model
$r.squared  # [1] 0.2001391
$adj.r.squared  # [1] 0.1021859

# Perform ANOVA on the climatic model to test the significance of the climatic variables
anova(pRDAclim)

# Permutation test for the RDA under the reduced model
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = AllFreq[, -1] ~ BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15 + Condition(longitude + latitude + PC1 + PC2 + PC3), data = Variables)
         Df Variance      F Pr(>F)
Model     6  182.807 1.4886   0.11
Residual  3   61.402

# Structural Model: Include only the PCs, conditioned on geographic and environmental variables
pRDAstruct <- rda(AllFreq[, -1] ~ PC1 + PC2 + PC3 + Condition(longitude + latitude + BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15), Variables)

# Calculate the adjusted R-squared for the structural model
RsquareAdj(pRDAstruct)

# Output the R-squared and adjusted R-squared for the structural model
$r.squared  # [1] 0.1300918
$adj.r.squared  # [1] 0.1466913

# Perform ANOVA on the structural model to test the significance of the PCs
anova(pRDAstruct)

# Permutation test for the RDA under the reduced model
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = AllFreq[, -1] ~ PC1 + PC2 + PC3 + Condition(longitude + latitude + BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15), data = Variables)
         Df Variance      F Pr(>F)
Model     3  118.826 1.9352    0.1
Residual  3   61.402 

## Pure Geography Model: Only include geographic coordinates, conditioned on environmental and PCA variables
pRDAgeog <- rda(AllFreq[, -1] ~ longitude + latitude + Condition(BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15 + PC1 + PC2 + PC3), Variables)

# Calculate the adjusted R-squared for the geography-only model
RsquareAdj(pRDAgeog)

# Output the R-squared and adjusted R-squared for the geography-only model
$r.squared  # [1] 0.04353836
$adj.r.squared  # [1] -0.003577554

# Perform ANOVA on the geography-only model to test the significance of the geographic factors
anova(pRDAgeog)

# Permutation test for the RDA under the reduced model
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = AllFreq[, -1] ~ longitude + latitude + Condition(BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15 + PC1 + PC2 + PC3), data = Variables)
         Df Variance      F Pr(>F)
Model     2   39.768 0.9715  0.528
Residual  3   61.402   
 



######################################################################################################################################
###################core
# Load allele frequency data for core loci
AllFreq <- fread("core_loci2_AF.csv")
### Variance Partitioning: Decompose the factors driving genetic variation
# Perform redundancy analysis (RDA) to analyze the impact of principal components (PC1, PC2, PC3), geographic coordinates (longitude, latitude), and environmental variables (BIO02, BIO03, BIO04, BIO10, BIO12, BIO15) on allele frequencies.
pRDAfull <- rda(AllFreq[,-1] ~ PC1 + PC2 + PC3 + longitude + latitude + BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15, Variables)
# Calculate the adjusted R-squared for the full model (including all factors)
RsquareAdj(pRDAfull)

# Output the R-squared and adjusted R-squared for the full model
$r.squared  # [1] 0.7489621
$adj.r.squared  # [1] 0.6598842

# Perform ANOVA for the full model to test the significance of each factor
anova(pRDAfull)

# Permutation test for the RDA under the reduced model (999 permutations)
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = AllFreq[, -1] ~ PC1 + PC2 + PC3 + longitude + latitude + BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15, data = Variables)
         Df Variance      F Pr(>F)    
Model    11   92.786 8.4079  0.001 ***
Residual 31   31.100                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Climatic Model: Include only the environmental variables, conditioned on geographic and PCA factors
pRDAclim <- rda(AllFreq[,-1] ~ BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15 + Condition(longitude + latitude + PC1 + PC2 + PC3), Variables)

# Calculate the adjusted R-squared for the climatic model
RsquareAdj(pRDAclim)

# Output the R-squared and adjusted R-squared for the climatic model
$r.squared  # [1] 0.1949505
$adj.r.squared  # [1] 0.1661412

# Perform ANOVA on the climatic model to test the significance of the climatic variables
anova(pRDAclim)

# Permutation test for the RDA under the reduced model
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = AllFreq[, -1] ~ BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15 + Condition(longitude + latitude + PC1 + PC2 + PC3), data = Variables)
         Df Variance      F Pr(>F)   
Model     6   24.152 4.0123  0.003 **
Residual 31   31.100                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Structural Model: Include only the PCs, conditioned on geographic and environmental variables
pRDAstruct <- rda(AllFreq[,-1] ~ PC1 + PC2 + PC3 + Condition(longitude + latitude + BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15), Variables)

# Calculate the adjusted R-squared for the structural model
RsquareAdj(pRDAstruct)

# Output the R-squared and adjusted R-squared for the structural model
$r.squared  # [1] 0.06773618
$adj.r.squared  # [1] 0.05366389

# Perform ANOVA on the structural model to test the significance of the PCs
anova(pRDAstruct)

# Permutation test for the RDA under the reduced model
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = AllFreq[, -1] ~ PC1 + PC2 + PC3 + Condition(longitude + latitude + BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15), data = Variables)
         Df Variance      F Pr(>F)   
Model     3   8.3916 2.7882  0.008 **
Residual 31  31.1002                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Pure Geography Model: Only include geographic coordinates, conditioned on environmental and PCA variables
pRDAgeog <- rda(AllFreq[,-1] ~ longitude + latitude + Condition(BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15 + PC1 + PC2 + PC3), Variables)

# Calculate the adjusted R-squared for the geography-only model
RsquareAdj(pRDAgeog)

# Output the R-squared and adjusted R-squared for the geography-only model
$r.squared  # [1] 0.01415166
$adj.r.squared  # [1] -0.002601879

# Perform ANOVA on the geography-only model to test the significance of the geographic factors
anova(pRDAgeog)

# Permutation test for the RDA under the reduced model
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = AllFreq[, -1] ~ longitude + latitude + Condition(BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15 + PC1 + PC2 + PC3), data = Variables)
         Df Variance      F Pr(>F)
Model     2   1.7532 0.8738  0.478
Residual 31  31.1002





#####################################
########core_west
# Load geographic data (longitude and latitude)
geog = read.csv("28pop_longlat.csv", head=T)
# Load allele frequency data for core western loci
AllFreq = fread("core_west_AF.csv") 
# Load allele frequency data for neutral loci (for PCA)
Freq_neutral = fread("F:\\王聪颖\\景观_2\\中性位点Prda\\netural_west_AF.csv", head=T)
# Load environmental data
env = read.csv("Env_west.csv", head=T)
ENV = as.data.frame(env[, c(2:7)])
# Standardize environmental variables (center and scale)
ENV <- scale(ENV, center = TRUE, scale = TRUE)
write.csv(ENV, file="env.csv", quote=F)

library(vegan)
### Perform PCA on the neutral genetic loci to examine variation
# Conduct PCA using the neutral genetic data (excluding the first column, which is assumed to be population IDs)
pca <- rda(Freq_neutral[,-1], scale=T)
# Plot the PCA eigenvalues as a barplot to visualize the explained variance by each PC
screeplot(pca, type = "barplot", npcs = 10, main = "PCA Eigenvalues")
### Output the first three principal components (PCs) for further analysis
PCs <- scores(pca, choices = c(1:3), display = "sites", scaling = 0)
# Combine the population information with the PCA scores
PopStruct <- data.frame(POP = Freq_neutral[, 1], PCs)
colnames(PopStruct) <- c("POP", "PC1", "PC2", "PC3")
write.csv(PopStruct, file = "PopStruct.csv", quote = F)
### Combine all variables (geographic, PCA, environmental) into a single data frame for analysis
Variables <- data.frame(geog, PopStruct[,-1], ENV)
write.csv(Variables, file = "Variables.csv", quote = F)
### Variance Partitioning: Decompose the factors driving genetic variation
# Perform redundancy analysis (RDA) to assess the relative contribution of PCA (PC1, PC2, PC3), geographic factors (longitude, latitude), and environmental variables (BIO02, BIO03, BIO04, BIO10, BIO12, BIO15) on allele frequency data
pRDAfull <- rda(AllFreq[,-1] ~ PC1 + PC2 + PC3 + longitude + latitude + BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15, Variables)

# Calculate the adjusted R-squared value for the full model
RsquareAdj(pRDAfull)

# Output the R-squared and adjusted R-squared values for the full model
$r.squared  # [1] 0.9115053
$adj.r.squared  # [1] 0.8506652

# Perform ANOVA to assess the significance of the model
anova(pRDAfull)

# Permutation test for the RDA under the reduced model (999 permutations)
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = AllFreq[, -1] ~ PC1 + PC2 + PC3 + longitude + latitude + BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15, data = Variables)
         Df Variance      F Pr(>F)    
Model    11  147.079 14.982  0.001 ***
Residual 16   14.279                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Climatic Model: RDA with only environmental variables, conditioned on geographic and PCA variables
pRDAclim <- rda(AllFreq[,-1] ~ BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15 + Condition(longitude + latitude + PC1 + PC2 + PC3), Variables)

# Calculate the adjusted R-squared for the climatic model
RsquareAdj(pRDAclim)

# Output the R-squared and adjusted R-squared for the climatic model
$r.squared  # [1] 0.04762283
$adj.r.squared  # [1] 0.01771853

# Perform ANOVA to assess the significance of the climatic model
anova(pRDAclim)

# Permutation test for the RDA under the reduced model (999 permutations)
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = AllFreq[, -1] ~ BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15 + Condition(longitude + latitude + PC1 + PC2 + PC3), data = Variables)
         Df Variance     F Pr(>F)
Model     6   7.6843 1.435   0.15
Residual 16  14.2793  

# Structural Model: RDA with only PCA components, conditioned on geographic and environmental factors
pRDAstruct <- rda(AllFreq[,-1] ~ PC1 + PC2 + PC3 + Condition(longitude + latitude + BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15), Variables)

# Calculate the adjusted R-squared for the structural model
RsquareAdj(pRDAstruct)

# Output the R-squared and adjusted R-squared for the structural model
$r.squared  # [1] 0.1500628
$adj.r.squared  # [1] 0.189668

# Perform ANOVA to assess the significance of the structural model
anova(pRDAstruct)

# Permutation test for the RDA under the reduced model (999 permutations)
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = AllFreq[, -1] ~ PC1 + PC2 + PC3 + Condition(longitude + latitude + BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15), data = Variables)
         Df Variance      F Pr(>F)    
Model     3   24.214 9.0439  0.001 ***
Residual 16   14.279                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Pure Geography Model: RDA with only geographic coordinates (longitude and latitude), conditioned on environmental and PCA factors
pRDAgeog <- rda(AllFreq[,-1] ~ longitude + latitude + Condition(BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15 + PC1 + PC2 + PC3), Variables)

# Calculate the adjusted R-squared for the geography-only model
RsquareAdj(pRDAgeog)

# Output the R-squared and adjusted R-squared for the geography-only model
$r.squared  # [1] 0.01751558
$adj.r.squared  # [1] 0.009680624

# Perform ANOVA to assess the significance of the geography-only model
anova(pRDAgeog)

# Permutation test for the RDA under the reduced model (999 permutations)
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = AllFreq[, -1] ~ longitude + latitude + Condition(BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15 + PC1 + PC2 + PC3), data = Variables)
         Df Variance      F Pr(>F)
Model     2   2.8263 1.5834  0.119
Residual 16  14.2793






########################################################################################################
##########core east
# Load geographic data (longitude and latitude)
geog = read.csv("15pop_longlat.csv", head = T)
# Load allele frequency data for core eastern loci
AllFreq = fread("core_east_AF.csv") 
# Load allele frequency data for neutral loci (for PCA)
Freq_neutral = fread("F:\\王聪颖\\景观_2\\中性位点Prda\\netural_east_AF.csv", head = T)
# Load environmental data
env = read.csv("east_env.csv", head = T)
ENV = as.data.frame(env[, c(2:7)])
# Standardize environmental variables (center and scale)
ENV <- scale(ENV, center = TRUE, scale = TRUE)
write.csv(ENV, file = "env.csv", quote = F)


library(vegan)
### Perform PCA on the neutral genetic loci to examine variation
# Conduct PCA using the neutral genetic data (excluding the first column, which is assumed to be population IDs)
pca <- rda(Freq_neutral[,-1], scale = T)
# Plot the PCA eigenvalues as a barplot to visualize the explained variance by each PC
screeplot(pca, type = "barplot", npcs = 10, main = "PCA Eigenvalues")
### Output the first three principal components (PCs) for further analysis
PCs <- scores(pca, choices = c(1:3), display = "sites", scaling = 0)
# Combine the population information with the PCA scores
PopStruct <- data.frame(POP = Freq_neutral[, 1], PCs)
# Rename columns for clarity
colnames(PopStruct) <- c("POP", "PC1", "PC2", "PC3")
write.csv(PopStruct, file = "PopStruct.csv", quote = F)
### Combine all variables (geographic, PCA, environmental) into a single data frame for analysis
Variables <- data.frame(geog, PopStruct[,-1], ENV)
write.csv(Variables, file = "Variables.csv", quote = F)
### Variance Partitioning: Decompose the factors driving genetic variation
# Perform redundancy analysis (RDA) to assess the relative contribution of PCA (PC1, PC2, PC3), geographic factors (longitude, latitude), and environmental variables (BIO02, BIO03, BIO04, BIO10, BIO12, BIO15) on allele frequency data
pRDAfull <- rda(AllFreq[,-1] ~ PC1 + PC2 + PC3 + longitude + latitude + BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15 , Variables)

# Calculate the adjusted R-squared value for the full model
RsquareAdj(pRDAfull)

# Output the R-squared and adjusted R-squared values for the full model
$r.squared  # [1] 0.9350977
$adj.r.squared  # [1] 0.6971226

# Perform ANOVA to assess the significance of the model
anova(pRDAfull)

# Permutation test for the RDA under the reduced model (999 permutations)
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = AllFreq[, -1] ~ PC1 + PC2 + PC3 + longitude + latitude + BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15, data = Variables)
         Df Variance      F Pr(>F)    
Model    11  23.5184 3.9294  0.001 ***
Residual  3   1.6323                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Climatic Model: RDA with only environmental variables, conditioned on geographic and PCA variables
pRDAclim <- rda(AllFreq[,-1] ~ BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15 + Condition(longitude + latitude + PC1 + PC2 + PC3), Variables)

# Calculate the adjusted R-squared for the climatic model
RsquareAdj(pRDAclim)

# Output the R-squared and adjusted R-squared for the climatic model
$r.squared  # [1] 0.1357631
$adj.r.squared  # [1] 0.009268757

# Perform ANOVA to assess the significance of the climatic model
anova(pRDAclim)

# Permutation test for the RDA under the reduced model (999 permutations)
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = AllFreq[, -1] ~ BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15 + Condition(longitude + latitude + PC1 + PC2 + PC3), data = Variables)
         Df Variance      F Pr(>F)
Model     6   3.4145 1.0459  0.449
Residual  3   1.6323  

# Structural Model: RDA with only PCA components, conditioned on geographic and environmental factors
pRDAstruct <- rda(AllFreq[,-1] ~ PC1 + PC2 + PC3 + Condition(longitude + latitude + BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15), Variables)

# Calculate the adjusted R-squared for the structural model
RsquareAdj(pRDAstruct)

# Output the R-squared and adjusted R-squared for the structural model
$r.squared  # [1] 0.07997298
$adj.r.squared  # [1] 0.03516494

# Perform ANOVA to assess the significance of the structural model
anova(pRDAstruct)

# Permutation test for the RDA under the reduced model (999 permutations)
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = AllFreq[, -1] ~ PC1 + PC2 + PC3 + Condition(longitude + latitude + BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15), data = Variables)
         Df Variance      F Pr(>F)
Model     3   2.0114 1.2322  0.306
Residual  3   1.6323 

## Pure Geography Model: RDA with only geographic coordinates (longitude and latitude), conditioned on environmental and PCA factors
pRDAgeog <- rda(AllFreq[,-1] ~ longitude + latitude + Condition(BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15 + PC1 + PC2 + PC3), Variables)

# Calculate the adjusted R-squared for the geography-only model
RsquareAdj(pRDAgeog)

# Output the R-squared and adjusted R-squared for the geography-only model
$r.squared  # [1] 0.04684723
$adj.r.squared  # [1] 0.01002129

# Perform ANOVA to assess the significance of the geography-only model
anova(pRDAgeog)

# Permutation test for the RDA under the reduced model (999 permutations)
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = AllFreq[, -1] ~ longitude + latitude + Condition(BIO02 + BIO03 + BIO04 + BIO10 + BIO12 + BIO15 + PC1 + PC2 + PC3), data = Variables)
         Df Variance      F Pr(>F)
Model     2   1.1782 1.0827  0.447
Residual  3   1.6323
