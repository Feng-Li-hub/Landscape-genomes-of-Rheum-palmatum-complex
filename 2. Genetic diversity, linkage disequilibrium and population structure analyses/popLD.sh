#######
## LD decay can reflect the selection pressure on a population. 
## The r² in PopLDdecay is a key statistic for assessing linkage disequilibrium (LD) between different genetic loci.
## First step is to prepare the VCF file after LD filtering:

plink --bfile 213snp5 --recode vcf --out popla --allow-extra-chr
## Calculate LD decay for the "east" subgroup using the VCF file:
PopLDdecay -InVCF popla.vcf -OutStat LDdecay_east -SubPop east_sample.list -MaxDist 100

## Calculate LD decay for the "west1" subgroup using the VCF file:
PopLDdecay -InVCF popla.vcf -OutStat LDdecay_west1 -SubPop west1_sample.list -MaxDist 100

## Calculate LD decay for the "west2" subgroup using the VCF file:
PopLDdecay -InVCF popla.vcf -OutStat LDdecay_west2 -SubPop west2_sample.list -MaxDist 100 

## To calculate the subgroup LD decay in VCF files, 
## include the group sample names in the respective Group_sample.list files.

## Generate plots for multiple populations using Plot_MultiPop.pl script:
perl /home/wangcy/data/snpdmc/213/diversity/PopLDdecay/PopLDdecay/bin/Plot_MultiPop.pl -inList Pop_k3.list -output k3 -bin1 10 -bin2 5000
## Note: The `-bin1` and `-bin2` options define the binning parameters for plotting.

library(ggplot2)
library(cowplot)
library(data.table)
library(tidyr)

East_quantile = fread("east.stat", skip = 1, header = F)  # Load quantile data from the "east.stat" file, skipping the first row
East_line = read.table("pk_paper.East", header = F)  # Load East line data from the "pk_paper.East" file
East_line <- East_line[-1, ]  # Remove the first row (header)
East_line$V2 <- as.numeric(East_line$V2)  # Convert V2 column to numeric
str(East_line)
East_line$V1 <- as.numeric(East_line$V1)
str(East_line)

# Function to calculate the linkage disequilibrium (LD) quantiles
syp_LD <- function(line_file, quantile_file) {
  line_file = line_file[, c(1, 2)]  # Select only the first two columns (Dist and r2)
  colnames(line_file) = c("Dist", "r2")  # Rename columns for clarity
  
  colnames(quantile_file) = c("Dist", "r2", "count")  # Rename quantile columns
  quantile_file = subset(quantile_file, Dist <= 100000)  # Filter for distances <= 100,000
  line_file = subset(line_file, Dist <= 100000)  # Filter for distances <= 100,000
  
  # Initialize an empty data frame for output
  syp_output <- data.frame(Dist = 0, r2 = 0, quantile_05 = 0, quantile_95 = 0)
  syp_output = syp_output[-1, ]  # Remove the first (empty) row
  
  # Loop to calculate the 5th and 95th percentiles for each distance range
  i = 1
  temp_quantile = subset(quantile_file, Dist <= line_file$Dist[i])  # Filter quantiles for the current distance
  syp_output[i, 1] = line_file[i, 1]  # Store the distance
  syp_output[i, 2] = line_file[i, 2]  # Store the r2 value
  
  # Compute the 5th and 95th percentiles of r2
  value_temp = quantile(rep(temp_quantile$r2, temp_quantile$count), probs = c(0.05, 0.95))
  syp_output[i, 3] = value_temp[1]  # Store 5th percentile
  syp_output[i, 4] = value_temp[2]  # Store 95th percentile
  
  # Iterate through remaining distances
  for (i in 2:length(line_file$Dist)) {
    start = line_file$Dist[i-1]  # Start distance for the current range
    end = line_file$Dist[i]  # End distance for the current range
    temp_quantile = subset(quantile_file, Dist <= end & Dist > start)  # Get quantiles for the current range
    syp_output[i, 1] = line_file[i, 1]  # Store the distance
    syp_output[i, 2] = line_file[i, 2]  # Store the r2 value
    
    # Compute the 5th and 95th percentiles of r2 for this range
    value_temp = quantile(rep(temp_quantile$r2, temp_quantile$count), probs = c(0.05, 0.95))
    syp_output[i, 3] = value_temp[1]  # Store 5th percentile
    syp_output[i, 4] = value_temp[2]  # Store 95th percentile
  }
  
  return(syp_output)  # Return the computed LD decay data
}

# Call the function and store the output
plot_file = syp_LD(East_line, East_quantile)
write.table(plot_file, "plot_file_East.txt", quote = F, row.names = F)


### plot:
East_plot = ggplot() + 
  # Create a shaded ribbon to represent the 5th and 95th percentiles (quantiles)
  geom_ribbon(data = plot_file, aes(ymin = quantile_05, ymax = quantile_95, x = Dist / 1000), 
              fill = "coral1", alpha = 0.2) +  # Set the ribbon color to coral1 with transparency
  # Plot a line representing the r2 values over distance
  geom_line(data = plot_file, aes(x = Dist / 1000, y = r2), col = "coral1", size = 2) +  # Line color is coral1, size is 2
  xlab("Distance(kb)") + ylab(expression(r^2)) +  # X-axis label in kilobases (kb) and y-axis label with r² notation
  theme_bw() +  
  theme(
    text = element_text(family = "serif"),
    axis.text.x = element_text(hjust = 1, size = 13, family = "serif"),  
    axis.text.y = element_text(size = 13, family = "serif"),    
    plot.title = element_text(size = 12L, hjust = 0.5), 
    axis.ticks = element_line(colour = "black", size = .001, linetype = 1, lineend = .1),   
    axis.title.y = element_text(size = 14, family = "serif"),  
    plot.subtitle = element_text(hjust = 0, size = 12, family = "serif"), 
    legend.position = "none",  # Remove the legend
    axis.title.x = element_text(size = 16, family = "serif"), 
    panel.background = element_rect(fill = "white"),  
    plot.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    plot.margin = unit(c(1, 1, 1, 1), "mm"),  
    panel.border = element_blank(),  
    axis.line = element_line(color = "black", size = 0.6)  
  )


	
	
# Read west1 data and convert columns to numeric
west1_quantile = fread("west1.stat", skip = 1, header = F)
west1_line = read.table("pk_paper.west1", header = F)
west1_line$V2 <- as.numeric(west1_line$V2)  # Convert V2 to numeric
west1_line$V1 <- as.numeric(west1_line$V1)  # Convert V1 to numeric
# Check structure of the data
str(west1_line)
# Generate the plot data using syp_LD function
plot_file2 = syp_LD(west1_line, west1_quantile)
# Save the plot data to a file
write.table(plot_file2, "plot_file_west1.txt", quote = F, row.names = F)
# Create the plot for west1
west1_plot = ggplot() +
  geom_ribbon(data = plot_file2, aes(ymin = quantile_05, ymax = quantile_95, x = Dist / 1000), fill = "lightslateblue", alpha = 0.2) +
  geom_line(data = plot_file, aes(x = Dist / 1000, y = r2), col = "lightslateblue", size = 2) +
  xlab("Distance(kb)") + ylab(expression(r^2)) +
  theme_bw() +
  theme(
    text = element_text(family = "serif"),
    axis.text.x = element_text(hjust = 1, size = 13, family = "serif"),
    axis.text.y = element_text(size = 13, family = "serif"),
    plot.title = element_text(size = 12L, hjust = 0.5),
    axis.ticks = element_line(colour = "black", size = .001, linetype = 1, lineend = .1),
    axis.title.y = element_text(size = 14, family = "serif"),
    plot.subtitle = element_text(hjust = 0, size = 12, family = "serif"),
    legend.position = "none",
    axis.title.x = element_text(size = 16, family = "serif"),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "mm"),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.6)
  )

# Read west2 data and convert columns to numeric
west2_quantile = fread("west2.stat", skip = 1, header = F)
west2_line = read.table("pk_paper.west2", header = F)
west2_line$V2 <- as.numeric(west2_line$V2)  # Convert V2 to numeric
west2_line$V1 <- as.numeric(west2_line$V1)  # Convert V1 to numeric
# Check structure of the data
str(west2_line)
# Generate the plot data for west2 using syp_LD function
plot_file3 = syp_LD(west2_line, west2_quantile)
write.table(plot_file3, "plot_file_west2.txt", quote = F, row.names = F)
# Create the plot for west2
west2_plot = ggplot() +
  geom_ribbon(data = plot_file3, aes(ymin = quantile_05, ymax = quantile_95, x = Dist / 1000), fill = "olivedrab3", alpha = 0.2) +
  geom_line(data = plot_file, aes(x = Dist / 1000, y = r2), col = "olivedrab3", size = 2) +
  xlab("Distance(kb)") + ylab(expression(r^2)) +
  theme_bw() +
  theme(
    text = element_text(family = "serif"),
    axis.text.x = element_text(hjust = 1, size = 13, family = "serif"),
    axis.text.y = element_text(size = 13, family = "serif"),
    plot.title = element_text(size = 12L, hjust = 0.5),
    axis.ticks = element_line(colour = "black", size = .001, linetype = 1, lineend = .1),
    axis.title.y = element_text(size = 14, family = "serif"),
    plot.subtitle = element_text(hjust = 0, size = 12, family = "serif"),
    legend.position = "none",
    axis.title.x = element_text(size = 16, family = "serif"),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "mm"),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.6)
  )

	
library(ggplot2)
combined_plot <- ggplot(data = all_data, aes(x = Dist / 1000)) +  # Use 'all_data' dataset, scale distance by dividing by 1000 (convert to kb)
  geom_ribbon(aes(ymin = quantile_05, ymax = quantile_95, fill = region), alpha = 0.2) +  # Draw shaded region with quantile 5-95% for each region
  geom_line(aes(y = r2, color = region), size = 2) +  # Plot line for r2 values, color by region
  xlab("Distance (kb)") +  # Label for x-axis
  ylab(expression(r^2)) +  # Label for y-axis
  theme_bw() +  # Use a white background theme
  theme(
    text = element_text(family = "serif"),  # Set font to serif
    axis.text.x = element_text(hjust = 1, size = 13, family = "serif"),  # Customize x-axis text
    axis.text.y = element_text(size = 13, family = "serif"),  # Customize y-axis text
    plot.title = element_text(size = 12L, hjust = 0.5),  # Customize plot title
    axis.ticks = element_line(colour = "black", size = .001, linetype = 1, lineend = .1),  # Customize axis ticks
    axis.title.y = element_text(size = 14, family = "serif"),  # Customize y-axis title
    plot.subtitle = element_text(hjust = 0, size = 12, family = "serif"),  # Customize subtitle
    legend.position = "none",  # Hide the legend
    axis.title.x = element_text(size = 16, family = "serif"),  # Customize x-axis title
    panel.background = element_rect(fill = "white"),  # Set panel background to white
    plot.background = element_rect(fill = "white"),  # Set plot background to white
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    plot.margin = unit(c(1, 1, 1, 1), "mm"),  # Set plot margins
    panel.border = element_blank(),  # Remove panel border
    axis.line = element_line(color = "black", size = 0.6)  # Set axis line color and size
  ) +
  scale_fill_manual(values = c("East" = "coral1", "West 1" = "lightslateblue", "West 2" = "olivedrab3")) +  # Set fill colors for regions
  scale_color_manual(values = c("East" = "coral1", "West 1" = "lightslateblue", "West 2" = "olivedrab3"))  # Set line colors for regions
print(combined_plot)
