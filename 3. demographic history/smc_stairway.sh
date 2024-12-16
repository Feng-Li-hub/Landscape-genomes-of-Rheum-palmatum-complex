###SMCPP
# **Compress the VCF file**: 
bgzip top11chrAll_LD.vcf  # Compress the VCF file
# **Build an index for the VCF file**:
tabix top11chrAll_LD.vcf.gz  # Generate a tbi index for the compressed VCF file
# **SMCPP estimation command**:
smc++ estimate --spline cubic --knots 15 --timepoints 1000 1000000 --cores 32 -o /home/puyang/DH/DH_allsite/rh_dynamic/estimate/rh_west1 7.0e-09 /data/data01/DH_smcpp/rh_dynamic/west1/*.smc.gz  # Perform population history estimation for the west1 population using SMCPP
smc++ estimate --spline cubic --knots 15 --timepoints 1000 1000000 --cores 32 -o /home/puyang/DH/DH_allsite/rh_dynamic/estimate/rh_west2 7.0e-09 /data/data01/DH_smcpp/rh_dynamic/west2/*.smc.gz  # Perform population history estimation for the west2 population using SMCPP
smc++ estimate --spline cubic --knots 15 --timepoints 1000 1000000 --cores 32 -o /home/puyang/DH/DH_allsite/rh_dynamic/estimate/rh_east 7.0e-09 /data/data01/DH_smcpp/rh_dynamic/east/*.smc.gz  # Perform population history estimation for the east population using SMCPP
# **Plot the estimation results**:
smc++ plot --csv -g 5 --ylim 0 20000 --xlim 0 1000000 --linear rh_west1.pdf /home/puyang/DH/DH_allsite/rh_dynamic/estimate/rh_west1/model.final.json  # Plot the estimation result for the west1 population
smc++ plot --csv -g 5 --ylim 0 20000 --xlim 0 1000000 --linear rh_west2.pdf /home/puyang/DH/DH_allsite/rh_dynamic/estimate/rh_west2/model.final.json  # Plot the estimation result for the west2 population
smc++ plot --csv -g 5 --ylim 0 30000 --xlim 0 1000000 --linear rh_east.pdf /home/puyang/DH/DH_allsite/rh_dynamic/estimate/rh_east/model.final.json  # Plot the estimation result for the east population
smc++ plot --csv -g 5 --ylim 0 20000 --xlim 0 1000000 --linear rh_all.pdf /home/puyang/DH/DH_allsite/rh_dynamic/estimate/rh_all/model.final.json  # Plot the estimation result for all populations combined
smc++ plot --csv -g 5 --ylim 0 30000 --xlim 0 1000000 --linear rh_dynamic.pdf /home/puyang/DH/DH_allsite/rh_dynamic/estimate/rh_west1/model.final.json /home/puyang/DH/DH_allsite/rh_dynamic/estimate/rh_west2/model.final.json /home/puyang/DH/DH_allsite/rh_dynamic/estimate/rh_east/model.final.json  # Plot the combined estimation results for multiple populations

# **Plot in R**:
# Load necessary R libraries
library(dplyr)  # For data manipulation
library(ggplot2)  # For plotting

# Read in CSV data
data <- read.csv("rh_dynamic.csv")  # Read the CSV file containing data
# Custom color palette
col <- c("coral1", "lightslateblue", "olivedrab3")  # Define custom colors for plotting
# Filter data to include only entries where the x value is <= 1 million years
new_data <- data %>%
  filter(x <= 1000000)  # Filter out data where x (years) > 1,000,000
# Plot using ggplot2
ggplot(new_data, aes(x = x, y = y, color = label)) +  # Map x and y to axes and color by label
  geom_line(size = 1) +  # Draw a line with a width of 1
  scale_color_manual(values = col) +  # Use the custom color palette
  labs(x = "Years (in thousands)", y = "Estimated Ne/10^4", color = "Lineage") +  # Set axis labels and legend title
  theme_minimal() +  # Apply a minimalistic theme
  scale_x_continuous(expand = c(0, 0), limits = c(0, 200000), labels = function(x) x / 1000) +  # Set x-axis limits and convert labels to thousands of years
  scale_y_continuous(expand = c(0, 0), limits = c(0, 50000), labels = function(y) y / 10000) +  # Set y-axis limits and convert labels to 10,000s
  annotate("rect", xmin = 0, xmax = 12000, ymin = 80000, ymax = 90000, fill = "#A14462", alpha = 1) +  # Add a rectangle for the Holocene
  annotate("rect", xmin = 12000, xmax = 1000000, ymin = 80000, ymax = 90000, fill = "#5D90BA", alpha = 0.78) +  # Add a rectangle for the Pleistocene
  annotate("text", x = 500000, y = 85500, label = "Pleistocene", size = 5, color = "black", hjust = 0.5) +  # Add label for Pleistocene
  annotate("rect", xmin = 11700, xmax = 126000, ymin = 0, ymax = 80000, fill = "gray", alpha = 0.23) +  # Add a rectangle for the Late Pleistocene
  annotate("rect", xmin = 190000, xmax = 265000, ymin = 0, ymax = 80000, fill = "gray", alpha = 0.23) +  # Add a rectangle for the LGM (Last Glacial Maximum)
  theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill = NA)) +  # Customize legend and grid lines
  coord_cartesian(clip = "on")  # Ensure all plot elements remain visible within the plot area



###stairway
#  **Running EasySFS for SFS generation**:
# Run EasySFS script to generate SFS files
nohup python /home/wangcy/data/snpdmc/213/diversity/stairway/easySFS-master/easySFS-master/easySFS.py -i top11chrAll_LD.vcf -p pops_file.txt --preview -a > easySFS.log 2>&1 &  # Generate SFS preview
nohup python /home/wangcy/data/snpdmc/213/diversity/stairway/easySFS-master/easySFS-master/easySFS.py -i top11chrAll_LD.vcf -p pops_file.txt -a --proj=96,108,74  # Generate SFS projection with specific populations

# Generate Stairway Plot
java -cp stairway_plot_es Stairbuilder rh_east_fold.blueprint  # Generate Stairway Plot for east population
java -cp stairway_plot_es Stairbuilder rh_west1_fold.blueprint  # Generate Stairway Plot for west1 population
java -cp stairway_plot_es Stairbuilder rh_west2_fold.blueprint  # Generate Stairway Plot for west2 population

#  **Running Bash scripts**:
nohup bash east_bash.sh &  # Execute the bash script for east population
nohup bash west1_bash.sh &  # Execute the bash script for west1 population
nohup bash west2_bash.sh &  # Execute the bash script for west2 population


data <- read.csv("stair.csv") 
col <- c("coral1", "lightslateblue", "olivedrab3")  # Define custom colors for plotting
data$group <- as.factor(data$group)  # Convert the 'group' column to a factor for categorical grouping

library(ggplot2) 
ggplot(data, aes(x = log(year, 10), y = Ne_median / 1000, color = group)) + 
    geom_line(size = 1) +  # Plot the median Ne curve with line thickness of 1
    scale_color_manual(values = col, labels = c("E", "NW", "SW")) +  # Apply custom colors and labels for the groups
    scale_y_log10() +  # Set y-axis to log scale

    # Add upper and lower confidence interval lines
    geom_line(aes(y = Ne_2.5. / 1000, color = group), linetype = "dashed", size = 0.8) +  # Lower confidence interval
    geom_line(aes(y = Ne_97.5. / 1000, color = group), linetype = "dashed", size = 0.8) +  # Upper confidence interval

    # Customize x-axis
    scale_x_continuous(
        breaks = c(log(100, 10), log(500, 10), log(1000, 10), log(5000, 10), log(10000, 10), 
                   log(50000, 10), log(100000, 10), log(500000, 10), log(1000000, 10), 
                   log(5000000, 10), log(10000000, 10), log(50000000, 10)),  # Define x-axis breakpoints
        labels = c("0.1", "0.5", "1", "5", "10", "50", "100", "500", "1000", "5000", "10000", "50000"),  # Custom labels
        limits = c(log(100, 10), log(max(data$year), 10) + 0.1)  # Set x-axis limits
    ) + 

    # Add plot labels
    labs(title = "Demography of different lineages", x = "kilo years ago", y = "Ne (1k individual)") +  # Title and axis labels

    # Set plot theme
    theme_minimal() +  # Use a minimal theme for the plot
    theme(axis.title = element_text(size = 16),  # Set axis title font size
          axis.text = element_text(size = 14),  # Set axis label font size
          axis.text.x = element_text(angle = 30, hjust = 1),  # Rotate x-axis labels
          panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank())  # Remove minor grid lines






