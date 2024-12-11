！！！！！！！！！！！！！！！！！！！！！！！！ Structure Pre-analysis ！！！！！！！！！！！！！！！！！！！！！！！！！

# Note: Do not use conda to install Structure v2.3.4 because mainparams and extraparams are not found
####### Software Setup
# Download structure v2.3.4 from the website (Download package without front end):
wget https://web.stanford.edu/group/pritchardlab/structure_software/release_versions/v2.3.4/release/structure_linux_console.tar.gz
tar zxvf structure_linux_console.tar.gz 
cd console
# In the console directory, the following files are generated: extraparams, mainparams, structure, structure_doc.pdf, testdata1

###### 1. Modify parameters in mainparams
#define BURNIN  change to 10000
#define NUMREPS   change to 30000
# Change define INFILE   infile       to  define INFILE   /home/wangcy/data/snpdmc/213/201snpfinal.recode.strct_in     ## For the West: /home/wangcy/data/snpdmc/213/strcuture_311084/west.str
# Change define OUTFILE  outfile  to  define OUTFILE   /home/wangcy/data/snpdmc/213/strcuture/201_1_1  ## For the West: /home/wangcy/data/snpdmc/213/strcuture_311084/west201_1_1
#define NUMINDS    201    // (int) number of diploid individuals in data file ## For the West: 135
#define NUMLOCI    311084    // (int) number of loci in data file  
#define POPDATA   0     // (B) Input file contains a population identifier
#define POPFLAG   0        change to #define POPFLAG   1
#define MARKERNAMES      0  // (B) data file contains row of marker names

###### 2. Modify parameters in extraparams
# Change #define LOCISPOP      1         to  define LOCISPOP      0    //(B) use POPDATA for location information

nohup structure -m ./params/mainparams_1_4 -e ./params/extraparams_4 
nohup structure -m ./params/mainparams_1_5 -e ./params/extraparams_5 
nohup structure -m ./params/mainparams_1_1 -e ./params/extraparams_1 
nohup structure -m ./params/mainparams_1_2 -e ./params/extraparams_2 
nohup structure -m ./params/mainparams_1_3 -e ./params/extraparams_3

K=1-10
### 3. 
python structureHarvester.py --dir=/home/wangcy/data/snpdmc/213/strcuture_311084_1/strcuture --out=/home/wangcy/data/snpdmc/213/strcuture_311084_1/strcuture/out --evanno --clumpp


###### 4. Resampling Analysis for Clusters
####### Software Setup
# Download CLUMPP version 1.1.2 from the website
wget https://web.stanford.edu/group/rosenberglab/software/CLUMPP_Linux64.1.1.2.tar.gz
tar zxvf CLUMPP_Linux64.1.1.2.tar.gz # Unzip and use directly
cd CLUMPP_Linux64.1.1.2/
# In the CLUMPP_Linux64.1.1.2 directory, the following files are generated: arabid.outfile, arabid.permutationfile, CLUMPP, example, arabid.perm_datafile, arabid.popfile, CLUMPP_Manual.pdf, paramfile
###### Modify parameters in the paramfile
# Copy the Kx.indfile obtained from the previous step to the CLUMPP directory and open the paramfile in a text editor (e.g., Notepad++)
# Modify the following parameters:
# DATATYPE: file type to read, 0 means reading *.indfile, 1 means reading *.popfile. I set it to 0
# INDFILE: the name of the individual population file, K2.indfile;
# OUTFILE: the output individual file name, K2_i.outfile;
# K: the number of clusters, K is 2;
# C: total number of individuals, 201;
# R: number of runs, the number of repetitions when Structure runs for a specific K, 5;
###### Run ######
./CLUMPP

###### 5. Graphical Display of Results
####### Software Setup, can run on either server or Windows
# Download distruct v1.1 from the website
wget https://rosenberglab.stanford.edu/software/distruct1.1.tar.gz
tar zxvf distruct1.1.tar.gz      # Unzip and use directly
cd distruct1.1/
# In the distruct1.1 directory, the following files are generated: casia_f, casia.perm, distructLinux1.1, distructWindows1.1.exe, casia.indivq, casia.popq, distructMacOSX1.1, drawparams, casia.languages, casia.ps, distructManual.pdf, casia.names, ColorBrewer, distructWindows1.1.bat

###### Modify files:
# (1) casia.indivq - Individual Q-matrix file, can directly use the previous K2_i.outfile; K2west_i.outfile
# (2) casia.popq - Population Q-matrix file, can directly use the previous K2_poutfile; K2westcasia.popq
# (3) casia.names - Names for the labels under the Structure result chart, e.g., casiak2_2.names; casiawest.names
# (4) casia.languages - Names for the labels above the Structure result chart, e.g., casiak2_2.languages; casiawest.languages
# (5) casia.perm - Colors for different clusters in the Structure result chart
# (6) drawparams - Parameters for distruct plot settings. After defining the above 5 files, also modify K (best K value), NUMPOPS (number of populations), NUMINDS (number of individuals), and output file name (OUTFILE). The output file is a *.ps file, which is a Postscript file

###### Run ######
# 1) On the server
./distructLinux1.1

convert -density 300 k2.ps/k2west.ps k2.png/k2west.png





###### ln(K) ?K
data <- read.csv("data_all.csv", header = TRUE) ## For the West: data <- read.csv("west_group(1).csv", header = TRUE) 
library(ggplot2)
ggplot(data, aes(x = K, y = Delta_K)) +
  scale_x_continuous(breaks = seq(0, max(data$K), 1)) +  # Set the x-axis as integers
  scale_y_continuous(breaks = seq(0, 8000, 2000)) +
  geom_line(size = 1) +
  geom_point(size = 4, shape = 21, fill = "blue") +
  theme_minimal() +  # Use a minimal background theme
  theme(
    panel.grid.major = element_blank(),   # Remove major grid lines
    panel.grid.minor = element_blank(),   # Remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add border around the panel
  )


library(ggplot2)
ggplot(data, aes(x = K, y = MK)) +
  geom_errorbar(aes(ymin = MK - SD, ymax = MK + SD), width = 0.1) +
  geom_line(size = 1) +
  geom_point(size = 4, shape = 21, fill = "blue") +
  scale_x_continuous(breaks = seq(0, max(data$K), 1)) +  # Ensure the x-axis as integers
  theme_minimal() +  # Use a minimal theme
  theme(
    panel.grid.major = element_blank(),   # Remove major grid lines
    panel.grid.minor = element_blank(),   # Remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add a border around the panel
  )






