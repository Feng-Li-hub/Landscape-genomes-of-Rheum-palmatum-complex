### Site Annotation
# Pan Annotation for chr
awk '\\$3 == "mRNA" {print \\$0}' Rp.v1.1.gff > Rp.gene.gff  # Extract all mRNA lines from GFF file and save to a new file
tr -s ' ' '\t' < pan.list > pan.list.tab  # Convert spaces to tabs in pan.list and save the result
bedtools intersect -a pan.list.tab -b Rp.gene.gff -wo | awk '{print $(NF-1)}' | sed 's/;.+//' | sort -u > select_gene.intersect  # Find the intersection between pan.list and mRNA file, extract gene info, remove extra parts, and save unique results
bedtools intersect -a pan.list.tab -b ../Rp.gene.gff -wo > anno_pan.txt  # Annotate pan list with gene information from mRNA
awk -F';' '{print \$1}' select_gene.intersect > extracted_ids.txt  # Extract gene IDs from the intersection results and save to a text file

# Core Annotation
tr -s ',' '\t' < core_loci2.csv > core.list  # Convert commas to tabs in core_loci2.csv and save the result
bedtools intersect -a core.list -b ../Rp.gene.gff -wo | awk '{print $(NF-1)}' | sed 's/;.+//' | sort -u > select_gene.intersect  # Find the intersection between core list and mRNA file, extract gene info, remove extra parts, and save unique results
bedtools intersect -a core.list -b ../Rp.gene.gff -wo > anno_core.txt  # Annotate core list with gene information from mRNA
awk -F';' '{print \$1}' select_gene.intersect > extracted_ids.txt  # Extract gene IDs from the intersection results and save to a text file


samtools faidx /home/wangcy/dh_chr/Rh_pal.chr.fa  # Index the genome sequence file Rh_pal.chr.fa using samtools
bedtools makewindows -g Rh_pal.chr.fa.fai -w 1000000 > region.bed  # Create windows of size 1 million base pairs based on the genome file index
bedtools coverage -a region.bed -b /home/wangcy/data/snpdmc/213/jingguan_2/213snp3_ID.vcf -counts > result.txt  # Calculate coverage of SNP data over the defined regions and save results to result.txt

## SNP-Gene Distribution
## Get chromosome names and length information based on the genome
## Plotting guide: https://www.jianshu.com/p/07ae1fe18071

library(RIdeogram)
snp_density = read.table("result.txt", header = TRUE)  # Read SNP count data with column names
karyotype = read.csv("chr_stat.csv")  # Read chromosome statistics information
color = read.csv("color.csv")  # Read color information

# Create a basic chromosome ideogram and save as SVG
ideogram(karyotype = karyotype, 
         output = 'chr.svg')

# Overlay SNP density on the chromosome ideogram and save as SVG
ideogram(karyotype = karyotype, 
         overlaid = snp_density, 
         output = 'chr_gene.svg')

# Add marker labels to the chromosome ideogram and save as SVG
ideogram(karyotype = karyotype, 
         label = color, 
         label_type = "marker", 
         output = 'chr_marker.svg')

# Overlay SNP density and add marker labels, then save as SVG
ideogram(karyotype = karyotype, 
         overlaid = snp_density, 
         label = color, 
         label_type = "marker", 
         output = 'chr_snp_marker.svg')

# Adjust the legend position and output as SVG
ideogram(karyotype = karyotype, 
         overlaid = snp_density, 
         label = color, 
         label_type = "marker", 
         width = 100, Lx = 80, Ly = 25, 
         output = 'chr_snp_marker.svg')

# Convert the generated SVG file to PNG format
convertSVG("chromosome.svg", device = "png")


### GO Annotation
library(clusterProfiler)
gene_GO <- read.delim('gene_GO.txt', stringsAsFactors = FALSE)  # Read GO annotations
gene_KEGG <- read.delim('gene_KEGG.txt', stringsAsFactors = FALSE)  # Read KEGG annotations

# Read gene list with gene names
genes <- read.csv('gene_snp.csv', stringsAsFactors = FALSE)$gene_id

# Enrich KEGG pathways with the gene list
kegg_rich <- enricher(gene = genes,  # Genes to be enriched
                       TERM2GENE = gene_KEGG[c('Pathway', 'gene_id')],  # Background gene set
                       TERM2NAME = gene_KEGG[c('Pathway', 'Description')],
                       pAdjustMethod = 'BH',  # Adjust p-values using BH method
                       pvalueCutoff = 1,  # p-value threshold (set to 1 to include all)
                       qvalueCutoff = 1)  # q-value threshold (set to 1 to include all)

write.table(kegg_rich, 'kegg_rich.txt', sep = '\t', row.names = FALSE, quote = FALSE)  # Save KEGG enrichment results

# Enrich GO terms with the gene list
go_rich <- enricher(gene = genes,  # Genes to be enriched
                    TERM2GENE = gene_GO[c('ID', 'gene_id')],  # Background gene set
                    TERM2NAME = gene_GO[c('ID', 'Description')],
                    pAdjustMethod = 'BH',  # Adjust p-values using BH method
                    pvalueCutoff = 1,  # p-value threshold (set to 1 to include all)
                    qvalueCutoff = 1)  # q-value threshold (set to 1 to include all)

write.table(go_rich, 'go_rich.txt', sep = '\t', row.names = FALSE, quote = FALSE)  # Save GO enrichment results

# Plot GO enrichment as barplot and dotplot
barplot(go_rich)  # Bar plot for GO enrichment
dotplot(go_rich)  # Dot plot for GO enrichment


### File Processing
library(tidyverse)
data <- read.csv("go.csv")

# Split the 'geneID' column by '/' into multiple rows
processed_data <- data %>%

# Summarize data: keep first ID and pvalue, merge descriptions
processed_data2 <- processed_data %>%
  group_by(geneID) %>%
  summarise(
    ID = first(ID),  # Keep first ID for each geneID
    pvalue = first(pvalue),  # Keep first pvalue for each geneID
    Description = paste(Description, collapse = " / "),  # Merge descriptions
    .groups = 'drop'  # Drop the grouping structure
  )
write.table(processed_data2, 'processed_go_summarized.csv', sep = '\t', row.names = FALSE, quote = FALSE)

# View the summarized data
View(processed_data2)

# Save the data as CSV
write.csv(processed_data2, 'processed_go2.csv')

# Split 'geneID' by '/' and summarize
processed_data <- data %>%
  separate_rows(geneID, sep = "/") %>%
  group_by(geneID) %>%
  summarise(
    ONTOLOGY = paste(ONTOLOGY, collapse = "/"),  # Merge ontology terms
    ID = paste(ID, collapse = "/"),  # Merge IDs
    Description = paste(Description, collapse = "/"),  # Merge descriptions
    pvalue = first(pvalue)  # Keep first pvalue for each geneID
  ) %>%
  ungroup()
write.csv(processed_data, "go.csv")

  