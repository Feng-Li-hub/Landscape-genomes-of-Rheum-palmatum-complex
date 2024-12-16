### SNP Genetic Load

perl make-SIFT-db-all.pl -config ./test_files/Rheum_palmatum_config2.txt  # (SIFT4G_Create_Genomic_DB to build the database)

## Initial Directories
chr-src --- Rh_pal.chr.fa.gz
gene-annotation-src --- **.gtf.gz
dbSNP --- No file

nohup java -jar /home/zhangyi/241008/SIFT4G_Annotator/SIFT4G_Annotator.jar -c -i /home/zhangyi/241008/anno/213snp3_ID.vcf -d /home/zhangyi/241008/result/Rp.V1 -r /home/zhangyi/241008/anno -t > sift4g_213snp3ID.log 2>&1 &

####
# Loss of Function Sites
grep "START\|STOP" 213snp3_ID_SIFTannotations.xls > rh.cds.loss_of_function.xls

# Deleterious Sites
awk '\$9=="NONSYNONYMOUS"' 213snp3_ID_SIFTannotations.xls | awk '\$17=="DELETERIOUS"' > rh.cds.deleterious.xls

# Tolerated Sites
awk '\$9=="NONSYNONYMOUS"' 213snp3_ID_SIFTannotations.xls | grep -v "DELETERIOUS" | awk '\$13!="NA"' > rh.cds.tolerated.xls

# Synonymous Sites
awk '\$9=="SYNONYMOUS"' 213snp3_ID_SIFTannotations.xls | grep -v "DELETERIOUS" | awk '\$13!="NA"' > rh.cds.synonymous.xls  ## (Filter using Excel as per command)

# Extract the first and second columns from the input files and replace tabs with colons, save as SNP files
cut -f 1,2 rh.cds.loss_of_function.xls | sed 's/\t/:/g' > rh.cds.loss_of_function.snps
cut -f 1,2 rh.cds.deleterious.xls | sed 's/\t/:/g' > rh.cds.deleterious.snps
cut -f 1,2 rh.cds.tolerated.xls | sed 's/\t/:/g' > rh.cds.tolerated.snps
cut -f 1,2 rh.cds.synonymous.xls | sed 's/\t/:/g' > rh.cds.synonymous.snps

# Extract the first two columns from 213snp3_ID_SIFTannotations.xls and replace tabs with colons
cut -f 1,2 213snp3_ID_SIFTannotations.xls | sed 's/\t/:/g' > 213snp3_SIFTannotations.all.snps

awk 'BEGIN {OFS="\t"} !/^#/ { \$3 = \$1":"\$2 } { print }' 213snp3_ID.vcf > 213snp3_ID2.vcf


# Extract records related to specific SNPs from VCF file
vcftools --vcf 213snp3_ID2.vcf --snps rh.cds.deleterious.snps --recode --recode-INFO-all --out rh.snp.deleterious
vcftools --vcf 213snp3_ID2.vcf --snps rh.cds.loss_of_function.snps --recode --recode-INFO-all --out rh.snp.loss_of_function
vcftools --vcf 213snp3_ID2.vcf --snps rh.cds.synonymous.snps --recode --recode-INFO-all --out rh.snp.synonymous
vcftools --vcf 213snp3_ID2.vcf --snps rh.cds.tolerated.snps --recode --recode-INFO-all --out rh.snp.tolerated

# Compress VCF files and remove original files
bcftools view rh.snp.deleterious.recode.vcf -O z -o rh.snp.deleterious.recode.vcf.gz && rm rh.snp.deleterious.recode.vcf
bcftools view rh.snp.loss_of_function.recode.vcf -O z -o rh.snp.loss_of_function.recode.vcf.gz && rm rh.snp.loss_of_function.recode.vcf
bcftools view rh.snp.synonymous.recode.vcf -O z -o rh.snp.synonymous.recode.vcf.gz && rm rh.snp.synonymous.recode.vcf
bcftools view rh.snp.tolerated.recode.vcf -O z -o rh.snp.tolerated.recode.vcf.gz && rm rh.snp.tolerated.recode.vcf

# Remove FORMAT information from VCF files and compress them
bcftools annotate -x FORMAT -O z -o rh.snp.deleterious.simple.vcf.gz rh.snp.deleterious.recode.vcf.gz && rm rh.snp.deleterious.recode.vcf.gz
bcftools annotate -x FORMAT -O z -o rh.snp.loss_of_function.simple.vcf.gz rh.snp.loss_of_function.recode.vcf.gz && rm rh.snp.loss_of_function.recode.vcf.gz
bcftools annotate -x FORMAT -O z -o rh.snp.synonymous.simple.vcf.gz rh.snp.synonymous.recode.vcf.gz && rm rh.snp.synonymous.recode.vcf.gz
bcftools annotate -x FORMAT -O z -o rh.snp.tolerated.simple.vcf.gz rh.snp.tolerated.recode.vcf.gz && rm rh.snp.tolerated.recode.vcf.gz

# Extract records related to intersected SNPs
vcftools --gzvcf rh.snp.deleterious.simple.vcf.gz --snps deleterious.intersect2.snp.txt --recode --recode-INFO-all --out rh.deleterious.intersect
vcftools --gzvcf rh.snp.loss_of_function.simple.vcf.gz --snps loss_of_function.intersect2.snp.txt --recode --recode-INFO-all --out rh.loss_of_function.intersect
vcftools --gzvcf rh.snp.synonymous.simple.vcf.gz --snps synonymous.intersect2.snp.txt --recode --recode-INFO-all --out rh.synonymous.intersect
vcftools --gzvcf rh.snp.tolerated.simple.vcf.gz --snps tolerated.intersect2.snp.txt --recode --recode-INFO-all --out rh.tolerated.intersect

# Extract intersected SNPs from another VCF file
vcftools --gzvcf Rp_OdRn2_213.result2.vcf.gz --snps deleterious.intersect2.snp.txt --recode --recode-INFO-all --out tahuang.deleterious.intersect2
vcftools --gzvcf Rp_OdRn2_213.result2.vcf.gz --snps loss_of_function.intersect2.snp.txt --recode --recode-INFO-all --out tahuang.loss_of_function.intersect2
vcftools --gzvcf Rp_OdRn2_213.result2.vcf.gz --snps synonymous.intersect2.snp.txt --recode --recode-INFO-all --out tahuang.synonymous.intersect2
vcftools --gzvcf Rp_OdRn2_213.result2.vcf.gz --snps tolerated.intersect2.snp.txt --recode --recode-INFO-all --out tahuang.tolerated.intersect2

# Compress and remove original intersect results
bcftools view rh.deleterious.intersect.recode.vcf -O z -o rh.deleterious.intersect.recode.vcf.gz && rm rh.deleterious.intersect.recode.vcf
bcftools view rh.loss_of_function.intersect.recode.vcf -O z -o rh.loss_of_function.intersect.recode.vcf.gz && rm rh.loss_of_function.intersect.recode.vcf
bcftools view rh.synonymous.intersect.recode.vcf -O z -o rh.synonymous.intersect.recode.vcf.gz && rm rh.synonymous.intersect.recode.vcf
bcftools view rh.tolerated.intersect.recode.vcf -O z -o rh.tolerated.intersect.recode.vcf.gz && rm rh.tolerated.intersect.recode.vcf

# Index compressed VCF files
bcftools index rh.deleterious.intersect.recode.vcf.gz
bcftools index rh.loss_of_function.intersect.recode.vcf.gz
bcftools index rh.synonymous.intersect.recode.vcf.gz
bcftools index rh.tolerated.intersect.recode.vcf.gz


# Compress and delete original intersect results
bcftools view tahuang.deleterious.intersect2.recode.vcf -O z -o tahuang.deleterious.intersect2.recode.vcf.gz && rm tahuang.deleterious.intersect2.recode.vcf
bcftools view tahuang.loss_of_function.intersect2.recode.vcf -O z -o tahuang.loss_of_function.intersect2.recode.vcf.gz && rm tahuang.loss_of_function.intersect2.recode.vcf
bcftools view tahuang.synonymous.intersect2.recode.vcf -O z -o tahuang.synonymous.intersect2.recode.vcf.gz && rm tahuang.synonymous.intersect2.recode.vcf
bcftools view tahuang.tolerated.intersect2.recode.vcf -O z -o tahuang.tolerated.intersect2.recode.vcf.gz && rm tahuang.tolerated.intersect2.recode.vcf

# Index compressed VCF files
bcftools index tahuang.deleterious.intersect2.recode.vcf.gz
bcftools index tahuang.loss_of_function.intersect2.recode.vcf.gz
bcftools index tahuang.synonymous.intersect2.recode.vcf.gz
bcftools index tahuang.tolerated.intersect2.recode.vcf.gz



# Extract header and modify VCF file header
bcftools view -h tahuang.deleterious.intersect2.recode.vcf.gz > header.txt
##FILTER=<ID=SnpCluster,Description="Variant in a SNP cluster">
##FILTER=<ID=lowQualFilter,Description="Low quality variant">
bcftools reheader -h header.txt tahuang.deleterious.intersect2.recode.vcf.gz -o tahuang.deleterious.intersect3.recode.vcf.gz
bcftools index tahuang.deleterious.intersect3.recode.vcf.gz

bcftools view -h tahuang.loss_of_function.intersect2.recode.vcf.gz > header.txt
bcftools reheader -h header.txt tahuang.loss_of_function.intersect2.recode.vcf.gz -o tahuang.loss_of_function.intersect3.recode.vcf.gz
bcftools index tahuang.loss_of_function.intersect3.recode.vcf.gz

bcftools view -h tahuang.synonymous.intersect2.recode.vcf.gz > header.txt
bcftools reheader -h header.txt tahuang.synonymous.intersect2.recode.vcf.gz -o tahuang.synonymous.intersect3.recode.vcf.gz
bcftools index tahuang.synonymous.intersect3.recode.vcf.gz

bcftools view -h tahuang.tolerated.intersect2.recode.vcf.gz > header.txt
bcftools reheader -h header.txt tahuang.tolerated.intersect2.recode.vcf.gz -o tahuang.tolerated.intersect3.recode.vcf.gz
bcftools index tahuang.tolerated.intersect3.recode.vcf.gz

# Merge SNP results with outgroup SNP results
bcftools merge -m snps rh.deleterious.intersect.recode.vcf.gz tahuang.deleterious.intersect3.recode.vcf.gz -O z -o rh.deleterious.outgroup.vcf.gz
bcftools merge -m snps rh.loss_of_function.intersect.recode.vcf.gz tahuang.loss_of_function.intersect3.recode.vcf.gz -O z -o rh.loss_of_function.outgroup.vcf.gz
bcftools merge -m snps rh.synonymous.intersect.recode.vcf.gz tahuang.synonymous.intersect3.recode.vcf.gz -O z -o rh.synonymous.outgroup.vcf.gz
bcftools merge -m snps rh.tolerated.intersect.recode.vcf.gz tahuang.tolerated.intersect3.recode.vcf.gz -O z -o rh.tolerated.outgroup.vcf.gz

# Extract frequency information for specific species
vcftools --gzvcf ../rh.deleterious.outgroup.vcf.gz --keep pop.tsv --freq --out rh.deleterious.est_sfs
vcftools --gzvcf ../rh.loss_of_function.outgroup.vcf.gz --keep pop.tsv --freq --out rh.loss_of_function.est_sfs
vcftools --gzvcf ../rh.synonymous.outgroup.vcf.gz --keep pop.tsv --freq --out rh.synonymous.est_sfs
vcftools --gzvcf ../rh.tolerated.outgroup.vcf.gz --keep pop.tsv --freq --out rh.tolerated.est_sfs

# Extract frequency information for outgroup
vcftools --gzvcf ../tahuang.tolerated.intersect3.recode.vcf.gz --freq --out tahuang.tolerated.outgroup1.est_sfs
vcftools --gzvcf ../tahuang.loss_of_function.intersect3.recode.vcf.gz --freq --out tahuang.loss_of_function.outgroup1.est_sfs
vcftools --gzvcf ../tahuang.deleterious.intersect3.recode.vcf.gz --freq --out tahuang.deleterious.outgroup1.est_sfs
vcftools --gzvcf ../tahuang.synonymous.intersect3.recode.vcf.gz --freq --out tahuang.synonymous.outgroup1.est_sfs



## est-sfs
/home/zhangyi/241008/anno/est-sfs/est-sfs-release-2.04
est_sfs_input.rh.pl (just modify the file names in the script)

cut -d $'\t' -f 2- rh_synonymous_input.txt > rh_synonymous_input2.txt
cut -d $'\t' -f 2- rh_tolerated_input.txt > rh_tolerated_input2.txt
cut -d $'\t' -f 2- rh_deleterious_input.txt > rh_deleterious_input2.txt
cut -d $'\t' -f 2- rh_loss_of_function_input.txt > rh_loss_of_function_input2.txt

gen <- read.csv("rh_synonymous_input.csv") # If the first column data lacks a header, use row.names=1, otherwise do not use this function to read allele frequency data from all populations
dim(gen) # Check the file's integrity
sum(is.na(gen)) # Check how many genotypes are missing
gen.imp <- apply(gen, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x)))))) # Replace missing data with the most common genotype for each site
sum(is.na(gen.imp)) # Check if there is still missing data
write.csv(gen.imp, "rh_synonymous_input_no_NA.csv") # Save data with no missing values

gen <- read.csv("rh_deleterious_input2.csv") # If the first column data lacks a header, use row.names=1, otherwise do not use this function to read allele frequency data from all populations
dim(gen) # Check the file's integrity
sum(is.na(gen)) # Check how many genotypes are missing
gen.imp <- apply(gen, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x)))))) # Replace missing data with the most common genotype for each site
sum(is.na(gen.imp)) # Check if there is still missing data
write.csv(gen.imp, "rh_deleterious_input2_no_NA.csv") # Save data with no missing values

gen <- read.csv("rh_loss_of_function_input2.csv") # If the first column data lacks a header, use row.names=1, otherwise do not use this function to read allele frequency data from all populations
dim(gen) # Check the file's integrity
sum(is.na(gen)) # Check how many genotypes are missing
gen.imp <- apply(gen, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x)))))) # Replace missing data with the most common genotype for each site
sum(is.na(gen.imp)) # Check if there is still missing data
write.csv(gen.imp, "rh_loss_of_function_input2_no_NA.csv") # Save data with no missing values

gen <- read.csv("rh_tolerated_input2.csv") # If the first column data lacks a header, use row.names=1, otherwise do not use this function to read allele frequency data from all populations
dim(gen) # Check the file's integrity
sum(is.na(gen)) # Check how many genotypes are missing
gen.imp <- apply(gen, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x)))))) # Replace missing data with the most common genotype for each site
sum(is.na(gen.imp)) # Check if there is still missing data
write.csv(gen.imp, "rh_tolerated_input2_no_NA.csv") # Save data with no missing values

awk '{print \$1","\$2","\$3","\$4 "\t" \$5","\$6","\$7","\$8}' rh_synonymous_input_no_NA.txt > rh_synonymous_input_no_NA2.txt 
awk '{print \$1","\$2","\$3","\$4 "\t" \$5","\$6","\$7","\$8}' rh_synonymous_input_no_NA.txt > rh_synonymous_input_no_NA2.txt 
awk '{print \$1","\$2","\$3","\$4 "\t" \$5","\$6","\$7","\$8}' rh_synonymous_input_no_NA.txt > rh_synonymous_input_no_NA2.txt 
awk '{print \$1","\$2","\$3","\$4 "\t" \$5","\$6","\$7","\$8}' rh_synonymous_input_no_NA.txt > rh_synonymous_input_no_NA2.txt 

#rate6
/home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/est-sfs /home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/config-rate6.txt rh_deleterious_input2.txt /home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/seedfile.txt rh.deleterious.est_sfs.output.txt rh.deleterious.est_sfs.output.p_anc.txt
/home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/est-sfs /home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/config-rate6.txt rh_loss_of_function_input2.txt /home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/seedfile.txt rh.loss_of_function.est_sfs.output.txt rh.loss_of_function.est_sfs.output.p_anc.txt
/home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/est-sfs /home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/config-rate6.txt rh_synonymous_input2.txt /home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/seedfile.txt rh.synonymous.est_sfs.output.txt rh.synonymous.est_sfs.output.p_anc.txt
/home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/est-sfs /home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/config-rate6.txt rh_tolerated_input2.txt /home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/seedfile.txt rh.tolerated.est_sfs.output.txt rh.tolerated.est_sfs.output.p_anc.txt



#jc
/home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/est-sfs /home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/config-JC.txt rh_deleterious_input2.txt /home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/seedfile.txt rh.deleterious.est_sfs.jc.output.txt rh.deleterious.est_sfs.jc.output.p_anc.txt
/home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/est-sfs /home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/config-JC.txt rh_loss_of_function_input2.txt /home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/seedfile.txt rh.loss_of_function.est_sfs.jc.output.txt rh.loss_of_function.est_sfs.jc.output.p_anc.txt
/home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/est-sfs /home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/config-JC.txt rh_synonymous_input2.txt /home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/seedfile.txt rh.synonymous.est_sfs.jc.output.txt rh.synonymous.est_sfs.jc.output.p_anc.txt
/home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/est-sfs /home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/config-JC.txt rh_tolerated_input2.txt /home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/seedfile.txt rh.tolerated.est_sfs.jc.output.txt rh.tolerated.est_sfs.jc.output.p_anc.txt


#kimura
/home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/est-sfs /home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/config-kimura.txt rh_deleterious_input2.txt /home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/seedfile.txt rh.deleterious.est_sfs.kimura.output.txt rh.deleterious.est_sfs.kimura.output.p_anc.txt
/home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/est-sfs /home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/config-kimura.txt rh_loss_of_function_input2.txt /home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/seedfile.txt rh.loss_of_function.est_sfs.kimura.output.txt rh.loss_of_function.est_sfs.kimura.output.p_anc.txt
/home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/est-sfs /home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/config-kimura.txt rh_synonymous_input2.txt /home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/seedfile.txt rh.synonymous.est_sfs.kimura.output.txt rh.synonymous.est_sfs.kimura.output.p_anc.txt
/home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/est-sfs /home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/config-kimura.txt rh_tolerated_input2.txt /home/wangcy/data/snpdmc/213/geneticload2/est-sfs/est-sfs-release-2.04/seedfile.txt rh.tolerated.est_sfs.kimura.output.txt rh.tolerated.est_sfs.kimura.output.p_anc.txt





perl est_sfs_ancestral_derived.pl rh.deleterious.est_sfs.frq rh.deleterious.est_sfs.output.p_anc.txt 
perl est_sfs_ancestral_derived.pl rh.loss_of_function.est_sfs.frq rh.loss_of_function.est_sfs.output.p_anc.txt
perl est_sfs_ancestral_derived.pl rh.synonymous.est_sfs.frq rh.synonymous.est_sfs.output.p_anc.txt
perl est_sfs_ancestral_derived.pl rh.tolerated.est_sfs.frq rh.tolerated.est_sfs.output.p_anc.txt
##
##for each indivdiual, extract the homozygotes and heterzygotes snps for various snp types
#the output format can be: individual, derived_heterozygotes, derived_homozygotes
awk '$3=="ALT"' rh.deleterious.est_sfs.ances_derive.txt |cut -f 1,2 |sed 's/\t/:/g' > rh_derived_hom_het/rh.deleterious.derived.snps
awk '$3=="ALT"' rh.loss_of_function.est_sfs.ances_derive.txt |cut -f 1,2 |sed 's/\t/:/g' > rh_derived_hom_het/rh.loss_of_function.derived.snps
awk '$3=="ALT"' rh.synonymous.est_sfs.ances_derive.txt |cut -f 1,2 |sed 's/\t/:/g' > rh_derived_hom_het/rh.synonymous.derived.snps
awk '$3=="ALT"' rh.tolerated.est_sfs.ances_derive.txt |cut -f 1,2 |sed 's/\t/:/g' > rh_derived_hom_het/rh.tolerated.derived.snps

vcftools --gzvcf rh.deleterious.outgroup.vcf.gz --snps rh.deleterious.derived.snps --recode --recode-INFO-all --out rh.deleterious.derived
vcftools --gzvcf rh.loss_of_function.outgroup.vcf.gz --snps rh.loss_of_function.derived.snps --recode --recode-INFO-all --out rh.loss_of_function.derived
vcftools --gzvcf rh.synonymous.outgroup.vcf.gz --snps rh.synonymous.derived.snps --recode --recode-INFO-all --out rh.synonymous.derived
vcftools --gzvcf rh.tolerated.outgroup.vcf.gz --snps rh.tolerated.derived.snps --recode --recode-INFO-all --out rh.tolerated.derived

##Calculate genotype counts for each sample, including the number of (homozygous) and (heterozygous) per sample
plink2 --vcf rh.deleterious.derived.recode.vcf --sample-counts --allow-extra-chr --out rh.deleterious.derived.hom_het
plink2 --vcf rh.loss_of_function.derived.recode.vcf --sample-counts --allow-extra-chr --out rh.loss_of_function.derived.hom_het
plink2 --vcf rh.synonymous.derived.recode.vcf --sample-counts --allow-extra-chr --out rh.synonymous.derived.hom_het
plink2 --vcf rh.tolerated.derived.recode.vcf --sample-counts --allow-extra-chr --out rh.tolerated.derived.hom_het


## Plotting
# Correlation scatter plot with line for GF-RDA-load

library(ggplot2)
library(ggpubr)  # Load package
library(ggplot2)
library(ggpubr)
data <- read.csv("dataset.csv")

# Plot for Deleterious/Syn Ratio vs GO_5085
ggplot(data, mapping = aes(x = GO_5085, y = average_deleterious)) + 
  geom_point(size = 5, shape = 16, color = "blue") +  # Solid blue circles
  geom_smooth(method = 'lm', formula = 'y ~ x', color = "black",  # Change line color to black
              fill = "gray", alpha = 0.5) +  # Change fill color to gray
  stat_cor(method = 'spearman', 
           label.sep = "\n", 
           size = 5, 
           digits = 2) +  # Set digits to 2 for r and p values
  theme(panel.grid = element_blank(),  # Hide grid lines
        panel.background = element_rect(fill = 'white', color = 'black'),  # Background settings
        plot.background = element_rect(fill = 'white', color = NA),  # Entire plot background
        text = element_text(size = 12)) +  # Set font size to 12
  labs(y = "Deleterious/Syn Ratio")  # Y-axis label

# Plot for Lof/Syn Ratio vs GO_5085
ggplot(data, mapping = aes(x = GO_5085, y = average_loss_of_function)) +
  geom_point(size = 5, shape = 16, color = "blue") +  # Solid blue circles
  geom_smooth(method = 'lm', formula = 'y ~ x', color = "black",  # Change line color to black
              fill = "gray", alpha = 0.5) +  # Change fill color to gray
  stat_cor(method = 'spearman', 
           label.sep = "\n", 
           size = 5, 
           digits = 2) +  # Set digits to 2 for r and p values
  theme(panel.grid = element_blank(),  # Hide grid lines
        panel.background = element_rect(fill = 'white', color = 'black'),  # Background settings
        plot.background = element_rect(fill = 'white', color = NA),  # Entire plot background
        text = element_text(size = 12)) +  # Set font size to 12
  labs(y = "Lof/Syn Ratio")  # Y-axis label

# Plot for Tolerated/Syn Ratio vs GO_5085
ggplot(data, mapping = aes(x = GO_5085, y = average_tolerated)) +
  geom_point(size = 5, shape = 16, color = "blue") +  # Solid blue circles
  geom_smooth(method = 'lm', formula = 'y ~ x', color = "black",  # Change line color to black
              fill = "gray", alpha = 0.5) +  # Change fill color to gray
  stat_cor(method = 'spearman', 
           label.sep = "\n", 
           size = 5, 
           digits = 2) +  # Set digits to 2 for r and p values
  theme(panel.grid = element_blank(),  # Hide grid lines
        panel.background = element_rect(fill = 'white', color = 'black'),  # Background settings
        plot.background = element_rect(fill = 'white', color = NA),  # Entire plot background
        text = element_text(size = 12)) +  # Set font size to 12
  labs(y = "Tolerated/Syn Ratio")  # Y-axis label

# Plot for Pi vs GO_5085
ggplot(data, mapping = aes(x = GO_5085, y = pi)) +
  geom_point(size = 5, shape = 16, color = "blue") +  # Solid blue circles
  geom_smooth(method = 'lm', formula = 'y ~ x', color = "black",  # Change line color to black
              fill = "gray", alpha = 0.5) +  # Change fill color to gray
  stat_cor(method = 'spearman', 
           label.sep = "\n", 
           size = 5, 
           digits = 2) +  # Set digits to 2 for r and p values
  theme(panel.grid = element_blank(),  # Hide grid lines
        panel.background = element_rect(fill = 'white', color = 'black'),  # Background settings
        plot.background = element_rect(fill = 'white', color = NA),  # Entire plot background
        text = element_text(size = 12)) +  # Set font size to 12
  labs(y = "π")  # Y-axis label

	

