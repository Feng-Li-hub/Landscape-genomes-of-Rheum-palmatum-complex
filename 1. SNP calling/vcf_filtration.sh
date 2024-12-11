# Convert VCF to PED and MAP format
plink -vcf top11chrAll.SNPable.filter.snp.vcf -recode -out 213snp1 --allow-extra-chr  # VCF to PED/MAP
plink -file 213snp1 -make-bed -out 213snp1 --allow-extra-chr  # PED/MAP to BIM/BED/FAM
# Check SNP and sample counts
wc -l 213snp1.map 213snp1.ped  # Check SNP and sample numbers （14299652，213）


# Check missing data (individuals and SNPs)
plink --file 213snp1 --missing --allow-extra-chr  # Missing data summary
# Outputs: plink.imiss (individual missing data), plink.lmiss (SNP missing data)

# Filter SNPs with >2% missing data
plink -bfile 213snp1 -geno 0.02 -make-bed -out 213snp2 --allow-extra-chr  # SNP missing >2% filter
plink -bfile 213snp2 -recode -out 213snp2 --allow-extra-chr  # Convert back to PED/MAP
wc -l 213snp2.map 213snp2.ped  # Check SNP and sample numbers （1432892，213）

# Filter samples with >2% missing data
plink -bfile 213snp2 -mind 0.02 -make-bed -out 213snp3 --allow-extra-chr  # Sample missing >2% filter
plink -bfile 213snp3 -recode -out 213snp3 --allow-extra-chr  # Convert back to PED/MAP
wc -l 213snp3.map 213snp3.ped  # Check SNP and sample numbers （1432892，201）
# SNP quality control: First filter SNPs, then filter samples. Run "--geno" first, then "--mind"

plink -bfile 213snp3 -freq -out MAF_check --allow-extra-chr  # Calculate MAF for each SNP
plink -bfile 213snp3 -maf 0.02 -make-bed -out 213snp4 --allow-extra-chr  # Filter SNPs with MAF < 0.02
plink -bfile 213snp4 -recode -out 213snp4 --allow-extra-chr  # Convert back to PED/MAP
wc -l 213snp4.map 213snp4.ped  # Check SNP and sample numbers （1432892，201）

# Perform LD (Linkage Disequilibrium) analysis: prune SNPs based on pairwise LD
plink --file 213snp3 --indep-pairwise 50 10 0.2 -aec --out 213snp_LD  # LD analysis with pruning (50 SNP window, 10 SNP step, r^2 threshold 0.2) 
# Extract the pruned SNPs and create a new binary file with only the selected SNPs
plink --file 213snp3 --extract 213snp_LD.prune.in --make-bed --out 213snp_LD --allow-extra-chr  # Extract pruned SNPs and convert to binary format (BED/BIM/FAM)
wc -l 213snp_LD.map 213snp_LD.ped  # Check SNP and sample numbers （378620，201）
# Convert the pruned data back to VCF format
plink --bfile 213snp_LD --recode vcf --out 213snp_LD --allow-extra-chr  # Convert binary format to VCF
# Convert the VCF file back to PED/MAP format
plink -vcf 213snp_LD.vcf -recode -out 213snp_LD --allow-extra-chr  # Convert VCF to PED/MAP format






