# Loop through each population name in the 'pop.tsv' file
for popname in $(cat ./pop.tsv)
do
    # Extract the information for each population from the VCF file based on IDs in a population-specific file ($popname.txt)
    # --gzvcf: Specifies the input compressed VCF file
    # --keep: Keeps only the individuals listed in the $popname.txt file
    # --recode: Outputs the data in VCF format
    # --recode-INFO-all: Rewrites all INFO fields in the VCF
    # --out: Specifies the output file name prefix (e.g., p[popname].vcf)
    vcftools --gzvcf /Archive/data/Rheum/RPC_Resequencing/final_vcf/top11chrAll.SNPable.filter.snp.vcf.gz \
             --keep $popname.txt \
             --recode \
             --recode-INFO-all \
             --out p$popname

    # Calculate the PI (nucleotide diversity) for each population from the extracted VCF file
    # --vcf: Specifies the input VCF file (output of the previous vcftools command)
    # --out: Specifies the output file prefix for PI results (e.g., q[popname]_pi_500kb)
    # --window-pi: Sets the window size for calculating PI (500,000 base pairs in this case)
    # --window-pi-step: Sets the step size for the sliding window (50,000 base pairs in this case)
    vcftools --vcf p$popname.recode.vcf \
             --out q$popname\_pi_500kb \
             --window-pi 500000 \
             --window-pi-step 50000

    # Print the current date and time to show the script's progress
    date
    
    # Print the current population name being processed for tracking
    echo $popname
done
echo "done"
date
exit 0

