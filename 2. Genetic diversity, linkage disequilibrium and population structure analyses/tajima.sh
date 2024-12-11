# Loop through each population name (popname) read from the pop.tsv file
while read popname; do
    # Assume your VCF file names follow the format p<popname>.recode.vcf, e.g., pGSZN.recode.vcf
    vcf_file="p${popname}.recode.vcf"
    
    # Check if the VCF file exists
    if [ -f "$vcf_file" ]; then
        # Calculate Tajima's D using vcftools, and output the result to pop.tajimaD.txt
        # -c option calculates Tajima's D for each site, and 'sed' removes the header row
        # 'awk' averages the Tajima's D values and appends the result with the population name
        vcftools --vcf $vcf_file --TajimaD 100000 -c | sed '1d' | awk -v pop=$popname '{td+=\\$4};END{print pop"\t"td/NR};' >> pop.tajimaD.txt
    else
        # If the VCF file is not found, output a warning message to pop.tajimaD.txt
        echo "Warning: VCF file $vcf_file not found, skipping." >> pop.tajimaD.txt
    fi
done < ./pop.tsv

