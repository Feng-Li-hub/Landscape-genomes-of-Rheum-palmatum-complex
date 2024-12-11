#!/bin/bash

# Loop through each population name listed in pop.tsv
for popname in $(cat ./pop.tsv)
do
    # Use bcftools to list the sample names from the VCF file, then filter by the current population name
    bcftools query -l /Archive/data/Rheum/RPC_Resequencing/final_vcf/top11chrAll.SNPable.filter.snp.vcf.gz | grep "$popname" > $popname.txt

    date

 
    echo $popname
done
echo "done"
date
exit 0

