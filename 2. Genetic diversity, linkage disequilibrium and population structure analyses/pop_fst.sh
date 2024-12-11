#!/bin/bash

# Define the list of groups
groups=("CQJFS" "GSC" "GSHC" "GSLN" "GSM" "GSSD" "GSSN" "GSTZ" "GSW" "GSYC" "GSZN" "GSZQ" "HBXS" "HNFNS" "HNLC" "QHDR" "QHMQ" "QHMY" "QHTJ" "QHYS" "QHZK" "SCBZ" "SCKD" "SCLD" "SCMG" "SCMN" "SCMX" "SCMYL" "SCPW" "SCSPS" "SCSPY" "SCWC" "SCXL" "SNH" "SNMW" "SNPL" "SNZZL" "SNZZN" "SXQS" "XZBQ" "XZNML" "XZQND" "YNXZD")

output_file="FST_matrix.txt"

# Initialize the header row for the FST output matrix
echo -n "Group" > $output_file  # Start with the first column as "Group"
for group in "${groups[@]}"; do
    echo -n -e "\t$group" >> $output_file  # Add each group name as a header in the matrix
done
echo "" >> $output_file  # Add a newline after the header to start the data rows


# Loop through each pair of groups
for i in ${!groups[@]}; do
    group1=${groups[$i]}  # Get the first group
    
    # Output the first group name in the first column of the output file
    echo -n -e "$group1" >> $output_file
    
    # Loop through the remaining groups (starting from the next one)
    for j in $(seq $((i+1)) ${#groups[@]}); do
        group2=${groups[$j]}  # Get the second group
        
        echo "Running vcftools for $group1 vs $group2"  # Print the current comparison
        
        # Run vcftools to calculate FST between the two groups
        vcftools --vcf /Archive/data/Rheum/RPC_Resequencing/final_vcf/top11chrAll.SNPable.filter.snp.vcf \
            --weir-fst-pop ${group1}.txt \
            --weir-fst-pop ${group2}.txt \
            --out ${group1}_vs_${group2}_FST
        
    done
    
    # Add a newline after each group comparison to start a new line in the output file
    echo "" >> $output_file
done

echo "vcftools FST calculations have been completed, but no further processing has been done."

