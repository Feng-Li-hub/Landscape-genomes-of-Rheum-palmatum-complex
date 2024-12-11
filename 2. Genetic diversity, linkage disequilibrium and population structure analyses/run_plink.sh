#!/bin/bash

# Loop through all .list files in the current directory
for list_file in *.list; do
    # Remove the .list extension from the file name to use it as the output file prefix
    output_prefix="${list_file%.list}"
    
    # Run PLINK analysis with the specified options
    # --file: Specifies the input file (top11chrAll)
    # --keep: Keeps only the individuals listed in the .list file
    # --hardy: Performs Hardy-Weinberg equilibrium testing
    plink --file top11chrAll --keep "$list_file" --hardy --out "$output_prefix" --allow-extra-chr
    
    echo "Analysis complete: $list_file -> $output_prefix"
done


