import os

# Function to process each .hwe file and compute the average of the 7th and 8th columns
def process_hwe_file(file_name):
    sum_O_HET = 0  # Variable to store the sum of observed heterozygosity (O(HET))
    sum_E_HET = 0  # Variable to store the sum of expected heterozygosity (E(HET))
    line_count = 0  # Counter for valid data lines

    with open(file_name, 'r') as f:
        # Skip the header line of the file
        header = f.readline()
        
        # Loop through the remaining lines of the file
        for line in f:
            cols = line.strip().split()  # Split each line into columns
            try:
                # Try to retrieve and convert the 7th and 8th columns (O(HET) and E(HET)) into floats
                O_HET = float(cols[6])
                E_HET = float(cols[7])

                # If O(HET) or E(HET) is not a valid number between 0 and 1, skip this line
                if not (0 <= O_HET <= 1) or not (0 <= E_HET <= 1):
                    continue

                # Add valid O(HET) and E(HET) to the sums, and increment the line count
                sum_O_HET += O_HET
                sum_E_HET += E_HET
                line_count += 1
            except ValueError:
                # If the conversion to float fails, skip this line
                continue

    # If there were valid lines, calculate and return the averages of O(HET) and E(HET)
    if line_count > 0:
        avg_O_HET = sum_O_HET / line_count
        avg_E_HET = sum_E_HET / line_count
        return avg_O_HET, avg_E_HET
    else:
        # If no valid data, return None for both averages
        return None, None

# Function to process all .hwe files in the current directory and save the results to a summary file
def process_all_files(output_file):
    with open(output_file, 'w') as result_file:
        # Loop through all files in the current directory
        for file_name in os.listdir('.'):
            # Only process .hwe files
            if file_name.endswith('.hwe'):
                avg_O_HET, avg_E_HET = process_hwe_file(file_name)
                
                # If valid averages are calculated, write them to the result file
                if avg_O_HET is not None and avg_E_HET is not None:
                    result_file.write(f"{file_name}\t{avg_O_HET}\t{avg_E_HET}\n")
                    print(f"Processed {file_name}: O(HET) = {avg_O_HET}, E(HET) = {avg_E_HET}")
                else:
                    # If no valid data, print a message
                    print(f"Processed {file_name}: No valid data.")


if __name__ == "__main__":
    # Define the output file name where results will be saved
    output_file = "hwe_summary_results.txt"
    # Process all files and generate the summary
    process_all_files(output_file)
    # Notify that all results have been saved
    print(f"All results have been saved to {output_file}")