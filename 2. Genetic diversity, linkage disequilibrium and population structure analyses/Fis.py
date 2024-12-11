import pandas as pd

# Read the data file (assuming the file uses whitespace as delimiter)
data = pd.read_csv("top11chrAll.het", delim_whitespace=True)

# Calculate expected heterozygosity He (assumed as He = 1 - E(HOM))
data['He'] = 1 - data['E(HOM)']

# Calculate observed heterozygosity Ho (Ho = (1 - F) * He)
data['Ho'] = (1 - data['F']) * data['He']

# Calculate FIS = (He - Ho) / He
data['FIS'] = (data['He'] - data['Ho']) / data['He']

# Calculate the average FIS for each population (grouped by FID)
average_FIS = data.groupby('FID')['FIS'].mean().reset_index()

# Print the result
print(average_FIS)

# Save the result to a file with tab-separated values
average_FIS.to_csv('average_fis_by_population.txt', sep='\t', index=False)


