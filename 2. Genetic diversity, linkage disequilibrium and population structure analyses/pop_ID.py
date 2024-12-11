# pop_ID.py

from collections import defaultdict

# Read the original pop.list file
with open("pop.list", "r") as f:
    lines = f.readlines()

# Create a dictionary where the key is FID, and the value is the corresponding list of samples
fid_dict = defaultdict(list)

# Process each line, grouping by FID
for line in lines:
    fid, iid = line.strip().split()  # Split each line into FID and IID
    fid_dict[fid].append(line.strip())  # Add the line to the list of samples for that FID

# Write the data for each FID to a separate file
for fid, data in fid_dict.items():
    with open(f"{fid}.list", "w") as f:
        f.write("\n".join(data) + "\n")  # Write each sample for this FID to the file

print("Splitting complete, list files for each FID have been generated.")


