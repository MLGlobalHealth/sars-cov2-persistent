import os
import subprocess

# Define the path to the IQ-TREE executable
iqtree_path = ""

# Ensure the IQ-TREE executable exists at the specified path
if not os.path.isfile(iqtree_path):
    raise FileNotFoundError(f"IQ-TREE executable not found at {iqtree_path}. Please check the path.")

# Define the model to use
model = "GTR+F"

# Define the current working directory
main_dir = os.getcwd()

# Define the outgroup file
outgroup_file = os.path.join(main_dir, "hCoV-19_Wuhan_WIV04_2019.fasta")
if not os.path.isfile(outgroup_file):
    raise FileNotFoundError(f"Outgroup file not found at {outgroup_file}. Please ensure it exists.")

# Extract the taxon name and sequence from the outgroup file
with open(outgroup_file, 'r') as og_file:
    outgroup_lines = og_file.readlines()
    outgroup_taxon = outgroup_lines[0].strip().lstrip('>')
    outgroup_sequence = ''.join(line.strip() for line in outgroup_lines[1:])

# Define the directory containing the case pair trees
case_pair_trees_dir = os.path.join(main_dir, "case_pair_trees")
if not os.path.isdir(case_pair_trees_dir):
    raise NotADirectoryError(f"Case pair trees directory not found at {case_pair_trees_dir}. Please ensure it exists.")

# Function to prepend the outgroup to a FASTA file
def add_outgroup_to_fasta(fasta_path, outgroup_taxon, outgroup_sequence):
    with open(fasta_path, 'r') as fasta_file:
        original_content = fasta_file.read()

    with open(fasta_path, 'w') as fasta_file:
        # Write the outgroup as the first taxon
        fasta_file.write(f">{outgroup_taxon}\n{outgroup_sequence}\n")
        # Append the original content
        fasta_file.write(original_content)

# Loop through each FASTA file in the case pair trees directory
for fasta_file in os.listdir(case_pair_trees_dir):
    if fasta_file.endswith(".fasta"):
        fasta_path = os.path.join(case_pair_trees_dir, fasta_file)
        print(f"Adding outgroup to: {fasta_path}")

        # Add the outgroup to the FASTA file
        add_outgroup_to_fasta(fasta_path, outgroup_taxon, outgroup_sequence)

        print(f"Running IQ-TREE on: {fasta_path}")
        
        # Run IQ-TREE with the specified parameters
        try:
            subprocess.run([
                iqtree_path,
                "-s", fasta_path,         # Input FASTA file
                "-m", model,             # Model to use
                "-nt", "16",            # Number of threads
                "-o", outgroup_taxon,   # Specify the outgroup
                "-redo"                  # Overwrite previous results
            ], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error occurred while processing {fasta_path}: {e}")

print("IQ-TREE processing completed for all FASTA files.")
