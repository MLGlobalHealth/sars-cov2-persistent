from Bio import SeqIO
import pandas as pd

# Step 1: Load the case-to-control mapping
case_to_controls = {}

with open("case_to_controls.txt", "r") as f:
    for line in f:
        parts = line.strip().split("\t")
        case_strain = parts[0]
        control_strains = parts[1:]
        case_to_controls[case_strain] = control_strains

# Step 2: Load the FASTA file and store sequences in a dictionary
fasta_file = ""
sequences = {}

# Parse the FASTA file and store the sequences in a dictionary
for record in SeqIO.parse(fasta_file, "fasta"):
    sequences[record.id] = str(record.seq)

# Helper function to calculate Hamming distance, considering only A, T, C, G, and lowercase
def hamming_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        return float('inf')  # In case the sequences are of different lengths

    # Convert sequences to uppercase for comparison
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    
    # Only count positions where both sequences have valid nucleotides (A, T, C, G)
    valid_bases = {'A', 'T', 'C', 'G'}
    differences = 0
    
    for base1, base2 in zip(seq1, seq2):
        if base1 in valid_bases and base2 in valid_bases:  # Ignore positions with N or any invalid characters
            if base1 != base2:
                differences += 1
    
    return differences

# Step 3: Compare case strains with their corresponding control strains
matches = []
zero_diff_cases = []
one_diff_cases = []
two_diff_cases = []

for case_strain, control_strains in case_to_controls.items():
    # If the case strain is in the FASTA file
    if case_strain in sequences:
        case_sequence = sequences[case_strain]
        
        for control_strain in control_strains:
            # If the control strain is in the FASTA file
            if control_strain in sequences:
                control_sequence = sequences[control_strain]
                
                # Calculate the Hamming distance between the case and control sequences
                diff_count = hamming_distance(case_sequence, control_sequence)
                
                # Classify based on the number of differences
                if diff_count == 0:
                    zero_diff_cases.append({
                        'case_strain': case_strain,
                        'control_strain': control_strain,
                        'differences': diff_count
                    })
                elif diff_count == 1:
                    one_diff_cases.append({
                        'case_strain': case_strain,
                        'control_strain': control_strain,
                        'differences': diff_count
                    })
                elif diff_count == 2:
                    two_diff_cases.append({
                        'case_strain': case_strain,
                        'control_strain': control_strain,
                        'differences': diff_count
                    })
                
                # Append the result (including matches and mismatches)
                matches.append({
                    'case_strain': case_strain,
                    'control_strain': control_strain,
                    'differences': diff_count
                })

# Step 4: Create DataFrames for each category of differences (0, 1, 2)
matches_df = pd.DataFrame(matches)
zero_diff_df = pd.DataFrame(zero_diff_cases)
one_diff_df = pd.DataFrame(one_diff_cases)
two_diff_df = pd.DataFrame(two_diff_cases)

# Step 5: Output the number of cases with 0, 1, or 2 differences
num_zero_diff_cases = len(zero_diff_df)
num_one_diff_cases = len(one_diff_df)
num_two_diff_cases = len(two_diff_df)

print(f"Number of cases with 0 differences: {num_zero_diff_cases}")
print(f"Number of cases with 1 difference: {num_one_diff_cases}")
print(f"Number of cases with 2 differences: {num_two_diff_cases}")

# Step 6: Save the DataFrames to CSV files
matches_df.to_csv("case_control_matches.csv", index=False)
zero_diff_df.to_csv("zero_diff_matches.csv", index=False)
one_diff_df.to_csv("one_diff_matches.csv", index=False)
two_diff_df.to_csv("two_diff_matches.csv", index=False)

print("CSV files have been saved successfully.")
