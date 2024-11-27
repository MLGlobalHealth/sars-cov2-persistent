from Bio import SeqIO
import pandas as pd

# Paths to your files
fasta_file = 'combined_fasta.fasta'
alpha_sequences_file = 'alpha_sequences.csv'
output_csv_file = 'alpha_nucleotide_counts_per_position.csv'

# Read the CSV file to get the list of sequence names
alpha_sequences_df = pd.read_csv(alpha_sequences_file)
alpha_sequences_set = set(alpha_sequences_df['strain'].tolist())

# Initialize a dictionary to store counts
nucleotides = ['A', 'T', 'C', 'G', 'N', '-', 'Other']
genome_length = 29891
position_counts = {nuc: [0] * genome_length for nuc in nucleotides}

# Define a mapping for lowercase to uppercase nucleotides
nucleotide_mapping = {'a': 'A', 't': 'T', 'c': 'C', 'g': 'G'}

# Read sequences from the FASTA file and update counts for the filtered sequences
for record in SeqIO.parse(fasta_file, "fasta"):
    if record.id in alpha_sequences_set:
        sequence = str(record.seq)
        for i, nucleotide in enumerate(sequence):
            if nucleotide in nucleotide_mapping:
                nucleotide = nucleotide_mapping[nucleotide]
            if nucleotide in position_counts:
                position_counts[nucleotide][i] += 1
            else:
                position_counts['Other'][i] += 1

# Create a DataFrame from the counts dictionary
df = pd.DataFrame(position_counts)
df.index.name = 'Position'
df.reset_index(inplace=True)

# Save the DataFrame to CSV
df.to_csv(output_csv_file, index=False)
