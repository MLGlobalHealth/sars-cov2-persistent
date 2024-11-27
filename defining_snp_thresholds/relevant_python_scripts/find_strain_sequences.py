import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Load the strain pairs TXT file
strain_pairs = pd.read_table("strain_pairs.txt")

# Path to the large FASTA file
fasta_path = 'combined_fasta.fasta'

# Function to extract sequences
def extract_sequences(fasta_path, strain_list):
    sequences = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        if record.id in strain_list:
            sequences[record.id] = str(record.seq)
    return sequences

# Function to save sequences into a FASTA file
def save_fasta(sequences, strains, filename):
    records = [SeqRecord(Seq(sequences[s]), id=s, description="") for s in strains if s in sequences]
    SeqIO.write(records, filename, "fasta")

# Group the data by variant and process each variant separately
for variant in strain_pairs['variant'].unique():
    # Filter strains by variant
    variant_strains = pd.concat([
        strain_pairs.loc[strain_pairs['variant'] == variant, 'strain_1'],
        strain_pairs.loc[strain_pairs['variant'] == variant, 'strain_2']
    ]).unique()
    
    # Extract sequences for these strains
    sequences = extract_sequences(fasta_path, variant_strains)
    
    # Save the sequences into a FASTA file
    save_fasta(sequences, variant_strains, f"{variant}_sequences.fasta")
