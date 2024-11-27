from Bio import SeqIO

# Load the list of strains
with open("all_possible_persistent_infections_strain_list_without_CT.txt") as f:
    strain_set = set(line.strip() for line in f)

# Define the input and output FASTA files
input_fasta = "combined_fasta.fasta"
output_fasta = "all_possible_persistent_infections_persistent_sequences_without_CT.fasta"

# Extract relevant sequences
with open(output_fasta, "w") as output_handle:
    for record in SeqIO.parse(input_fasta, "fasta"):
        if record.id in strain_set:
            SeqIO.write(record, output_handle, "fasta")

print(f"Extracted sequences saved to {output_fasta}")
