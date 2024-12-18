import os

# Paths
input_fasta = ""
sequence_list_folder = "./sequence_names/"
output_folder = "./case_pair_trees/"

# Ensure output directory exists
os.makedirs(output_folder, exist_ok=True)

# Function to parse a FASTA file and extract specific sequences
def extract_sequences(fasta_path, sequence_ids):
    sequences = {}
    with open(fasta_path, "r") as fasta_file:
        current_id = None
        current_seq = []
        for line in fasta_file:
            if line.startswith(">"):
                # Save the current sequence before moving to the next
                if current_id and current_seq:
                    sequences[current_id] = "".join(current_seq)
                # Start a new sequence
                current_id = line[1:].strip()  # Remove ">" and whitespace
                current_seq = []
            else:
                # Add line to the current sequence
                current_seq.append(line.strip())
        # Save the last sequence in the file
        if current_id and current_seq:
            sequences[current_id] = "".join(current_seq)
    # Return only the sequences with matching IDs
    return {seq_id: sequences[seq_id] for seq_id in sequence_ids if seq_id in sequences}

# Process each event_id file
for file in os.listdir(sequence_list_folder):
    if file.endswith("_sequences.txt"):
        event_id = file.split("_sequences.txt")[0]
        # Read the sequence IDs for the current event_id
        with open(os.path.join(sequence_list_folder, file), "r") as f:
            sequence_ids = set(line.strip() for line in f)

        # Extract the sequences
        extracted_sequences = extract_sequences(input_fasta, sequence_ids)

        # Write to a new FASTA file
        output_fasta = os.path.join(output_folder, f"{event_id}.fasta")
        with open(output_fasta, "w") as out_fasta:
            for seq_id, sequence in extracted_sequences.items():
                out_fasta.write(f">{seq_id}\n")
                out_fasta.write(f"{sequence}\n")

