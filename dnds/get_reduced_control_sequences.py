# Read the unique strains from the text file
with open("reduced_control_strains.txt", "r") as file:
    reduced_control_strains = {line.strip() for line in file}  # Use a set for faster lookups

# Input and output FASTA files
input_fasta = ""
output_fasta = "reduced_control_sequences.fasta"

# Open the output FASTA file in write mode
with open(output_fasta, "w") as out_fasta:
    # Read the input FASTA file
    with open(input_fasta, "r") as fasta_file:
        current_sequence = ""
        current_header = ""

        for line in fasta_file:
            line = line.strip()

            # If the line starts with '>', it's a header line
            if line.startswith(">"):
                # If we already have a sequence, process the previous one
                if current_header and current_sequence:
                    strain = current_header.split()[0]  # Extract the first word from the header as the strain
                    if strain in reduced_control_strains:
                        out_fasta.write(f">{current_header}\n")
                        out_fasta.write(f"{current_sequence}\n")

                # Reset for the new sequence
                current_header = line[1:]  # Remove the '>' from the header
                current_sequence = ""
            else:
                # Otherwise, this is a sequence line
                current_sequence += line

        # Don't forget to process the last sequence after the loop
        if current_header and current_sequence:
            strain = current_header.split()[0]
            if strain in reduced_control_strains:
                out_fasta.write(f">{current_header}\n")
                out_fasta.write(f"{current_sequence}\n")

print(f"Filtered sequences written to {output_fasta}")
