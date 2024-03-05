import csv

# Path to your CSV file
csv_file_path = 'GATA4_anti-GST/GATA4_anti-GST_8mers_top_enrichment.txt'
# Path to the output FASTA file
fasta_file_path = 'output_sequences.fasta'

with open(csv_file_path, mode='r') as csv_file:
    # Create a CSV reader object
    csv_reader = csv.reader(csv_file, delimiter='\t')  # Assuming tab-separated values; adjust if necessary
    with open(fasta_file_path, mode='w') as fasta_file:
        for i, row in enumerate(csv_reader):
            sequence = row[0].replace(".", "")  # Extract the sequence and remove any dots
            # Write to the FASTA file
            fasta_file.write(f'>Sequence{i+1}\n{sequence}\n')
