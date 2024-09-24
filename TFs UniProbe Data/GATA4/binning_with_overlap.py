import os
from Bio import SeqIO

def split_into_bins(fasta_file, bin_size_percentage=50, overlap_percentage=5):
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    total_sequences = len(sequences)
    
    sequences_per_bin = max(int(bin_size_percentage / 100 * total_sequences), 1)
    overlap_amount = max(int(overlap_percentage / 100 * sequences_per_bin), 1)
    
    output_folders = []
    seq_count = 0
    i = 0
    bin_count = 0
    while i < total_sequences:
        bin_sequences = sequences[i:min(i + sequences_per_bin, total_sequences)]
        bin_count += 1
        output_folder = f"binning_1_bin_{bin_count}"
        os.makedirs(output_folder, exist_ok=True)
        
        for seq_record in bin_sequences:
            seq_count += 1
            seq_record.id = f"Sequence{seq_count}"
            seq_record.description = ""
        
        seq_count = 0
        
        output_filename = os.path.join(output_folder, f"binning_1_bin_{bin_count}.fasta")
        SeqIO.write(bin_sequences, output_filename, "fasta")
        output_folders.append(output_folder)
        
        i += sequences_per_bin - overlap_amount
    
    return output_folders

fasta_file = '/Users/bristi/Desktop/Design Project/Working-with-TF/GATA4/bin_1/bin_1.fasta'
output_folders = split_into_bins(fasta_file, 10, 5)  

print("Sequences have been divided into the following folders:", output_folders)
