import os
from Bio import SeqIO


def split_into_bins(fasta_file, bin_size_percentage=50):
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    total_sequences = len(sequences)
    
    sequences_per_bin = max(int(bin_size_percentage / 100 * total_sequences), 1)
    
    output_folders = []
    seq_count = 0  
    for i in range(0, total_sequences, sequences_per_bin):
        bin_sequences = sequences[i:i + sequences_per_bin]
        output_folder = f"mutated_bin_{i//sequences_per_bin + 1}"
        os.makedirs(output_folder, exist_ok=True)
        
        for seq_record in bin_sequences:
            seq_count += 1
            seq_record.id = f"Sequence{seq_count}"
            seq_record.description = ""
        seq_count = 0
        
        output_filename = os.path.join(output_folder, f"mutated_bin_{i//sequences_per_bin + 1}.fasta")
        SeqIO.write(bin_sequences, output_filename, "fasta")
        output_folders.append(output_folder)
    
    return output_folders

fasta_file = '/Users/bristi/Desktop/Design Project/Working-with-TF/HOXD13_I297V/mutated_HOXD13_I297V_output_sequences.fasta'
output_folders = split_into_bins(fasta_file, 10)  

print("Sequences have been divided into the following folders:", output_folders)
