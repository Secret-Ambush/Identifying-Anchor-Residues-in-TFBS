from Bio import SeqIO

def split_into_bins(fasta_file, bin_size_percentage=50):
    # Read sequences from the FASTA file
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    total_sequences = len(sequences)
    
    # Calculate the number of sequences per bin
    sequences_per_bin = max(int(bin_size_percentage / 100 * total_sequences), 1)
    
    # Create bins and write to separate files
    output_files = []
    seq_count = 0  
    for i in range(0, total_sequences, sequences_per_bin):
        bin_sequences = sequences[i:i + sequences_per_bin]
        output_filename = f"bin_{i//sequences_per_bin + 1}.fasta"
        
        # Reset sequence numbering for each new file
        for seq_record in bin_sequences:
            seq_count += 1
            seq_record.id = f"Sequence{seq_count}"
            seq_record.description = ""
        seq_count = 0
        
        SeqIO.write(bin_sequences, output_filename, "fasta")
        output_files.append(output_filename)
    
    return output_files

fasta_file = '/Users/bristi/Desktop/Design Project/Binding sites/GATA4/output_sequences.fasta'
output_files = split_into_bins(fasta_file, 50)  
print("Sequences have been divided into the following files:", output_files)
