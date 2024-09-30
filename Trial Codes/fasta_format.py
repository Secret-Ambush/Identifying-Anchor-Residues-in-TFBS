import csv

csv_file_path = '/Users/bristi/Desktop/Design Project/Working-with-TF/TFs UniProbe Data/GATA4/GATA4_anti-GST_8mers_top_enrichment.txt'
fasta_file_path = 'gata(2).fasta'
            
sequences = []
with open(csv_file_path, 'r') as input_file, open(fasta_file_path, 'w') as output_file:
        next(input_file)
        for i, line in enumerate(input_file, start=1):
            sequence = line.strip().split('\t')[0]
            sequence = sequence.replace('.','')
            sequences.append(sequence)
        
        sequences = set(sequences)
        for i in sequences:      
            output_file.write(f"{i}\n")
