from collections import Counter

consensus_sequences = ["ATGC", "ATGT", "ATGG", "ATGA", "ATGC"]

for position in range(len(consensus_sequences[0])):
    residues_at_position = [seq[position] for seq in consensus_sequences]
    count = Counter(residues_at_position)
    most_common = count.most_common(1)[0] 
    
    threshold = 0.8 * len(consensus_sequences)
    if most_common[1] >= threshold:
        print(f"Position {position+1}: Anchor residue {most_common[0]} appears {most_common[1]} times")
    else:
        print(f"Position {position+1}: No anchor residue meets the threshold.")
