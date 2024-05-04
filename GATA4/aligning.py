import pandas as pd

def read_sequences_from_file(file_path):
    sequences = []
    with open(file_path, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                sequences.append(line.strip())
    return sequences

def align_sequences_to_table(sequences, base_sequence):
    # Calculate the start position of 'GATAA' in the base sequence
    base_index = base_sequence.find('GATAA')
    if base_index == -1:
        raise ValueError("The base sequence does not contain 'GATAA'.")

    # Create a DataFrame to hold the alignment
    columns = [i for i in range(-7, 15)]  # From -7 to +14
    df = pd.DataFrame(columns=columns)
    
    # Add the base sequence to the DataFrame first
    df.loc[0] = pd.Series(list(base_sequence), index=[x + base_index for x in range(8)])

    # Align other sequences
    for sequence in sequences:
        current_index = sequence.find('GATAA')
        if current_index != -1:
            # Create a series from the sequence with the correct offset
            offset = base_index - current_index
            series_index = [x + offset for x in range(8)]
            new_row = pd.Series(list(sequence), index=series_index)
            df = df.append(new_row, ignore_index=True)

    return df

# Paths to the files
input_file_path = 'bin_1/bin_1.fasta'
output_file_path = 'aligned_sequences.csv'

# Read the base sequence and other sequences
base_sequences = read_sequences_from_file(input_file_path)
if not base_sequences:
    raise ValueError("No sequences found in the base sequence file.")
base_sequence = base_sequences[0]  # First sequence is the base

# Read other sequences
sequences = read_sequences_from_file(input_file_path)

# Align the sequences
df_aligned = align_sequences_to_table(sequences, base_sequence)

# Save the DataFrame to a CSV file
df_aligned.to_csv(output_file_path, index=False)

print(f"Sequences have been aligned and saved to {output_file_path}")
