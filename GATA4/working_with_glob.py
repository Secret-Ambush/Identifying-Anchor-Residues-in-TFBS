import glob
import os

# Pattern to match files of the format "binning_binX_Y/binning_binX_Y.fasta"
pattern = 'binning_bin*_*/*.fasta'

# Use glob to find all files matching the pattern
fasta_files = glob.glob(pattern)

# Preparing the list to include only the needed format
output_files = [os.path.join(os.path.basename(os.path.dirname(f)), os.path.basename(f)) for f in fasta_files]

print(output_files)
