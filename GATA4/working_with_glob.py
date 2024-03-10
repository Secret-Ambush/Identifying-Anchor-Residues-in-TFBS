import glob
import os

pattern = 'binning_bin*_*/*.fasta'
fasta_files = glob.glob(pattern)
output_files = [os.path.join(os.path.basename(os.path.dirname(f)), os.path.basename(f)) for f in fasta_files]

print(output_files)
