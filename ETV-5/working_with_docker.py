import os
import subprocess
import glob

def meme_analysis(input_files):
    if not input_files:
        print("No input files provided.")
        return
    
    for input_file in input_files:
        print(f"Processing file: {input_file}")
        input_abs_path = os.path.abspath(input_file)
        input_dir, input_basename = os.path.split(input_abs_path)
        print(input_basename)
        if not os.path.exists(input_dir):
            print(f"Directory does not exist: {input_dir}")
            continue
        
        print(f"Executing Docker command for: {input_file}")
        try:
            subprocess.run(["docker", "run", "--rm", 
                "-v", f"{input_dir}:/data",  
                "memesuite/memesuite:latest", 
                "meme", f"/data/{input_basename}", 
                "-dna", 
                "-o",
                "-nostatus",
                "-maxw", "6", 
                "-minw", "6", 
                "-nmotifs", "1", 
                "-mod", "zoops", 
                "-objfun", "classic", 
                "-revcomp", 
                "-markov_order", "0", 
                "-o", f"/data/{input_basename}_6mers"],
               check=True)
        except subprocess.CalledProcessError as e:
            print(f"Command failed with error: {e}")

pattern = 'ETV-5/bin_*/*.fasta'
fasta_files = glob.glob(pattern)

# Use absolute paths for input files
absolute_fasta_files = [os.path.abspath(f) for f in fasta_files]

meme_analysis(absolute_fasta_files)