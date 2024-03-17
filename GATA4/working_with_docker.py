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
                "-maxw", "8", 
                "-minw", "4", 
                "-nmotifs", "1", 
                "-mod", "zoops", 
                "-objfun", "classic", 
                "-revcomp", 
                "-markov_order", "0", 
                "-o", f"/data/{input_basename}_results_new"],
               check=True)
        except subprocess.CalledProcessError as e:
            print(f"Command failed with error: {e}")

pattern = 'bin_*/*.fasta'
fasta_files = glob.glob(pattern)
output_files = [os.path.join(os.path.basename(os.path.dirname(f)), os.path.basename(f)) for f in fasta_files]

meme_analysis(output_files)
