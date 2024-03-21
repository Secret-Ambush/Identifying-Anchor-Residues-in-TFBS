from bs4 import BeautifulSoup
import json
import re
import glob
from collections import Counter

file_pattern = 'bin_*/*.fasta_4mers/meme.html'
files = sorted(glob.glob(file_pattern, recursive=True), key=lambda x: int(re.search(r'bin_(\d+)', x).group(1)))

pwm_sections = []
consensus_sequences = []

def calculate_consensus(pwm):
    consensus = ''
    residues = ['A', 'C', 'G', 'T']
    
    for position in pwm:
        max_index = position.index(max(position))
        consensus += residues[max_index]
    
    return consensus

for file_path in files:
    with open(file_path, 'r') as file:
        html_content = file.read()

    soup = BeautifulSoup(html_content, 'html.parser')

    script_tags = soup.find_all("script")
    data_script = None
    for script in script_tags:
        if 'var data' in script.text:
            data_script = script.text
            break

    if data_script:
        json_str_match = re.search(r'var data = ({.*?});', data_script, re.DOTALL)
        if json_str_match:
            json_str = json_str_match.group(1)
            data = json.loads(json_str)
            if 'motifs' in data and data['motifs']:
                pwm_section = data['motifs'][0].get('pwm', 'PWM data not available')
                pwm_sections.append(pwm_section)
                print(f'Analysing PWM from {file_path}:')
            else:
                print(f"No 'motifs' data found in {file_path}.")
        else:
            print(f"JSON data not found in the script tag of {file_path}.")
    else:
        print(f"Script tag containing 'var data' was not found in {file_path}.")

for i in pwm_sections:
    consensus = calculate_consensus(i)
    consensus_sequences.append(consensus)
    
with open('8mer_motif_total_10_consensus_analysis.txt', 'w') as output_file:
    output_file.write("8mer motif - Discrete 10% binning of ENTIRE dataset\n\n")
    n = 1
    output_file.write("Consesus sequences of bins: \n")
    for i in consensus_sequences:
        output_file.write(f"Bin {n}: {i} \n")
        n += 1 
    
    output_file.write("\n\n")    
    output_file.write("Analysis\n") 
    for position in range(len(consensus_sequences[0])):
        residues_at_position = [seq[position] for seq in consensus_sequences]
        count = Counter(residues_at_position)
        
        most_common_residue, most_common_count = count.most_common(1)[0]
        
        output_file.write(f"Position {position+1}:\n")
        for residue, residue_count in count.items():
            output_file.write(f"  {residue}: {residue_count} times\n")
        output_file.write(f"Most common: {most_common_residue} appears {most_common_count} times\n\n")
