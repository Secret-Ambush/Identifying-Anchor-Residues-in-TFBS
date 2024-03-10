from bs4 import BeautifulSoup
import json
import re
import glob
from collections import Counter

# Use glob.glob to find all files that match the pattern
file_pattern = 'binning_1_bin_*/binning_1_bin_*.fasta_results/meme.html'  # Assuming your HTML files have a .html extension
files = glob.glob(file_pattern, recursive=True)

pwm_sections = []
consensus_sequences = []

def calculate_consensus(pwm):
    """
    Calculate the consensus sequence for a given PWM.
    
    Args:
    pwm (list of lists): The PWM as a list of lists, where each inner list represents
                         probabilities of residues at that position.
    
    Returns:
    str: The consensus sequence.
    """
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
                print(f'Found PWM from {file_path}:')
            else:
                print(f"No 'motifs' data found in {file_path}.")
        else:
            print(f"JSON data not found in the script tag of {file_path}.")
    else:
        print(f"Script tag containing 'var data' was not found in {file_path}.")

for i in pwm_sections:
    consensus = calculate_consensus(i)
    consensus_sequences.append(consensus)
    
for position in range(len(consensus_sequences[0])):
    residues_at_position = [seq[position] for seq in consensus_sequences]
    count = Counter(residues_at_position)
    most_common = count.most_common(1)[0]
    
    threshold = 0.8 * len(consensus_sequences)
    if most_common[1] >= threshold:
        print(f"Position {position+1}: Anchor residue {most_common[0]} appears {most_common[1]} times")
    else:
        print(f"Position {position+1}: No anchor residue meets the threshold.")