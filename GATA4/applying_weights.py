from bs4 import BeautifulSoup
import json
import re
import glob
from collections import Counter
import numpy as np

file_pattern = 'bin_*/*.fasta_8mers/meme.html'
files = sorted(glob.glob(file_pattern, recursive=True), key=lambda x: int(re.search(r'bin_(\d+)', x).group(1)))

pwm_sections = []
consensus_sequences = []

def calculate_consensus(position, pwm, base):    #return percentage of base in a particular position
    for position in pwm:
        print(position)
        max_index = position.index(max(position))


def read_data(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    data = [line.strip().split('\t') for line in lines[1:]]
    return data

def bin_escores(data, bin_percentage):
    escores = np.array([float(row[2]) for row in data])
    bin_count = int(100 / bin_percentage)
    bin_edges = np.linspace(escores.max(), escores.min(), bin_count + 1)
    binned_escores = np.digitize(escores, bin_edges, right=True)
    return escores, binned_escores, bin_edges

def calculate_stats(escores, binned_escores, bin_edges):
    stats = []
    for i in range(len(bin_edges)-1):
        indices = np.where(binned_escores == i + 1)[0]
        bin_escores = escores[indices]
        avg = np.mean(bin_escores) if bin_escores.size > 0 else 0
        std = np.std(bin_escores) if bin_escores.size > 0 else 0
        stats.append((bin_edges[i], bin_edges[i+1], avg, std))
    return stats

def main():
    file_path = 'GATA4/GATA4_anti-GST/GATA4_anti-GST_8mers_top_enrichment.txt'
    bin_percentage = 10
    data = read_data(file_path)
    escores, binned_escores, bin_edges = bin_escores(data, bin_percentage)
    stats = calculate_stats(escores, binned_escores, bin_edges)
    
    print(f"Binning {i}%")
    metrics = []
    weights = []
    for i, (bin_start, bin_end, avg, std) in enumerate(stats):
        metric = avg/std
        metrics.append(metric)
        
    for i, (bin_start, bin_end, avg, std) in enumerate(stats):
        weight = metrics[i]/sum(metrics)
        weights.append(weight)
        print(f"Bin {i+1}: {bin_start:.4f} to {bin_end:.4f} - Average: {avg:.4f}, Std Dev: {std:.4f}, Weight: {weight:.4f}")
    
    print("***************")
        
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

    for i in range(0,9):
        for base in range(0,9):
            for j in pwm_sections:
                consensus = calculate_consensus(i,j, base) #get percentage of base in a particular position
                consensus_sequences.append(consensus)   #multiply weights and 1/10
                print("**************")

    print("\n")

if __name__ == "__main__":
    main()




