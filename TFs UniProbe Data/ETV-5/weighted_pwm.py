from bs4 import BeautifulSoup
import json
import re
import glob
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px

pwm_sections = []
consensus_sequences = []

def calculate_consensus(position, pwm, base):    
    consensus = ''
    residues = ['A', 'C', 'G', 'T']
    
    for position in pwm:
        max_index = position.index(max(position))
        consensus += residues[max_index]
    
    return consensus

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

def read_html_pwm(files):
    for file_path in files:
        with open(file_path, 'r') as file:
            html_content = file.read()
        soup = BeautifulSoup(html_content, 'html.parser')
        script_tags = soup.find_all("script")
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
                else:
                    print(f"No 'motifs' data found in {file_path}.")
            else:
                print(f"JSON data not found in the script tag of {file_path}.")
        else:
            print(f"Script tag containing 'var data' was not found in {file_path}.")
    
    return pwm_sections

def apply_weights(pwms, weights):
    weighted_pwms = []
    for pwm, weight in zip(pwms, weights):
        weighted_pwm = [[value * weight for value in position] for position in pwm]
        weighted_pwms.append(weighted_pwm)
    return weighted_pwms

def sum_pwms(weighted_pwms):
    sum_pwm = np.sum(np.array(weighted_pwms), axis=0)
    return sum_pwm

def main():
    file_path = 'ETV-5/Datasets/Etv5_8mers_top_enrichment.txt'
    bin_percentage = 10
    data = read_data(file_path)
    escores, binned_escores, bin_edges = bin_escores(data, bin_percentage)
    stats = calculate_stats(escores, binned_escores, bin_edges)
    
    metrics = []
    weights = []
    for i, (bin_start, bin_end, avg, std) in enumerate(stats):
        metric = avg / std
        metrics.append(metric)
        
    for i, (bin_start, bin_end, avg, std) in enumerate(stats):
        weight = metrics[i] / sum(metrics)
        weights.append(weight)
        #print(f"Bin {i+1}: {bin_start:.4f} to {bin_end:.4f} - Average: {avg:.4f}, Std Dev: {std:.4f}, Weight: {weight:.4f}")
    
    file_pattern = 'ETV-5/bin_*/bin_*.fasta_8mers/meme.html'
    files = sorted(glob.glob(file_pattern, recursive=True), key=lambda x: int(re.search(r'bin_(\d+)', x).group(1)))
    
    pwms = read_html_pwm(files)
    weighted_pwms = apply_weights(pwms, weights)
    
    sum_pwm = sum_pwms(weighted_pwms)
    
    bases = ['A', 'C', 'G', 'T']
            
    for position in sum_pwm:
        position_sum = sum(position)  
        print(f"{position_sum}")
    
    n_positions = sum_pwm.shape[0]
    n_bases = sum_pwm.shape[1]
    bar_width = 0.2 
    x = np.arange(n_positions)  

    fig, ax = plt.subplots()

    for i in range(n_bases):
        ax.bar(x + i * bar_width, sum_pwm[:, i], width=bar_width, label=bases[i])

    ax.set_xlabel('Position')
    ax.set_ylabel('Values')
    ax.set_title('Summed PWM Values by Position and Base')
    ax.set_xticks(x + bar_width * (n_bases - 1) / 2) 
    ax.set_xticklabels(range(1, n_positions + 1))  
    ax.legend()

    plt.show()

if __name__ == "__main__":
    main()
    
                   