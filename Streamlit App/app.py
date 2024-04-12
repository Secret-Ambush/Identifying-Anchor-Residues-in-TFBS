import streamlit as st
import subprocess
import tempfile
import shutil
from bs4 import BeautifulSoup
import json
import re
from collections import Counter
from pathlib import Path
from Bio import SeqIO
from io import StringIO
import glob
from Bio.SeqRecord import SeqRecord

def calculate_consensus(pwm):
    consensus = ''
    residues = ['A', 'C', 'G', 'T']
    
    for position in pwm:
        max_index = position.index(max(position))
        consensus += residues[max_index]
    
    return consensus

def split_into_bins(sequences, bin_size_percentage=50, overlap_percentage=5):
    total_sequences = len(sequences)
    
    sequences_per_bin = max(int(bin_size_percentage / 100 * total_sequences), 1)
    overlap_amount = max(int(overlap_percentage / 100 * sequences_per_bin), 1)
    
    binned_sequences = []
    seq_count = 0
    i = 0
    bin_count = 0
    while i < total_sequences:
        bin_sequences = sequences[i:min(i + sequences_per_bin, total_sequences)]
        bin_count += 1
        # Reset the sequence ID to ensure uniqueness within each bin
        for seq_record in bin_sequences:
            seq_count += 1
            seq_record.id = f"Bin_{bin_count}_Sequence_{seq_count}"
            seq_record.description = ""
        
        seq_count = 0
        
        # Instead of writing to a file, add the bin_sequences to a list
        binned_sequences.append((f"Bin_{bin_count}", bin_sequences))
        
        i += sequences_per_bin - overlap_amount
    
    return binned_sequences

def run_meme_analysis_on_bins(binned_sequences, tmpdir_path):
    for bin_name, bin_seqs in binned_sequences:
        fasta_file_path = tmpdir_path / f"{bin_name}.fasta"
        
        # Save the sequences to a FASTA file
        with fasta_file_path.open('w') as fasta_file:
            for seq in bin_seqs:
                fasta_file.write(f">{seq.id}\n{seq.seq}\n")
        
        # Prepare the output directory for the MEME analysis
        output_dir_path = tmpdir_path / f"{bin_name}_memeout"
        if output_dir_path.exists():
            shutil.rmtree(output_dir_path)  # Remove the existing directory and its contents
        output_dir_path.mkdir(parents=True, exist_ok=True)
        
        # Execute MEME analysis using Docker
        print(f"Running MEME analysis on {fasta_file_path}")
        try:
            subprocess.run([
                "docker", "run", "--rm",
                "-v", f"{tmpdir_path}:/data",
                "memesuite/memesuite:latest",
                "meme", f"/data/{fasta_file_path.name}",
                "-dna",
                "-o", f"/data/{output_dir_path.name}",
                "-nostatus",
                "-maxw", "8",
                "-minw", "4",
                "-nmotifs", "1",
                "-mod", "zoops",
                "-objfun", "classic",
                "-revcomp",
                "-markov_order", "0"
            ], check=True)
            print(f"MEME analysis completed for {bin_name}")
        except subprocess.CalledProcessError as e:
            print(f"MEME analysis failed for {bin_name}: {e}")
            
    process_meme_results(tmpdir_path)
 
def process_meme_results(tmpdir_path):
    file_pattern = f'{tmpdir_path}/Bin_*_memeout/meme.html'
    
    files = sorted(glob.glob(file_pattern, recursive=True), key=lambda x: int(re.search(r'bin_(\d+)', x).group(1)))

    pwm_sections = []
    consensus_sequences = []

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
                    pwm_section = data['motifs'][0].get('pwm', [])
                    pwm_sections.append(pwm_section)
                else:
                    st.write(f"No 'motifs' data found in {file_path}.")
            else:
                st.write(f"JSON data not found in the script tag of {file_path}.")
        else:
            st.write(f"Script tag containing 'var data' was not found in {file_path}.")

    # Prepare the content to be written in the downloadable file
    output_content = StringIO()
    output_content.write(f"Binning of ENTIRE dataset with {overlap_percentage}% overlap \n\n")
    output_content.write("Consensus sequences of bins: \n")
    
    for i, pwm_section in enumerate(pwm_sections, start=1):
        consensus = calculate_consensus(pwm_section)
        consensus_sequences.append(consensus)
        output_content.write(f"Bin {i}: {consensus}\n")
    
    output_content.write("\n\nAnalysis\n")
    for position in range(len(consensus_sequences[0])):
        residues_at_position = [seq[position] for seq in consensus_sequences]
        count = Counter(residues_at_position)
        
        most_common_residue, most_common_count = count.most_common(1)[0]
        
        output_content.write(f"Position {position+1}:\n")
        for residue, residue_count in count.items():
            output_content.write(f"  {residue}: {residue_count} times\n")
        output_content.write(f"Most common: {most_common_residue} appears {most_common_count} times\n\n")
    
    output_content.seek(0)
    
    # Create a download button in Streamlit for the analysis file
    st.download_button(
        label="Download Analysis Results",
        data=output_content,
        file_name="consensus_analysis.txt",
        mime="text/plain"
    )
                
   
if 'sequences' not in st.session_state:
    st.session_state.sequences = []

st.title('Protein Anchor Residue Identifier')

# Placeholder for your GIF link
# st.markdown("![Protein Animation](https://link_to_your_gif.gif)")

st.markdown("""
This app analyses protein enrichment scores to identify possible anchor residues, 
which are crucial for protein functionality and interactions. 
Upload your text file with the scores, and use the slider to adjust the threshold 
for identifying anchor residues. The table below will update with the potential 
anchor residues based on your criteria.
""")

# User inputs for binning preferences
with st.form("binning_preferences_form"):
    st.subheader("Binning Preferences")
    st.caption("Here you can decide how much percentage of the dataset should be included in one bin")
    bin_size_percentage = st.number_input("Bin Size Percentage", min_value=1, max_value=100, value=50)
    overlap_choice = st.radio("Bin Overlap", ('Discrete', 'Overlap'))
    if overlap_choice == 'Overlap':
        overlap_percentage = st.number_input("Overlap Percentage", min_value=0, max_value=100, value=5)
    else:
        overlap_percentage = 0
    st.subheader("Upload the top-enrichment dataset")
    uploaded_file = st.file_uploader(" ", type=["txt", "csv"])
    submit_button = st.form_submit_button("Confirm")

# Process file only if submit button is pressed
if submit_button and uploaded_file is not None:
    sequences = [SeqRecord(seq, id=f"Sequence_{i+1}", description="") for i, seq in enumerate(uploaded_file.getvalue().decode("utf-8").split('\n')) if seq.strip()]
    st.session_state.sequences = sequences  # Store sequences in session state

# Display and process sequences only if they are in the session state
if st.session_state.sequences:
    st.markdown(f"Chosen Bin Size: `{bin_size_percentage}%` Overlap Percentage: `{overlap_percentage}%`")
    fasta_io = StringIO()
    SeqIO.write(st.session_state.sequences, fasta_io, "fasta")
    fasta_io.seek(0)
    
    with tempfile.TemporaryDirectory() as tmpdir_path:
        tmpdir_path = Path(tmpdir_path)
        binned_sequences = split_into_bins(st.session_state.sequences, bin_size_percentage, overlap_percentage)
        
        st.write("Download Bins here!")
        for bin_name, bin_seqs in binned_sequences:
            fasta_io = StringIO()
            SeqIO.write(bin_seqs, fasta_io, "fasta")
            fasta_content = fasta_io.getvalue()
            st.download_button(f"Download {bin_name}", data=fasta_content, file_name=f"{bin_name}.fasta", mime="text/plain")

        col1, col2, col3 = st.columns([1,1,1])
        confirm = col2.button("Do you want to continue with MEME Analysis?")
        if confirm:
            run_meme_analysis_on_bins(binned_sequences, tmpdir_path)

    
    