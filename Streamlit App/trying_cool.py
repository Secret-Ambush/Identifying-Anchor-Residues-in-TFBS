import uuid
from bs4 import BeautifulSoup
import streamlit as st
import subprocess
import tempfile
import shutil
import json
import re
from io import StringIO
from collections import Counter
from pathlib import Path
from Bio import SeqIO
import glob
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from datetime import datetime

min_motif = 0
max_motif = 0

def calculate_consensus(pwm):
    consensus = ''
    residues = ['A', 'C', 'G', 'T']
    
    for position in pwm:
        max_index = position.index(max(position))
        consensus += residues[max_index]
    
    return consensus

def split_into_bins(sequences, bin_size_percentage, overlap_percentage):
    total_sequences = len(sequences)
    
    sequences_per_bin = max(int(bin_size_percentage / 100 * total_sequences), 1)
    overlap_amount = max(int(overlap_percentage / 100 * sequences_per_bin), 0)
    
    binned_sequences = []
    seq_count = 0
    i = 0
    bin_count = 0
    while i < total_sequences:
        bin_sequences = sequences[i:min(i + sequences_per_bin, total_sequences)]
        bin_count += 1
        
        # Process each sequence in the bin
        for seq_record in bin_sequences:
            seq_count += 1
            seq_record.id = f"Sequence_{seq_count}"
            seq_record.description = ""
        
        binned_sequences.append((f"Bin_{bin_count}", bin_sequences))
        seq_count = 0
        i += sequences_per_bin - overlap_amount
    
    return binned_sequences

def run_meme_analysis_on_bins(binned_sequences, base_dir):
    base_dir = Path(base_dir).absolute()
    
    try:
        base_dir.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        print(f"Failed to create base directory {base_dir}: {e}")
        return
    
    with st.spinner('Please wait, accessing MEME tool'):

        for bin_name, bin_seqs in binned_sequences:
            input_dir = base_dir / f"{bin_name}"
            input_dir.mkdir(parents=True, exist_ok=True)  # Ensure directory is created

            input_basename = f"{bin_name}.fasta"
            fasta_file_path = input_dir / input_basename
            
            try:
                with open(fasta_file_path, 'w') as fasta_file:
                    for seq in bin_seqs:
                        fasta_file.write(f">{seq.id}\n{seq.seq}\n")
            except Exception as e:
                print(f"Failed to write to file {fasta_file_path}: {e}")
                continue
            
            print(f"Running Meme for {fasta_file_path}")

            command = [
                "docker", "run", "--rm",
                "-v", f"{input_dir}:/data",
                "memesuite/memesuite:latest",
                "meme", f"/data/{input_basename}",
                "-dna",
                "-o", f"/data/{input_basename}_8mers",
                "-nostatus",
                "-maxw", str(max_motif),
                "-minw", str(min_motif),
                "-nmotifs", "1",
                "-mod", "zoops",
                "-objfun", "classic",
                "-revcomp",
                "-markov_order", "0"
            ]

            try:
                subprocess.run(command, check=True)
                print(f"Done Running Meme for {fasta_file_path}")
                mymsg2.info(f"Done Running Meme for {input_basename}")
                                
            except Exception as e:
                print(f"Failed to run command for bin {bin_name}: {e}")           
    st.success('Done!', icon = "✅")
    mymsg2.write("")
 
def process_meme_results(tmpdir_path):
    print("*********************************************\n")
    file_pattern = f'{tmpdir_path}/Bin_*/Bin_*.fasta_8mers/meme.html'
    files = sorted(glob.glob(file_pattern, recursive=True), key=lambda x: int(re.search(r'/Bin_(\d+)/', x).group(1)))
    print("Found HTML files:", files)
    
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

    print("*********************************************\n")
    
    # Preparing downloadable file
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

    output_bytes = output_content.getvalue().encode()
    
    col1, col2, col3 = st.columns([1,1,1])
    col2.download_button(
        label="Download Analysis Results",
        data=output_bytes,
        file_name="consensus_analysis.txt",
        mime="text/plain"
    )

if 'sequences' not in st.session_state:
    st.session_state.sequences = []

st.title('Protein Anchor Residue Identifier')
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
    overlap_percentage = st.number_input("Overlap Percentage", min_value=0, max_value=100, value=0)
    st.caption("An overlap of 0 between bins means Discrete Binning")
    st.subheader("Motif Generation Details")
    min_motif = st.number_input("Minimum Motif Width", min_value=1, max_value=8, value=4)
    max_motif = st.number_input("Maximum Motif Width", min_value=1, max_value=8, value=4)
    st.caption("Keeping them same ensures a single length of motif is generated for each bin")
    st.subheader("Upload the top-enrichment dataset")
    uploaded_file = st.file_uploader(" ", type=["txt", "csv"])
    submit_button = st.form_submit_button("Confirm")

if submit_button and uploaded_file is not None:
    uploaded_content = uploaded_file.getvalue().decode("utf-8").strip()
    lines = uploaded_content.split('\n')

    sequences = []
    for i, line in enumerate(lines[1:]):  
        if line.strip(): 
            columns = line.split('\t')
            if columns: 
                sequence = columns[0].replace('.', '')  # Remove '.' from sequences
                seq_record = SeqRecord(Seq(sequence), id=f"Sequence{i+1}", description="")
                sequences.append(seq_record)
        st.session_state.sequences = sequences 

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
        mymsg = col2.empty() 
        confirm = mymsg.button("Do you want to continue with MEME Analysis?")
        mymsg2 = st.empty()
        if confirm:
            mymsg.write("")
            run_meme_analysis_on_bins(binned_sequences, tmpdir_path)
            process_meme_results(tmpdir_path)
            