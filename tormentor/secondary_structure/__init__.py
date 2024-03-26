import subprocess

def run_rnafold(fasta_file, output_file):
    subprocess.call(f'RNAfold --circ {fasta_file} > {output_file}', shell=True)

def compute_self_pairing_percent(rnafold_file):
    secondary_structure = open(rnafold_file).readlines()[2].strip("\n")
    return (secondary_structure.count('(') + secondary_structure.count(')')) / len(secondary_structure)