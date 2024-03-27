import subprocess
import sys
import os
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPDF


def run_rnafold(fasta_file, output_file, stdout=sys.stdout, stderr=sys.stderr):
    return subprocess.call(f'RNAfold --noLP --circ {fasta_file} > {output_file}', shell=True, stdout=stdout, stderr=stderr)

def run_rnaplot(rnafold_file, stdout=sys.stdout, stderr=sys.stderr):
    rnafold_file = os.path.abspath(rnafold_file)
    rnafold_file_basename = os.path.basename(rnafold_file)
    cwd = os.getcwd()
    rnafold_file_dir = os.path.dirname(rnafold_file)
    os.chdir(rnafold_file_dir)
    return_code = subprocess.call(f'RNAplot -o svg {rnafold_file}', shell=True, stdout=stdout, stderr=stderr)
    os.system(f"mv obelisk_ss.svg {rnafold_file_basename}.svg")
    drawing = svg2rlg(f"{rnafold_file_basename}.svg")
    renderPDF.drawToFile(drawing, f"{rnafold_file_basename}.pdf")
    os.chdir(cwd)
    return return_code


def compute_self_pairing_percent(rnafold_file):
    lines               = open(rnafold_file).readlines()
    sequence            = lines[1].strip("\n")
    secondary_structure = lines[2].strip("\n")
    secondary_structure = secondary_structure[0:len(sequence)]
    return (secondary_structure.count('(') + secondary_structure.count(')')) / len(secondary_structure)