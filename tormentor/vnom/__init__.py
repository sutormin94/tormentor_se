from Bio import SeqIO
import subprocess
import sys
import os

def run_vnom(fasta_preffix, max_length=2000, stdout=sys.stdout, stderr=sys.stderr):
    if os.path.isdir("0_non_singleton_clusters"):
        os.system('rm -r 0_non_singleton_clusters')
    command = (
        'python vnom/VNom.py '
        f'-i {os.getcwd()}/{fasta_preffix} '
        f'-max {max_length} '
        '-CF_k 10 '
        '-CF_simple 0 '
        '-CF_tandem 1 '
        '-USG_vs_all 1 '
    )
    return subprocess.call(command, shell=True, stdout=stdout, stderr=stderr)

def split_vnom_candidates(fasta_file, output_directory, min_length, max_length):
    files = []
    for r, record in enumerate(SeqIO.parse(fasta_file, 'fasta')):
        file = os.path.join(output_directory, f'{r+1}.fasta')
        record.id          = f'obelisk-candidate:{r+1}'
        record.description = f'obelisk-candidate:{r+1}'
        if min_length <= len(record.seq) <= max_length:
            SeqIO.write([record], open(file, 'w'), 'fasta')
        files.append(file)
    return files