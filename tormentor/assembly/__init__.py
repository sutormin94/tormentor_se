from Bio import SeqIO
import subprocess

def run_spades(forward_reads, reverse_reads, threads, output_directory, stranded_type):
    commandline = (
        'rnaspades.py '
        f'-1 {forward_reads} -2 {reverse_reads} '
        f'-t {threads} '
        f'--ss-{stranded_type} '
        f'-o {output_directory} '
    )
    return subprocess.call(commandline, shell=True)

def filter_and_rename_spades_transcripts(input_fasta, output_fasta):

    transcripts = SeqIO.parse(input_fasta, 'fasta')
    with open(output_fasta, 'w') as writer:
        for t, transcript in enumerate(transcripts):
            keep = True
            for base in str(transcript.seq):
                if base not in 'atcgATCG':
                    keep = False
            if not keep:
                continue    
            SeqIO.write([transcript], writer, 'fasta')