from Bio import SeqIO
import subprocess
import sys
import os

def run_spades(trimmed_reads, threads, output_directory, stranded_type, stdout=sys.stdout, stderr=sys.stderr):
    
    # Check if assembly already exists.
    spades_transcripts_path=os.path.join(output_directory, 'transcripts.fasta')
    
    if os.path.isfile(spades_transcripts_path):
        
        print(f'File with assembled transcripts already exists, skipping rnaSPades: {spades_transcripts_path}', file=sys.stdout)
        
        return   
    
    else:
        
        if stranded_type:
            stranded_type_arg = f'--ss-{stranded_type} '
        else:
            stranded_type_arg = ''
            
        if len(trimmed_reads)==2:
            forward_reads_trim=trimmed_reads[0]
            reverse_reads_trim=trimmed_reads[1]        
        
            commandline = (
                'rnaspades.py '
                f'-1 {forward_reads_trim} -2 {reverse_reads_trim} '
                f'-t {threads} '
                f'{stranded_type_arg} '
                f'-o {output_directory} '
            )
        
        elif len(trimmed_reads)==1:
            forward_reads_trim=trimmed_reads[0]      
        
            commandline = (
                'rnaspades.py '
                f'-s {forward_reads_trim} '
                f'-t {threads} '
                f'{stranded_type_arg} '
                f'-o {output_directory} '
            )   
            
        return subprocess.call(commandline, shell=True, stdout=stdout, stderr=stderr)

def filter_and_rename_spades_transcripts(input_fasta, output_fasta):
    
    # Check if a file with re-named transcripts already exists.
    if os.path.isfile(output_fasta):
        
        print(f'File with re-named transcripts already exists, skip renaming: {output_fasta}', file=sys.stdout)    

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