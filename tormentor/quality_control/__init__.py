import subprocess
import sys
import os

def run_fastp(input_reads, trimmed_reads, min_quality=30, stdout=sys.stdout, stderr=sys.stderr):
    if (len(input_reads)==2) and (len(trimmed_reads)==2):
        forward_reads=input_reads[0]
        reverse_reads=input_reads[1]
        forward_reads_trim=trimmed_reads[0]
        reverse_reads_trim=trimmed_reads[1]
        
        commandline = (
            'fastp '
            f'--in1={forward_reads} --in2={reverse_reads} '
            f'--out1={forward_reads_trim} --out2={reverse_reads_trim} '
            f'--average_qual={min_quality} --n_base_limit=0 --cut_front --cut_tail '
        )
        
    elif (len(input_reads)==1) and (len(trimmed_reads)==1):
        forward_reads=input_reads[0]
        forward_reads_trim=trimmed_reads[0]
        
        commandline = (
            'fastp '
            f'--in1={forward_reads} '
            f'--out1={forward_reads_trim} '
            f'--average_qual={min_quality} --n_base_limit=0 --cut_front --cut_tail '
        )        
         
    return subprocess.call(commandline, shell=True, stdout=stdout, stderr=stderr)