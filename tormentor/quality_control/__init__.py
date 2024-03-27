import subprocess
import sys
import os

def run_fastp(forward_reads, reverse_reads, forward_reads_trim, reverse_reads_trim, min_quality=30, stdout=sys.stdout, stderr=sys.stderr):
    commandline = (
        'fastp '
        f'--in1={forward_reads} --in2={reverse_reads} '
        f'--out1={forward_reads_trim} --out2={reverse_reads_trim} '
        f'--average_qual={min_quality} --n_base_limit=0 --cut_front --cut_tail '
    )
    return subprocess.call(commandline, shell=True, stdout=stdout, stderr=stderr)