from tormentor.assembly import run_spades, filter_and_rename_spades_transcripts
#from tormentor.vnom import run_vnom
from argparse import ArgumentParser
import os

def main():

    argument_parser = ArgumentParser(description="An obelisk prediction and annotation pipeline from stranded RNA-Seq data")
    argument_parser.add_argument('--reads', nargs=2, metavar='<FASTQ file>', help='forward and reverse FASTQ files from stranded RNA-Seq')
    argument_parser.add_argument('-o', '--output', help='output directory')
    argument_parser.add_argument('-s', '--minimum-self-pairing-percent', help='minimum percent in secondary structure')
    argument_parser.add_argument('--threads', help='number of CPU threads to use', default=4)
    argument_parser.add_argument('--stranded-type', choices=['fr', 'rf'], default='fr', help='stranded library layout (for rnaSPAdes)')
    arguments = argument_parser.parse_args()

    spades_directory              = os.path.join(arguments.output, 'step_1')
    spades_transcripts            = os.path.join(spades_directory, 'transcripts.fasta')
    spades_transcripts_clear      = os.path.join(spades_directory, 'transcripts-clear.fasta')
    vnom_directory                = os.path.join(arguments.output, 'step_2')
    secondary_structure_directory = os.path.join(arguments.output, 'step_3')
    annotation_directory          = os.path.join(arguments.output, 'step_4')

    # run assembly

    step_1_return_code = run_spades(arguments.reads[0], arguments.reads[1], arguments.threads, spades_directory, stranded_type=arguments.stranded_type) == 0
    step_1_return_code = 0
    if step_1_return_code != 0:
        print('Error while assembling RNA sequences')
        exit(1)
    
    filter_and_rename_spades_transcripts(spades_transcripts, spades_transcripts_clear)

    # run circRNA detection

    run_vnom(spades_transcripts_clear)

    # run secondary structure analysis
        
    step_2_return_code = 0
    
    # run annotation

    # 

if __name__ == '__main__':
    main()