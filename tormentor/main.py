from tormentor.quality_control import run_fastp
from tormentor.assembly import run_spades, filter_and_rename_spades_transcripts
from tormentor.vnom import run_vnom, split_vnom_candidates
from tormentor.annotation import run_prodigal, identify_cds, correct_contig, run_cmscan, add_rfam_sites_to_record
from tormentor.secondary_structure import run_rnafold, run_rnaplot, compute_self_pairing_percent
from tormentor.logo import logo
from argparse import ArgumentParser
from Bio import SeqIO
from BCBio import GFF
import multiprocessing
import os

CPU_COUNT = multiprocessing.cpu_count()

def main():

    print(logo)

    argument_parser = ArgumentParser(description="An obelisk prediction and annotation pipeline from stranded RNA-Seq data")
    argument_parser.add_argument('--reads', nargs='+', metavar='<FASTQ file>', help='forward (and reverse) FASTQ file(s) from stranded RNA-Seq')
    argument_parser.add_argument('-o', '--output', help='output directory')
    argument_parser.add_argument('-s', '--minimum-self-pairing-percent', help='minimum percent in secondary structure', default=0.5, type=float)
    argument_parser.add_argument('--threads', help='number of CPU threads to use', default=CPU_COUNT)
    argument_parser.add_argument('--min-qual', help='minimal quality of reads in Phred score', default=30)
    argument_parser.add_argument('--min-identity', help='minimal indetity for oblin identification using BLAST', default=0.7)
    argument_parser.add_argument('--evalue', help='e-value threshold used by BLAST and cmscan', default=1e-6)
    argument_parser.add_argument('--stranded-type', choices=['fr', 'rf'], default=None, help='stranded library layout (for rnaSPAdes)')
    argument_parser.add_argument('--min-len', default=900, type=int, help='min length for predicted obelisks')
    argument_parser.add_argument('--max-len', default=2000, type=int, help='max length for predicted obelisks')
    argument_parser.add_argument('--data-directory', help='directory containing the dataset files', required=True)
    arguments = argument_parser.parse_args()

    fastp_directory               = os.path.join(arguments.output, 'step_1')
    spades_directory              = os.path.join(arguments.output, 'step_2')
    spades_transcripts            = os.path.join(spades_directory, 'transcripts.fasta')
    spades_transcripts_clear      = os.path.join(spades_directory, 'transcripts-clear.fasta')
    vnom_directory                = os.path.join(arguments.output, 'step_3')
    vnom_directory_candidates     = os.path.join(arguments.output, 'step_3', 'candidates')
    prodigal_directory            = os.path.join(arguments.output, 'step_4')
    prodigal_directory_corrected  = os.path.join(arguments.output, 'step_4', 'corrected')

    os.system(f'mkdir -p {arguments.output}/logs/')
    step_1_log_handler = open(f'{arguments.output}/logs/step_1.log', 'w')
    step_2_log_handler = open(f'{arguments.output}/logs/step_2.log', 'w')
    step_3_log_handler = open(f'{arguments.output}/logs/step_3.log', 'w')
    step_4_log_handler = open(f'{arguments.output}/logs/step_4.log', 'w')
    
    if len(reads)==2:
        reads_1_raw  = arguments.reads[0]
        reads_2_raw  = arguments.reads[1]
        reads_1_trim = os.path.join(fastp_directory, 'reads_1.fastq')
        reads_2_trim = os.path.join(fastp_directory, 'reads_2.fastq')
        input_reads=[reads_1_raw, reads_2_raw]
        trimmed_reads=[reads_1_trim, reads_2_trim]
        
    elif len(reads)==1:
        reads_1_raw  = arguments.reads[0]
        reads_1_trim = os.path.join(fastp_directory, 'reads_1.fastq')
        input_reads=[reads_1_raw]
        trimmed_reads=[reads_1_trim]
        
    vnom_input   = os.path.join(vnom_directory, 'transcripts')
    vnom_output  = os.path.join(vnom_directory, 'transcripts_cir.fasta')
    
    # step_1 quality control
    
    print('Step 1: Running quality control ...')
    
    os.system(f'mkdir -p {fastp_directory}')
    step_1_return_code = run_fastp(
        input_reads, 
        trimmed_reads, 
        min_quality=arguments.min_qual,
        stdout=step_1_log_handler,
        stderr=step_1_log_handler
    )
    step_1_log_handler.close()
    if step_1_return_code != 0:
        print('Error while read quality control')
        exit(1)

    # step_2: run assembly

    print('Step 2: Running de novo meta-transcriptome assembly ...')
    
    os.system(f'mkdir -p {spades_directory}')
    step_2_return_code = run_spades(
        trimmed_reads, 
        arguments.threads, 
        spades_directory, 
        stranded_type=arguments.stranded_type,
        stdout=step_2_log_handler,
        stderr=step_2_log_handler
    )
    step_2_log_handler.close()
    
    if step_2_return_code != 0:
        print('Error while assembling RNA sequences')
        exit(1)
    
    filter_and_rename_spades_transcripts(spades_transcripts, spades_transcripts_clear)

    
    # step_3: run viroid circRNA detection

    print('Step 3: Running viroid-like circRNA prediction ...')

    os.system(f'mkdir -p {vnom_directory}')
    os.system(f'mkdir -p {vnom_directory_candidates}')
    os.system(f'cp {spades_transcripts_clear} {vnom_input}.fasta')
    
    step_3_return_code = run_vnom(vnom_input, max_length=arguments.max_len, stdout=step_3_log_handler, stderr=step_3_log_handler)
    step_3_log_handler.close()
    if step_3_return_code != 0:
        print('Error while identifying viroid-like sequences')
        exit(1)

    obelisks_candidates = split_vnom_candidates(vnom_output, vnom_directory_candidates, min_length=arguments.min_len, max_length=arguments.max_len)
    
    # step_4: run annotation: prodigal

    print('Step 4: Running annotation and secondary structure analysis ...')

    os.system(f'mkdir -p {prodigal_directory}')
    os.system(f'mkdir -p {prodigal_directory_corrected}')

    final_obelisk_count = 0

    for obelisk_id, obelisk_candidate in enumerate(obelisks_candidates):
        if not os.path.isfile(obelisk_candidate):
            continue
        prodigal_output = os.path.join(prodigal_directory, os.path.basename(obelisk_candidate) + '.gff')
        step_4_return_code = run_prodigal(obelisk_candidate, prodigal_output, stdout=step_4_log_handler, stderr=step_4_log_handler)

        if step_4_return_code != 0:
            print('Error while identifying CDSs in the sequences')
            exit(1)
        try:
            obelisk_record = next(GFF.parse(open(prodigal_output)) )
            obelisk_record.seq = next(SeqIO.parse(obelisk_candidate, 'fasta')).seq
        except:
            continue
        obelisk_candidate_corrected_fasta = os.path.join(prodigal_directory_corrected, os.path.basename(obelisk_candidate))
        obelisk_candidate_corrected_rnafold = os.path.join(prodigal_directory_corrected, os.path.basename(obelisk_candidate)) + '.txt'
        obelisk_candidate_corrected_genbank = os.path.join(prodigal_directory_corrected, os.path.basename(obelisk_candidate) + '.gbk')           
        record = correct_contig(obelisk_record, obelisk_candidate_corrected_genbank)

        record.id = 'obelisk'
        record.description = 'obelisk'
        SeqIO.write([record], obelisk_candidate_corrected_fasta, 'fasta')
        SeqIO.write([record], obelisk_candidate_corrected_genbank, 'genbank')

        if record is None:
            continue

        record = identify_cds(record, arguments.data_directory + '/proteins/oblins', prodigal_directory, evalue_treshold=arguments.evalue, min_identity=arguments.min_identity)

        oblin_count = 0

        for feature in record.features:
            if feature.type == 'CDS' and 'oblin' in feature.qualifiers.get('product',''):
                oblin_count += 1
        sites = run_cmscan(obelisk_candidate_corrected_fasta, prodigal_directory, arguments.data_directory + '/cms',  stdout=step_4_log_handler, stderr=step_4_log_handler, evalue_threshold=arguments.evalue)
        
        record = add_rfam_sites_to_record(record, sites)

        record.features = sorted(record.features, key=lambda x: int(x.location.start))

        SeqIO.write([record], obelisk_candidate_corrected_fasta, 'fasta')
        SeqIO.write([record], obelisk_candidate_corrected_genbank, 'genbank')

        # run secondary structure analysis
        
        run_rnafold(obelisk_candidate_corrected_fasta, obelisk_candidate_corrected_rnafold, stdout=step_4_log_handler, stderr=step_4_log_handler)
        run_rnaplot(obelisk_candidate_corrected_rnafold, stdout=step_4_log_handler, stderr=step_4_log_handler)

        self_pairing_percent = compute_self_pairing_percent(obelisk_candidate_corrected_rnafold)
        if self_pairing_percent >= arguments.minimum_self_pairing_percent and oblin_count > 0:
            os.system(f'cp {obelisk_candidate_corrected_fasta} {arguments.output}/obelisk_{obelisk_id+1}.fasta')
            os.system(f'cp {obelisk_candidate_corrected_fasta}.txt {arguments.output}/obelisk_{obelisk_id+1}.ss.txt')
            os.system(f'cp {obelisk_candidate_corrected_fasta}.txt.pdf {arguments.output}/obelisk_{obelisk_id+1}.ss.pdf')
            os.system(f'cp {obelisk_candidate_corrected_fasta}.txt.svg {arguments.output}/obelisk_{obelisk_id+1}.ss.svg')
            os.system(f'cp {obelisk_candidate_corrected_genbank} {arguments.output}/obelisk_{obelisk_id+1}.gbk')
            final_obelisk_count += 1
            print('obelisk: ', obelisk_id+1)
            print(' - self-pairing percent : ', self_pairing_percent)
            print(' - number of oblin genes: ', oblin_count)
    step_4_log_handler.close()
    print(f'Finished! {final_obelisk_count} obelisks identified!')
    
if __name__ == '__main__':
    main()
