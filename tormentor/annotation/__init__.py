from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Blast import NCBIXML
from BCBio import GFF
import subprocess
import glob
import sys
import os
import re

def run_prodigal(fasta_file, output_file, stdout=sys.stdout, stderr=sys.stderr):

    return subprocess.call(f'prodigal -i {fasta_file} -p meta -o {output_file} -f gff', shell=True, stdout=stdout, stderr=stderr)

def correct_contig(record, output_file):

    # Find the largest feature on the forward strand
    largest_feature = None
    for feature in record.features:
        if feature.strand == 1:  # Check if feature is on the forward strand
            if largest_feature is None or \
                    (feature.location.end - feature.location.start + 1) > \
                    (largest_feature.location.end - largest_feature.location.start + 1):
                largest_feature = feature

    # Adjust positions based on the largest feature
    if largest_feature is not None:
        offset = largest_feature.location.start - 1
        new_features = []
        for feature in record.features:
            new_location = FeatureLocation(
                feature.location.start - offset,
                feature.location.end - offset,
                feature.location.strand)
            new_feature = SeqFeature(new_location, type=feature.type)
            new_features.append(new_feature)
        record.features = new_features
        record.seq = record.seq[offset:] + record.seq[:offset]
    record.annotations['molecule_type'] = 'DNA'
    SeqIO.write([record], output_file, 'genbank')
    return record

def identify_cds(record, protein_database, output_directory, min_identity=0.8, evalue_treshold=1e-10):
    for f, feature in enumerate(record.features):
        blast_output_file = os.path.join(output_directory, f'cds_{f}.xml')
        feature_sequence  = str(feature.location.extract(record.seq))
        feature_sequence_file = f'{output_directory}/{f}_query.fa'
        with open(feature_sequence_file, 'w') as writer:
            writer.write(f'>{f}\n{feature_sequence}\n')
        
        if len(feature_sequence)>0:
            
            print(f'Running blastx for {feature_sequence_file}')
            
            subprocess.run(
                f'blastx -query {feature_sequence_file} -db {protein_database} -outfmt 5 -out {blast_output_file}',
                shell=True
            )
            
            for blast_record in NCBIXML.parse(open(blast_output_file)):
                for alignment in blast_record.alignments:
                    identity = alignment.hsps[0].identities / alignment.hsps[0].align_length
                    evalue = alignment.hsps[0].expect
                    if identity >= min_identity and evalue < evalue_treshold:
                        feature.qualifiers['product'] = alignment.hit_def
                        break
                else:
                    feature.qualifiers['product'] = 'hypothetical protein'
                
    return record
    
def run_cmscan(fasta_file, output_directory, cm_directory, evalue_threshold=1e-6, stdout=sys.stdout, stderr=sys.stderr):
    sites = []
    for cm_file in glob.glob(f'{cm_directory}/*.cm'):
        cm_family = os.path.splitext(os.path.basename(cm_file))[0]
        subprocess.call(f'cmpress -F {cm_file}', shell=True, stdout=stdout, stderr=stderr)
        subprocess.call(f'cmscan --tblout {output_directory}/{cm_family}.txt {cm_file} {fasta_file}', shell=True, stdout=stdout, stderr=stderr)
        for line in open(f'{output_directory}/{cm_family}.txt'):
            if line.startswith('#'):
                continue
            site = {}
            site['family'] = line[0:20].strip(' ')
            site['evalue'] = float(line[138:147].strip(' '))
            site['start']  = int(line[85:93].strip(' '))
            site['end']    = int(line[94:102].strip(' '))
            site['strand'] = int(line[107] + '1')
            if site['evalue'] < evalue_threshold:
                sites.append(site)
    return sites
    
def add_rfam_sites_to_record(record, sites):
    for s, site in enumerate(sites):
        location = FeatureLocation(start=site['start'], end=site['end'], strand=site['strand'])
        feature  = SeqFeature(type='misc_feature', location=location, qualifiers={'note': f'motif:{site["family"]}'})
        record.features.append(feature)
    return record