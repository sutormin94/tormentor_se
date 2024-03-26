from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import subprocess
import os

def run_prodigal(fasta_file, output_file):
    return subprocess.call(f'prodigal -i {fasta_file} -p meta -o {output_file} -f gbk')

def correct_contig(genbank_file):
    # Parse the GenBank file
    record = next(SeqIO.parse(genbank_file, 'genbank'))

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

    return record
    


def run_cmscan(fasta_file, output_file):
    return subprocess.call(f'prodigal -i {fasta_file} -p meta -o {output_file} -f gbk')