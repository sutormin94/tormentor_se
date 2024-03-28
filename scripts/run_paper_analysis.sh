#!/usr/bin/bash 

mkdir -p paper/data/results

for dataset in streptococcus hmp; do
    while read acc; do 
        echo $dataset, $acc
        tormentor \
        --reads paper/data/raw/$acc\_1.fastq paper/data/raw/$acc\_2.fastq \
        --output paper/data/results/$acc \
        --threads 4 \
        --data-directory data/ \
        --minimum-self-pairing-percent 0.5
    done < paper/data/accessions_$dataset.txt
done