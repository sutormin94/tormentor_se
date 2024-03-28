#!/usr/bin/bash 

mkdir -p paper/data/
cd paper/data/
mkdir -p raw

for dataset in streptococcus hmp; do
    while read acc; do 
        echo $dataset, $acc
        fasterq-dump $acc --split-3; 
        mv $acc\_1.fastq raw/
        mv $acc\_2.fastq raw/
        rm *.fastq
    done < accessions_$dataset.txt
done
