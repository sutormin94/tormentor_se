#!/usr/bin/bash 

mkdir -p paper/data/
cd paper/data/
mkdir -p raw

while read acc; do 
    prefetch $acc
    fasterq-dump $acc --split-3; 
    mv $acc\_1.fastq raw/
    mv $acc\_2.fastq raw/
done < accessions.txt

rm *.fastq

