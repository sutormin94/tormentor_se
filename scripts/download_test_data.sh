#!/usr/bin/env bash

cd tests
fasterq-dump SRR5949245
mv SRR5949245_1.fastq reads_1.fastq
mv SRR5949245_2.fastq reads_2.fastq
rm SRR5949245.fastq