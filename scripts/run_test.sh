#!/usr/bin/env bash

tormentor \
    --reads tests/reads_1.fastq tests/reads_2.fastq \
    -o tests/results/ \
    --data-directory data/