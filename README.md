# Tormentor

An [obelisk](https://www.biorxiv.org/content/10.1101/2024.01.20.576352v1.full.pdf) prediction and annotation pipeline

![](assets/tormentor.png)

**WARNING: this tool is still under development!**

## Setup

Configuring Tormentor requires `conda` (anaconda/miniconda/miniforge). After installation 
an environment named `tormentor` will be created.

```bash
(base) $ git clone git@github.com:/omixlab/tormentor
(base) $ cd tormentor
(base) $ make setup
(base) $ conda activate tormentor
```

## Running

```bash
(tormentor) $ tormentor --reads reads_1.fastq reads_2.fastq --output results/ --threads 4
```

## Results

Results are provided in FASTA (nucleotide and protein), GenBank and CSV format, along with
intermediate results from each program used in the pipeline.

## The name

The name is a reference to *Obelisk the Tormentor*, a card from the trading card game / anime / mangá "Yu-Gi-Oh!". 

# Cite us

Kremer, F (2024). *Tormentor: An obelisk prediction and annotation pipeline*.
