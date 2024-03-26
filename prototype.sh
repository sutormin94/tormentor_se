# create directories

mkdir -p tests/results/

# fastp

#fastp \
#    --in1=tests/reads_1.fastq --in2=tests/reads_2.fastq \
#    --out1=tests/results/reads_1.fastq --out2=tests/results/reads_2.fastq \
#    --average_qual=30 --n_base_limit=0 --cut_front --cut_tail

# spades

#rnaspades.py  \
#    -1 tests/results/reads_1.fastq \
#    -2 tests/results/reads_2.fastq \
#    -o tests/results/step_1/

# vnom

mkdir -p tests/results/step_2/
cp -p $PWD/tests/results/step_1/transcripts.fasta $PWD/tests/results/step_2/transcripts.fasta
python vnom/VNom.py \
    -i $PWD/tests/results/step_2/transcripts \
    -max 2000 \
    -CF_k 10 \
    -CF_simple 0 \
    -CF_tandem 1 \
    -USG_vs_all 1

# prodigal

mkdir -p tests/results/step_3/
prodigal -i tests/results/step_2/transcripts_cir.fasta -p meta -f gbk -o tests/results/step_3/annotation.gbk