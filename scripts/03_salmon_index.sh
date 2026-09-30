#!/usr/bin/env bash
# Build a Salmon index on GENCODE human protein-coding transcripts (fits 8 GB RAM)
set -euo pipefail
BASE=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release
mkdir -p data/ref
FA=$(curl -s "$BASE/" | grep -o 'gencode\.v[0-9]*\.pc_transcripts\.fa\.gz' | head -1)
echo "Using $FA" | tee data/ref/GENCODE_VERSION.txt
wget -q -O data/ref/$FA "$BASE/$FA"

# transcript -> gene table for tximport (from FASTA headers)
zcat data/ref/$FA | grep '^>' | sed 's/^>//' | \
  awk -F'|' 'BEGIN{OFS="\t"; print "TXNAME","GENEID","SYMBOL"} {print $1,$2,$6}' \
  > data/ref/tx2gene.tsv

salmon index -t data/ref/$FA -i data/ref/salmon_index_pc --gencode -k 31 -p 2
