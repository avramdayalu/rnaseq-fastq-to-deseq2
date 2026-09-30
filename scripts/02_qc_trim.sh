#!/usr/bin/env bash
# QC raw reads (FastQC), adapter/quality trim (fastp), summarise (MultiQC)
set -euo pipefail
mkdir -p results/fastqc_raw results/fastp data/trimmed
fastqc -t 2 -o results/fastqc_raw data/raw/*.fastq.gz
tail -n +2 samples.tsv | cut -f1 | while read -r run; do
  echo "== trimming $run =="
  fastp -i data/raw/${run}_1.fastq.gz -I data/raw/${run}_2.fastq.gz \
        -o data/trimmed/${run}_1.fastq.gz -O data/trimmed/${run}_2.fastq.gz \
        -w 2 -j results/fastp/${run}.json -h results/fastp/${run}.html
done
multiqc results -o results/multiqc
