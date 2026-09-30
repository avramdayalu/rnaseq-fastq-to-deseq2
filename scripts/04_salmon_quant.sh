#!/usr/bin/env bash
# Quantify trimmed paired-end reads against the protein-coding Salmon index
set -euo pipefail
mkdir -p results/salmon
tail -n +2 samples.tsv | cut -f1 | while read -r run; do
  echo "== quantifying $run =="
  salmon quant -i data/ref/salmon_index_pc -l A \
    -1 data/trimmed/${run}_1.fastq.gz -2 data/trimmed/${run}_2.fastq.gz \
    -p 2 --gcBias -o results/salmon/${run}
done
multiqc results -o results/multiqc -f
