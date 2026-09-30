#!/usr/bin/env bash
# Download the first N read pairs per run from SRA (subset to fit Codespaces disk/RAM)
set -euo pipefail
N=2000000
mkdir -p data/raw
tail -n +2 samples.tsv | cut -f1 | while read -r run; do
  echo "== $run =="
  fastq-dump --split-files --gzip -X "$N" -O data/raw "$run"
done
