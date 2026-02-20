#!/usr/bin/env bash
set -euo pipefail

module load samtools/1.16

bam="$1"
cram="$2"
ref="$3"
threads="${4:-4}"

echo "[cram_compress] bam=$bam"
echo "[cram_compress] cram=$cram"
echo "[cram_compress] ref=$ref"
echo "[cram_compress] threads=$threads"

samtools view -@ "$threads" -C -T "$ref" -o "$cram" "$bam"
samtools index -@ "$threads" "$cram"

