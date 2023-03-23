#! /bin/bash

set -e

module load samtools/1.16

samtools view \
    -s $(bc <<< "1337 + ${_QLUCORE_SUBSAMPLE_FRAC}") \
    --threads ${_QLUCORE_SUBSAMPLE_THREADS} \
    --bam --output ${_QLUCORE_SUBSAMPLE_OUTPUT_BAM} \
    ${_QLUCORE_SUBSAMPLE_INPUT_BAM}