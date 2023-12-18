#! /bin/bash

set -e

eval ${_QLUCORE_SUBSAMPLE_INIT}

input_base=$(basename ${_QLUCORE_SUBSAMPLE_INPUT_BAM} .bam)
input_root=$(dirname ${_QLUCORE_SUBSAMPLE_INPUT_BAM})
input_temp=${TMPDIR}/${input_base}.bam
output_temp=${TMPDIR}/${input_base}.downsampled.bam

cp ${_QLUCORE_SUBSAMPLE_INPUT_BAM} ${input_temp}

samtools view \
    -s $(bc <<< "1337 + ${_QLUCORE_SUBSAMPLE_FRAC}") \
    --threads ${_QLUCORE_SUBSAMPLE_THREADS} \
    --bam --output ${output_temp} \
    ${input_temp}

mv ${output_temp} ${input_root}/