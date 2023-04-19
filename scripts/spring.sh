#!/bin/bash

set -e -o pipefail

module load spring

echo "spring \
--decompress \
--gzipped-fastq \
--working-dir $TMPDIR \
--output-file $EXTRACT_PATH \
--input-file $COMPRESSED_PATH \
--num-threads $THREADS"

spring \
    --decompress \
    --gzipped-fastq \
    --working-dir $TMPDIR \
    --output-file $EXTRACT_PATH \
    --input-file $COMPRESSED_PATH \
    --num-threads $THREADS
