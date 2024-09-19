#!/bin/bash

set -e -o pipefail

eval "${UNPACK_INIT}"

echo "spring \
--decompress \
--gzipped-fastq \
--working-dir $TMPDIR \
--output-file $EXTRACTED_PATH \
--input-file $COMPRESSED_PATH \
--num-threads $THREADS"

spring \
    --decompress \
    --gzipped-fastq \
    --working-dir $TMPDIR \
    --output-file $EXTRACTED_PATH \
    --input-file $COMPRESSED_PATH \
    --num-threads $THREADS

eval "${UNPACK_EXIT}"