#!/bin/bash

set -e -o pipefail

module load petasuite

echo "CMD: petasuite \
--decompress \
--numthreads $THREADS \
$COMPRESSED_PATH"

petasuite \
    --decompress \
    --numthreads $THREADS \
    $COMPRESSED_PATH
