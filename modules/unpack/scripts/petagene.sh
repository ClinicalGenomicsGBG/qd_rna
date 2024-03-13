#!/bin/bash

set -e -o pipefail

echo -e "UNPACK_INIT:\n${UNPACK_INIT}\n"

eval ${UNPACK_INIT}

/bin/env

echo "CMD: petasuite \
--decompress \
--numthreads $THREADS \
--validate off \ 
$COMPRESSED_PATH"

petasuite \
    --decompress \
    --numthreads $THREADS \
    $COMPRESSED_PATH

eval ${UNPACK_EXIT}
