#!/bin/bash

set -e -o pipefail

echo -e "UNPACK_INIT:\n${UNPACK_INIT}\n"

eval ${UNPACK_INIT}

/bin/env

echo "CMD: petasuite \
--decompress \
--numthreads $THREADS \
--validate off \
--dstpath $(dirname $EXTRACTED_PATH) \
$COMPRESSED_PATH"

petasuite \
    --decompress \
    --numthreads $THREADS \
    --dstpath $(dirname $EXTRACTED_PATH) \
    $COMPRESSED_PATH

eval ${UNPACK_EXIT}
