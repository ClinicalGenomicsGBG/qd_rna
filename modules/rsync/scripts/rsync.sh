#!/bin/bash
set -e -o pipefail

echo "Copying files from ${SRC} to ${DST}"

rsync -a ${SRC} ${DST}
