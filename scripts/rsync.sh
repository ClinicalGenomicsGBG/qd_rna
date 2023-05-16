#!/bin/bash
set -e -o pipefail

SRC_ARRAY=( $SRC )
DST_ARRAY=( $DST )

for sidx in ${!SRC_ARRAY[@]}; do
    src=${SRC_ARRAY[$sidx]}
    dst=${DST_ARRAY[$sidx]}
    echo "Syncing files from ${src} to ${dst}"
    rsync -a $src $dst
done
