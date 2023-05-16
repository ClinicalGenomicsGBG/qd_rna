#!/bin/bash
set -e -o pipefail

SRC_ARRAY=( $SRC )
DST_ARRAY=( $DST )

for didx in ${!DST_ARRAY[@]}; do
    src=$(tr ",", " " <<< ${SRC_ARRAY[$didx]})
    dst=${DST_ARRAY[$didx]}
    echo "Syncing ${src} to ${dst}"
    rsync -a $src $dst
done
