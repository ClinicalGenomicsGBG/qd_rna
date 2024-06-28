#!/bin/bash
set -e -o pipefail

while IFS= read -r line; do
  printf 'rsync %s\n' "$line"
  rsync -a --no-compress $line
done < ${MANIFEST}
