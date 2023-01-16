#!/bin/bash
set -e -o pipefail

source ${_MODULES_INIT}
module load nextflow
HOME="$(pwd)" nextflow $@ & _nf_pid=$!
trap "echo 'Killing nextflow...'; kill -TERM $_nf_pid; wait $_nf_pid" 2 3 6 10 12 14 15
wait $_nf_pid
ß