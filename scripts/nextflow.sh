#!/bin/bash
set -e -o pipefail

source ${_MODULES_INIT}
module load nextflow
HOME="$(pwd)" nextflow $@ & _nf_pid=$!

trap "kill -0 $_nf_pid && echo 'Killing nextflow...' && kill -TERM $_nf_pid; wait $_nf_pid" EXIT

wait $_nf_pid
exit $?