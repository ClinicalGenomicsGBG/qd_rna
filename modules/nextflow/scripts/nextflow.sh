#!/bin/bash

set -e -o pipefail

_clean () {
    code=$?
    ps -p $_nxf_pid &> /dev/null && {
        echo "Killing process..."
        kill -TERM $_nxf_pid
        wait $_nxf_pid
        exit $?
    } || exit $code
}

eval ${_NXF_INIT}

NXF_HOME="${TMPDIR}/.nextflow" nextflow $@ & _nxf_pid=$!
trap _clean EXIT
wait $_nxf_pid
