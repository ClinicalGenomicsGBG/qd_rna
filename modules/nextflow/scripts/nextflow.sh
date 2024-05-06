#!/bin/bash

set -meo pipefail

_clean () {
    code=$?
    kill -0 $_nxf_pid && {
        echo "Killing process ${_nxf_pid}..."
        kill -TERM $_nxf_pid
        wait $_nxf_pid
    } || exit $code
}

trap _clean HUP INT QUIT ABRT USR1 USR2 ALRM TERM

eval "${_NXF_INIT}"
NXF_HOME="${TMPDIR}/.nextflow" nextflow $@ & _nxf_pid=$!
wait $_nxf_pid