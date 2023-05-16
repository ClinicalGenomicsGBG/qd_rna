#!/bin/bash

set -e -o pipefail

module load $_JAVA_MODULE
module load $_NXF_MODULE
NXF_HOME="${TMPDIR}/.nextflow" nextflow $@
