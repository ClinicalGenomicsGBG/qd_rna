#! /bin/bash

set -veo pipefail

eval ${_SUBSAMPLE_INIT}

input_fq1_name=$(basename "${_SUBSAMPLE_INPUT_FQ1}")
input_fq2_name=$(basename "${_SUBSAMPLE_INPUT_FQ2}")
input_fq1_temp="${TMPDIR}/${input_fq1_name}"
input_fq2_temp="${TMPDIR}/${input_fq2_name}"
output_fq1_temp="${TMPDIR}/r1.fq.gz"
output_fq2_temp="${TMPDIR}/r2.fq.gz"

cp "${_SUBSAMPLE_INPUT_FQ1}" "${input_fq1_temp}"
cp "${_SUBSAMPLE_INPUT_FQ2}" "${input_fq2_temp}"

fq subsample -s 1337 -p ${_SUBSAMPLE_FRAC} --r1-dst "${output_fq1_temp}" --r2-dst "${output_fq2_temp}" "${input_fq1_temp}" "${input_fq2_temp}"

mv "${output_fq1_temp}" "${_SUBSAMPLE_OUTPUT_FQ1}"
mv "${output_fq2_temp}" "${_SUBSAMPLE_OUTPUT_FQ2}"
