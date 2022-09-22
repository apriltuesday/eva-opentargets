#!/bin/bash

set -Eeuo pipefail
export LC_COLLATE=C

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export CODE_ROOT="$(dirname $(dirname $(dirname "${SCRIPT_DIR}")))"

# TODO put this in one place only
OT_SCHEMA_VERSION=2.2.8

# TODO should come from GH actions
export PYTHON_BIN=${CODE_ROOT}/env/bin/python
export BATCH_ROOT_BASE=${SCRIPT_DIR}/resources/end2end

CWD=${PWD}
BATCH_ROOT=${BATCH_ROOT_BASE}/test_batch
mkdir -p ${BATCH_ROOT}
cd ${BATCH_ROOT}

nextflow run ${CODE_ROOT}/eva_cttv_pipeline/evidence_string_generation/pipeline.nf \
  --batch_root ${BATCH_ROOT} \
  --schema ${OT_SCHEMA_VERSION} \
  --clinvar ${BATCH_ROOT_BASE}/input.xml.gz \
  -resume

diff ${BATCH_ROOT}/gene_mapping/consequences_snp.tsv ${BATCH_ROOT_BASE}/expected/consequences_snp.tsv
diff ${BATCH_ROOT}/gene_mapping/consequences_repeat.tsv ${BATCH_ROOT_BASE}/expected/consequences_repeat.tsv
diff ${BATCH_ROOT}/gene_mapping/consequences_structural.tsv ${BATCH_ROOT_BASE}/expected/consequences_structural.tsv
diff ${BATCH_ROOT}/evidence_strings/evidence_strings.json ${BATCH_ROOT_BASE}/expected/evidence_strings.json

cd ${CWD}
rm -r ${BATCH_ROOT}
