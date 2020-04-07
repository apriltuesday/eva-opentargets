#!/bin/bash
set -euxo pipefail

# A wrapper script to run the repeat_expansion_variants script.
# Filters out only "NT expansion" variants from the ClinVar data file to avoid loading it all into memory

# Usage: bash run_repeat_expansion_variants.sh variant_summary.txt.gz output_consequences.tsv
export VARIANT_SUMMARY="$1"
export OUTPUT_CONSEQUENCES="$2"

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
# Awk expression prints the first line, and also all lines where column 2 (variant type) is "NT expansion"
python3 "${DIR}/repeat_expansion_variants/pipeline.py" \
  --clinvar-summary-tsv <(zcat "${VARIANT_SUMMARY}" | awk -F$'\t' '(NR==1) || ($2 == "NT expansion")') \
  --output-consequences "${OUTPUT_CONSEQUENCES}"
