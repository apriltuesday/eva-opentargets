#!/bin/bash
set -euxo pipefail

# A wrapper script for running the consequence mapping pipeline concurrently and based on a VCF input.
# See /README.md for detailed information.

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

bcftools query "$1" -f '%CHROM:%POS:%REF:%ALT\n' \
  | sort -u \
  | parallel \
    --halt now,fail=1    `# If any job fails, kill the remaining ones immediately and report failure`    \
    --pipe               `# Input is read from STDIN and split by chunks`                                \
    -j 20                `# Number of concurrent workers`                                                \
    -N 200               `# Number of records (lines) per worker`                                        \
    --tmpdir .           `# Store temporary files in the current directory to avoid /tmp overflow`       \
    python3 $DIR/vep_mapping_pipeline/consequence_mapping.py \
  | sort -u > "$2"
