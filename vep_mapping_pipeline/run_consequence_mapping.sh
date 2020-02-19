#!/bin/bash
set -o pipefail

# A wrapper script for running the consequence mapping pipeline concurrently and based on a VCF input.
# See /README.md for detailed information.

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Explanation of GNU parallel options:
# --halt now,fail=1: if any job fails, kill the remaining ones immediately and report failure
# --pipe: input is read from STDIN and split by chunks
# -j 20: number of concurrent workers
# -N 200: number of records (lines) per worker
bcftools query "$1" -f '%CHROM:%POS:%REF:%ALT\n' \
  | sort -u \
  | parallel --halt now,fail=1 --pipe -j 20 -N 200 python3 $DIR/consequence_mapping.py \
  | sort -u > "$2"
