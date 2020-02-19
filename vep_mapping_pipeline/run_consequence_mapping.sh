#!/bin/bash

# A wrapper script for running the consequence mapping pipeline concurrently and based on a VCF input.
# See /README.md for detailed information.

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

bcftools query "$1" -f '%CHROM:%POS:%REF:%ALT\n' \
  | sort -u \
  | parallel --pipe -j 20 -N 200 python3 $DIR/consequence_mapping.py \
  | sort -u > "$2"
