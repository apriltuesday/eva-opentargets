#!/usr/bin/env python3
"""A wrapper script for running the repeat expansion pipeline."""

import argparse
import repeat_expansion_variants.pipeline

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    '--clinvar-summary-tsv', required=True,
    help='ClinVar summary TSV dump file (variant_summary.txt.gz)'
)
parser.add_argument(
    '--output-consequences', required=True,
    help='File to output functional consequences to. Format is compatible with the main VEP mapping pipeline.'
)
parser.add_argument(
    '--output-dataframe', required=True,
    help='File to output full dataframe for subsequent analysis and debugging.'
)
args = parser.parse_args()
repeat_expansion_variants.pipeline.main(args.clinvar_summary_tsv, args.output_consequences, args.output_dataframe)
#
#
# DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
# # Awk expression prints the first line, and also all lines where column 2 (variant type) is "NT expansion"
# python3 "${DIR}/repeat_expansion_variants/pipeline.py" \
#   --clinvar-summary-tsv <(zcat "${VARIANT_SUMMARY}" | awk -F$'\t' '(NR==1) || ($2 == "NT expansion")') \
#   --output-consequences "${OUTPUT_CONSEQUENCES}" \
#   --output-dataframe "${OUTPUT_DATAFRAME}"
