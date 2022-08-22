#!/usr/bin/env python3
"""A wrapper script for running the repeat expansion pipeline."""

import argparse
import structural_variants.pipeline

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    '--clinvar-xml', required=True,
    help='ClinVar XML dump file (ClinVarFullRelease_00-latest.xml.gz)'
)
parser.add_argument(
    '--output-consequences', required=True,
    help='File to output functional consequences to. Format is compatible with the main VEP mapping pipeline.'
)

args = parser.parse_args()
structural_variants.pipeline.main(args.clinvar_xml, args.output_consequences)
