#!/usr/bin/env python3

import argparse
from cmat.output_generation import clinvar_to_evidence_strings

parser = argparse.ArgumentParser('Generates Open Targets evidence strings from ClinVar data and trait mappings')
parser.add_argument('--clinvar-xml',  help='ClinVar XML release',                    required=True)
parser.add_argument('--efo-mapping',  help='Disease string to ontology mappings',    required=True)
parser.add_argument('--gene-mapping', help='Variant to gene & consequence mappings', required=True)
parser.add_argument('--ot-schema',    help='OpenTargets schema JSON',                required=True)
parser.add_argument('--out',          help='Output directory',                       required=True)
parser.add_argument('--start',        help='Start index (inclusive)',                required=False, type=int)
parser.add_argument('--end',          help='End index (exclusive)',                  required=False, type=int)


if __name__ == '__main__':
    args = parser.parse_args()
    clinvar_to_evidence_strings.launch_pipeline(
        clinvar_xml_file=args.clinvar_xml, efo_mapping_file=args.efo_mapping, gene_mapping_file=args.gene_mapping,
        ot_schema_file=args.ot_schema, dir_out=args.out, start=args.start, end=args.end)
