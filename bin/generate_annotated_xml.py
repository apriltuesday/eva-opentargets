#!/usr/bin/env python3

import argparse
from cmat.output_generation.annotated_clinvar import generate_annotated_clinvar_xml

parser = argparse.ArgumentParser('Generates annotated ClinVar XML from ClinVar data and trait mappings')
parser.add_argument('--clinvar-xml',  help='ClinVar XML release',                    required=True)
parser.add_argument('--trait-mapping', help='Disease string to ontology mappings',   required=True)
parser.add_argument('--gene-mapping', help='Variant to gene & consequence mappings', required=True)
parser.add_argument('--output-xml',   help='Output XML file',                        required=True)
parser.add_argument('--eval-gene-file', help='Gene mappings for evaluation',         required=False)
parser.add_argument('--eval-xref-file', help='Trait xref mappings for evaluation',   required=False)
parser.add_argument('--eval-latest-file', help='Latest mapping annotations for evaluation', required=False)


if __name__ == '__main__':
    args = parser.parse_args()
    generate_annotated_clinvar_xml(
        clinvar_xml_file=args.clinvar_xml, trait_mapping_file=args.trait_mapping, gene_mapping_file=args.gene_mapping,
        output_xml_file=args.output_xml, eval_gene_file=args.eval_gene_file, eval_xref_file=args.eval_xref_file,
        eval_latest_file=args.eval_latest_file)
