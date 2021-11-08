#!/usr/bin/env python3

import argparse
import re

from eva_cttv_pipeline.clinvar_xml_utils.clinvar_xml_utils import clinvar_xml_utils

parser = argparse.ArgumentParser('Processes ClinVar XML dump and extract all variants, in CHROM:POS:REF:ALT format,'
                                 'for processing by the VEP mapping pipeline.')
parser.add_argument('--clinvar-xml', required=True, help='Path to the ClinVar XML file')
args = parser.parse_args()

# A regular expression to detect alleles with IUPAC ambiguity bases.
IUPAC_AMBIGUOUS_SEQUENCE = re.compile(r'[^ACGT]')

for clinvar_record in clinvar_xml_utils.ClinVarDataset(args.clinvar_xml):
    if clinvar_record.measure is None or not clinvar_record.measure.has_complete_coordinates:
        continue
    m = clinvar_record.measure
    if IUPAC_AMBIGUOUS_SEQUENCE.search(m.vcf_ref + m.vcf_alt):
        continue
    print(f'{m.chr}:{m.vcf_pos}:{m.vcf_ref}:{m.vcf_alt}')
