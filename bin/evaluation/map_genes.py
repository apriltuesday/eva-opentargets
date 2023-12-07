#!/usr/bin/env python3
import argparse

import pandas as pd

from cmat.consequence_prediction.repeat_expansion_variants.pipeline import annotate_ensembl_gene_info
from cmat import clinvar_xml_io


def load_clinvar_data(clinvar_xml):
    """Load ClinVar data, process for gene symbols and HGNC IDs, and return it as a Pandas dataframe.
     Modified from similar functionality in the repeat expansion pipeline."""
    variant_data = []  # To populate the return dataframe (see columns below)
    for clinvar_record in clinvar_xml_io.ClinVarDataset(clinvar_xml):
        # Skip a record if it does not contain variant information
        if not clinvar_record.measure:
            continue
        measure = clinvar_record.measure

        # Extract gene symbol(s). Here and below, dashes are sometimes assigned to be compatible with the variant
        # summary format which was used previously.
        gene_symbols = measure.preferred_gene_symbols
        if not gene_symbols:
            gene_symbols = ['-']

        # Extract HGNC ID
        hgnc_ids = measure.hgnc_ids
        hgnc_id = hgnc_ids[0] if len(hgnc_ids) == 1 and len(gene_symbols) == 1 else '-'

        # Append data strings (including unused columns for compatibility)
        for gene_symbol in gene_symbols:
            variant_data.append([
                measure.preferred_or_other_name,
                clinvar_record.accession,
                gene_symbol,
                hgnc_id,
                None,
                None
            ])

    variants = pd.DataFrame(variant_data, columns=('Name',
                                                   'RCVaccession',
                                                   'GeneSymbol',
                                                   'HGNC_ID',
                                                   'TranscriptID',
                                                   'RepeatType'))

    # Since the same record can have coordinates in multiple builds, it can be repeated. Remove duplicates
    variants = variants.drop_duplicates()
    return variants


def main(clinvar_xml, output_file):
    """Load ClinVar XML, map to Ensembl gene IDs, and dump results to TSV."""
    variants = load_clinvar_data(clinvar_xml)
    # Don't include transcripts for evaluation
    annotated_variants = annotate_ensembl_gene_info(variants, include_transcripts=False)
    annotated_variants[['RCVaccession', 'EnsemblGeneID']].to_csv(output_file, sep='\t', index=False, header=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script to convert HGNC IDs and gene symbols in ClinVar to Ensembl')
    parser.add_argument('--clinvar-xml', required=True, help='ClinVar XML dump file')
    parser.add_argument('--output-file', required=True, help='File to output dataframe')
    args = parser.parse_args()
    main(args.clinvar_xml, args.output_file)
