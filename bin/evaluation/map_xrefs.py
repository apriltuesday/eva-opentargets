#!/usr/bin/env python3
import argparse
import csv
import multiprocessing

from cmat import clinvar_xml_io
from cmat.output_generation.evaluation.ols_utils import fetch_eval_data


def main(clinvar_xml, output_file):
    """Load ClinVar XML, map trait xrefs identifiers to synonyms in OLS, and dump results to TSV."""
    traits = set()
    for record in clinvar_xml_io.ClinVarDataset(clinvar_xml):
        for trait in record.traits_with_valid_names:
            traits.update([(db, iden) for db, iden, _ in trait.current_efo_aligned_xrefs])

    traits = list(traits)
    process_pool = multiprocessing.Pool(processes=24)
    annotated_traits = [
        process_pool.apply(fetch_eval_data, kwds={'db_iden': (db, iden), 'include_neighbors': True})
        for db, iden in traits
    ]
    with open(output_file, 'w+') as f:
        csv.writer(f, delimiter="\t").writerows(annotated_traits)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script to map trait xrefs in ClinVar to synonyms in OLS')
    parser.add_argument('--clinvar-xml', required=True, help='ClinVar XML dump file')
    parser.add_argument('--output-file', required=True, help='File to output dataframe')
    args = parser.parse_args()
    main(args.clinvar_xml, args.output_file)
