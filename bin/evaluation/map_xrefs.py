#!/usr/bin/env python3
import argparse
import csv
import multiprocessing
from functools import lru_cache

from eva_cttv_pipeline.clinvar_xml_io import clinvar_xml_io
from eva_cttv_pipeline.clinvar_xml_io.clinvar_xml_io.ontology_uri import OntologyUri
from eva_cttv_pipeline.trait_mapping.ols import build_ols_query
from eva_cttv_pipeline.trait_mapping.utils import json_request


@lru_cache
def get_synonyms(db, iden):
    """Find synonyms (replacement terms or exact matches) for this ontology identifier using OLS."""
    synonyms = set()
    ontology_uri = OntologyUri(iden, db)
    url = build_ols_query(str(ontology_uri))
    json_response = json_request(url)
    if json_response and '_embedded' in json_response:
        for term in json_response['_embedded']['terms']:
            # Get only EFO terms (even if imported)
            if term['ontology_name'] == 'efo':
                synonyms.add(term['iri'])
                # Check whether current term is obsolete
                if term['is_obsolete'] and term['term_replaced_by']:
                    synonyms.add(term['term_replaced_by'])
                # Also add exact matches
                if 'exactMatch' in term['annotation']:
                    synonyms.update(term['annotation']['exactMatch'])

        # Synonyms contains current EFO-included URIs, convert to DB:ID style
        synonyms = {OntologyUri.uri_to_curie(s) for s in synonyms}
        # Filter out Nones and sort lexicographically
        synonyms = {s for s in synonyms if s is not None}

    if synonyms:
        return ontology_uri.curie, synonyms
    # If no synonyms, just return the original identifier as is
    return ontology_uri.curie, {ontology_uri.curie}


def main(clinvar_xml, output_file):
    """Load ClinVar XML, map trait xrefs identifiers to synonyms in OLS, and dump results to TSV."""
    traits = set()
    for record in clinvar_xml_io.ClinVarDataset(clinvar_xml):
        for trait in record.traits_with_valid_names:
            traits.update([(db, iden) for db, iden, _ in trait.current_efo_aligned_xrefs])

    traits = list(traits)
    process_pool = multiprocessing.Pool(processes=24)
    annotated_traits = [
        process_pool.apply(get_synonyms, args=(db, iden))
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
