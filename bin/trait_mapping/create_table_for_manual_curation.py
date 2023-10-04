#!/usr/bin/env python3

import argparse

import pandas as pd

from cmat.output_generation.clinvar_to_evidence_strings import load_efo_mapping
from cmat.trait_mapping.ols import (
    get_ontology_label_from_ols, is_current_and_in_efo, is_in_efo, get_replacement_term,
)


def previous_and_replacement_mappings(trait_name, previous_mappings):
    if trait_name not in previous_mappings:
        yield '', ''
        return
    for uri, label in previous_mappings[trait_name]:
        trait_status = get_trait_status(uri)
        trait_string = '|'.join([uri, label, 'NOT_SPECIFIED', 'previously-used', trait_status])
        replacement_string = find_replacement_mapping(uri)
        yield trait_string, replacement_string


def find_replacement_mapping(previous_uri):
    replacement_uri = get_replacement_term(previous_uri)
    if not replacement_uri:
        return ''
    label = get_ontology_label(replacement_uri)
    trait_status = get_trait_status(replacement_uri)
    trait_string = '|'.join([replacement_uri, label, 'NOT_SPECIFIED', 'replacement', trait_status])
    return trait_string


def find_exact_mapping(trait_name, mappings):
    for mapping in mappings:
        if mapping.lower().split('|')[1] == trait_name:
            return mapping
    return ''


def get_ontology_label(uri):
    label = get_ontology_label_from_ols(uri)
    return label if label is not None else ''


def get_trait_status(uri):
    uri_is_current_and_in_efo = is_current_and_in_efo(uri)
    uri_in_efo = is_in_efo(uri)
    if uri_in_efo:
        trait_status = 'EFO_CURRENT' if uri_is_current_and_in_efo else 'EFO_OBSOLETE'
    else:
        trait_status = 'NOT_CONTAINED'
    return trait_status


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-t', '--traits-for-curation',
        help='Table with traits for which the pipeline failed to make a confident prediction')
    parser.add_argument(
        '-m', '--previous-mappings',
        help='Table with all mappings previously issued by EVA. TSV with columns: ClinVar trait name; ontology URI; '
             'ontology label (not used)')
    parser.add_argument(
        '-c', '--previous-comments',
        help='Table with last round of curator comments. TSV with columns: ClinVar trait name; comments')
    parser.add_argument(
        '-o', '--output',
        help='Output TSV to be loaded in Google Sheets for manual curation')
    args = parser.parse_args()

    # Load all previous mappings: ClinVar trait name to ontology URI
    previous_mappings = load_efo_mapping(args.previous_mappings)

    # Load previous curator comments: ClinVar trait name to comment string
    try:
        previous_comments = pd.read_csv(args.previous_comments, sep='\t', header=None)
        previous_comments = dict(zip(previous_comments[0], previous_comments[1]))
    except pd.errors.EmptyDataError:
        previous_comments = {}

    # Process all mappings which require manual curation
    rows = []
    for line in open(args.traits_for_curation):
        fields = line.split('\t')
        fields[-1] = fields[-1].rstrip()  # To avoid stripping the entire field if it's empty
        trait_name, trait_freq, notes = fields[:3]
        # Add previous curator comment if present
        if trait_name in previous_comments:
            notes = f'"{notes}\n{previous_comments[trait_name]}"'
        # Use maximum of 50 mappings to improve Google Sheets performance
        mappings = fields[3:53]
        exact_mapping = find_exact_mapping(trait_name, mappings)
        for previous_mapping, replacement_mapping in previous_and_replacement_mappings(trait_name, previous_mappings):
            rows.append([trait_name, trait_freq, notes, previous_mapping, exact_mapping, replacement_mapping]
                        + mappings)

    rows.sort(key=lambda x: (x[2], int(x[1])), reverse=True)
    with open(args.output, 'w') as outfile:
        for row in rows:
            out_line = '\t'.join(row) + '\n'
            outfile.write(out_line)
