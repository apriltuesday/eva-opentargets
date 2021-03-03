#!/usr/bin/env python3

import argparse

from eva_cttv_pipeline.trait_mapping.ols import (
    get_ontology_label_from_ols, is_current_and_in_efo, is_in_efo,
)


def find_previous_mapping(trait_name, previous_mappings):
    if trait_name not in previous_mappings:
        return ''
    uri = previous_mappings[trait_name]
    label = get_ontology_label_from_ols(uri)
    uri_is_current_and_in_efo = is_current_and_in_efo(uri)
    uri_in_efo = is_in_efo(uri)
    if uri_in_efo:
        trait_status = 'EFO_CURRENT' if uri_is_current_and_in_efo else 'EFO_OBSOLETE'
    else:
        trait_status = 'NOT_CONTAINED'
    trait_string = '|'.join([uri, label, 'NOT_SPECIFIED', 'previously-used', trait_status])
    return trait_string


def find_exact_mapping(trait_name, mappings):
    for mapping in mappings:
        if mapping.lower().split('|')[1] == trait_name:
            return mapping
    return ''


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
        '-o', '--output',
        help='Output TSV to be loaded in Google Sheets for manual curation')
    args = parser.parse_args()
    outfile = open(args.output, 'w')

    # Load all previous mappings: ClinVar trait name to ontology URI
    previous_mappings = dict(line.rstrip().split('\t')[:2] for line in open(args.previous_mappings))

    # Process all mappings which require manual curation
    for line in open(args.traits_for_curation):
        fields = line.split('\t')
        fields[-1] = fields[-1].rstrip()  # To avoid stripping the entire field if it's empty
        trait_name, trait_freq, notes = fields[:3]
        mappings = fields[3:]
        previous_mapping = find_previous_mapping(trait_name, previous_mappings)
        exact_mapping = find_exact_mapping(trait_name, mappings)
        out_line = '\t'.join(
            [trait_name, trait_freq, notes, previous_mapping, exact_mapping] + mappings
        ) + '\n'
        outfile.write(out_line)

    outfile.close()
