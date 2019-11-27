#!/usr/bin/env python3

import argparse
import re
import requests

# Name of ontology in OLS url, e. g. https://www.ebi.ac.uk/ols/ontologies/ordo/terms?iri=...
ontology_to_ols = {
    'HP': 'hp',
    'MONDO': 'mondo',
    'Orphanet': 'ordo',
}

# OLS url to query for a term details
ols_url_template = 'https://www.ebi.ac.uk/ols/api/ontologies/{ontology}/terms?iri={term}'

# OxO url to query ontology cross-references
oxo_url_template = 'https://www.ebi.ac.uk/spot/oxo/api/search?ids={curie}&distance=1&size=500'

# List of fields in the current version of Webulous submission template
webulous_fields = [
    'disease', 'child_of', 'definition', 'synonyms', 'located_in_organ', 'located_in_cell',
    'biological_process', 'msh_def_cite', 'ncit_def_cite', 'snomedct_def_cite', 'icd9_def_cite',
    'icd10_def_cite', 'omim_def_cite', 'doid_def_cite', 'meddra_def_cite', 'umls_def_cite',
    'wikipedia_def_cite', 'comments', 'ordo_def_cite', 'definition_editor', 'definition_citation',
    'mondo_def_cite'
]
webulous_format_string = '\t'.join('{' + f + '}' for f in webulous_fields) + '\n'

# String to join multiple values in a Webulous template
webulous_joiner = ' || '


def get_parent_terms(url):
    return [term['label'] for term in requests.get(url).json()['_embedded']['terms']]


def uri_to_curie(uri):
    """Converts URI to curie (short identifier).

    Args:
        uri: a full URI of an ontology term, e. g. http://purl.obolibrary.org/obo/MONDO_0009796. URIs are globally
            unique among all ontologies (and even other internet resources).

    Returns:
        curie: a short identifier (Compact URI) which contains an ontology prefix and identifier in that ontology.
            Example: MONDO:0009796. See also: http://www.obofoundry.org/docs/Citation.html
    """
    return uri.split('/')[-1].replace('#', '').replace('_', ':')


def get_cross_references(curie):
    """Queries OxO to return the list of cross-references for a given term curie."""
    url = oxo_url_template.format(curie=curie)
    mappings = requests.get(url).json()['_embedded']['searchResults'][0]['mappingResponseList']
    return [m['curie'] for m in mappings]


def get_ols_details(ontology, term):
    """Queries OLS and returns the details necessary for the EFO import table construction."""
    url = ols_url_template.format(ontology=ontology, term=term)
    data = requests.get(url).json()['_embedded']['terms'][0]
    label = data['label']
    parents = get_parent_terms(data['_links']['parents']['href'])

    # Definition: first, try to get from annotation field
    definition = data['annotation'].get('definition', [''])[0]
    # If no luck, simply use a description
    if not definition and data['description']:
        definition = data['description'][0]

    synonyms = data['synonyms'] or []

    # Cross-references
    term_curie = uri_to_curie(term)
    xrefs = {}
    for x in get_cross_references(term_curie):
        xref_ontology, xref_id = re.split('[_:]', x)
        xrefs.setdefault(xref_ontology, set()).add('{}:{}'.format(xref_ontology, xref_id))

    # If a term comes from either Orphanet or MONDO, we need to add these as xrefs as well
    # (since they won't be present in the normal list of xrefs).
    if ontology == 'mondo':
        xrefs.setdefault('MONDO', set()).add(term.split('/')[-1].replace('_', ':'))
    elif ontology == 'ordo':
        xrefs.setdefault('Orphanet', set()).add(term.split('/')[-1].replace('_', ':'))

    return label, parents, definition, synonyms, xrefs


def format_xref(xrefs, ontology):
    if ontology not in xrefs:
        return ''
    return webulous_joiner.join(sorted(xrefs[ontology]))


def format_output_string(ontology, term):
    label, parents, definition, synonyms, xrefs = get_ols_details(ontology, term)

    return webulous_format_string.format(
        disease=label,
        child_of=webulous_joiner.join(parents),
        definition=definition,
        synonyms=webulous_joiner.join(synonyms),
        located_in_organ='',
        located_in_cell='',
        biological_process='',
        msh_def_cite=format_xref(xrefs, 'MeSH'),
        ncit_def_cite=format_xref(xrefs, 'NCIT'),
        snomedct_def_cite='',
        icd9_def_cite='',
        icd10_def_cite='',
        omim_def_cite=format_xref(xrefs, 'OMIM'),
        doid_def_cite=format_xref(xrefs, 'DOID'),
        meddra_def_cite='',
        umls_def_cite=format_xref(xrefs, 'UMLS'),
        wikipedia_def_cite='',
        comments='',
        ordo_def_cite=format_xref(xrefs, 'Orphanet'),
        definition_editor='',
        definition_citation='',
        mondo_def_cite=format_xref(xrefs, 'MONDO'),
    )


def create_efo_table(input_file_path, output_file_path):
    with open(input_file_path) as infile, open(output_file_path, 'w') as outfile:
        for line in infile:
            term = line.rstrip()
            print('Processing ' + term)
            ontology = ontology_to_ols[re.split('[:_]', term.split('/')[-1])[0]]
            result = format_output_string(ontology, term)
            outfile.write(result)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-mappings', required=True,
                        help='Input file with ontology mappings')
    parser.add_argument('-o', '--output', required=True,
                        help='Output table for EFO import')
    args = parser.parse_args()
    create_efo_table(args.input_mappings, args.output)
