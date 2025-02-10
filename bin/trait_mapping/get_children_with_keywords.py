#!/usr/bin/env python3

import argparse
from collections import defaultdict

from cmat.clinvar_xml_io.ontology_uri import OntologyUri
from cmat.trait_mapping.ols import build_ols_query
from cmat.trait_mapping.utils import json_request


def append_embedded(results, json_response):
    if json_response and '_embedded' in json_response:
        for key in json_response['_embedded']:
            results[key].extend(json_response['_embedded'][key])


def query_and_depaginate(url):
    json_response = json_request(url)
    results = defaultdict(list)
    append_embedded(results, json_response)
    while 'next' in json_response['_links']:
        json_response = json_request(json_response['_links']['next']['href'])
        append_embedded(results, json_response)
    return results


def search_in(keywords, text):
    return set((keyword for keyword in keywords if keyword in text))


def main():
    parser = argparse.ArgumentParser('Search OLS for children of a term that match certain keywords in their label, description or synonyms')
    parser.add_argument('--ontology', type=str, default='MONDO', help='Name of the Ontology to find the parent and children')
    parser.add_argument('--parent_curie', type=str, help='Curie of the parent term', required=True)
    parser.add_argument('--keywords', type=str, nargs='+', help="Words that must be present in the child's ontology label, description or synonyms to be reported")
    args = parser.parse_args()
    keywords = set(args.keywords)

    db = args.ontology
    parent_curie = args.parent_curie
    url = build_ols_query(OntologyUri(parent_curie, db).uri)
    results = query_and_depaginate(url)
    for term in results['terms']:
        if term['ontology_prefix'] == db:
            children_results = query_and_depaginate(term['_links']['children']['href'])
            for child_term in children_results['terms']:
                    if child_term['ontology_prefix'] == db:
                        keyword_found = set()
                        keyword_found.update(search_in(keywords, child_term['label']))
                        keyword_found.update(search_in(keywords, child_term['description']))
                        for synonym in child_term['synonyms']:
                            keyword_found.update(search_in(keywords, synonym))
                        if keyword_found == keywords:
                            print(child_term['iri'], child_term['label'])


if __name__ == '__main__':
    main()

