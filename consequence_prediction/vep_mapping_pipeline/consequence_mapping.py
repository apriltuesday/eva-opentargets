#!/usr/bin/env python3
"""Pipeline for mapping variants to the genes they affect and their functional consequences, using Ensembl VEP API. For
documentation, refer to /README.md"""

import argparse
import itertools
import json
import logging
from collections import defaultdict
from functools import lru_cache

import requests
import sys

from retry import retry

parser = argparse.ArgumentParser(description=__doc__)

logging.basicConfig()
logger = logging.getLogger('consequence_mapping')
logger.setLevel(logging.INFO)

# The "distance to the nearest gene" parameters, used to query VEP.
VEP_SHORT_QUERY_DISTANCE = 5000


def deduplicate_list(lst):
    """Removes duplicates from a list containing arbitrary (possibly unhashable) values."""
    return [element for element, _ in itertools.groupby(sorted(lst))]


def colon_based_id_to_vep_id(colon_id):
    """Converts a colon-based identifier to VEP compatible one. Example: '15:7237571:C:T' → '15 7237571 . C T'"""
    id_fields = colon_id.split(':')
    assert len(id_fields) == 4, 'Invalid colon-based identifier supplied (should contain exactly 4 fields)'
    return '{} {} . {} {}'.format(*id_fields)


def vep_id_to_colon_id(vep_id):
    """Converts a specific type of VEP compatible identifier to colon-based one. VEP supports several types of variant
    identifiers. This function only converts a single type, of the form 'CHROM POS . REF ALT', delimited by spaces and
    compatible with the first five columns of VCF."""
    vep_id_fields = vep_id.split(' ')
    return ':'.join([vep_id_fields[0], vep_id_fields[1], vep_id_fields[3], vep_id_fields[4]])


@retry(tries=10, delay=5, backoff=1.2, jitter=(1, 3), logger=logger)
def query_vep(variants, search_distance=VEP_SHORT_QUERY_DISTANCE):
    """Query VEP and return results in JSON format. Upstream/downstream genes are searched up to a given distance."""
    ensembl_request_url = 'https://rest.ensembl.org/vep/human/region'
    headers = {'Content-Type': 'application/json', 'Accept': 'application/json'}
    result = requests.post(ensembl_request_url, headers=headers, data=json.dumps({
        'variants': variants, 'distance': search_distance, 'shift_3prime': 0,
    }))

    if result.status_code == 400:
        logger.error('Bad request for the following variants:')
        logger.error(variants)

    # If there was an HTTP error, raise an exception. This will be caught by @retry
    result.raise_for_status()
    return result.json()


@lru_cache
@retry(tries=10, delay=5, backoff=1.2, jitter=(1, 3), logger=logger)
def query_consequence_types():
    url = 'https://rest.ensembl.org/info/variation/consequence_types?content-type=application/json&rank=1'
    result = requests.get(url)
    result.raise_for_status()
    return result.json()


def get_severity_ranking():
    consequence_type_results = query_consequence_types()
    # Some terms have the same rank, for these we sort lexicographically within a rank to get a stable ordering.
    ranking_dict = defaultdict(list)
    for conseq in consequence_type_results:
        ranking_dict[int(conseq['consequence_ranking'])].append(conseq['SO_term'])
    severity_ranking = []
    for rank in sorted(ranking_dict.keys()):
        severity_ranking.extend(sorted(ranking_dict[rank]))
    return severity_ranking


def load_consequence_severity_rank():
    """Loads severity rankings for consequence terms."""
    return {term: index for index, term in enumerate(get_severity_ranking())}


def extract_consequences(vep_results, acceptable_biotypes):
    """Given VEP results, return a list of consequences matching certain criteria.

    Args:
        vep_results: results obtained from VEP for a list of variants, in JSON format.
        acceptable_biotypes: a list of transcript biotypes to consider (as defined in Ensembl documentation, see
            https://www.ensembl.org/info/genome/genebuild/biotypes.html). Consequences for other transcript biotypes
            will be ignored.
    """
    consequence_term_severity_rank = load_consequence_severity_rank()
    results_by_variant = defaultdict(list)
    for result in vep_results:
        variant_identifier = result['input']
        consequences = result.get('transcript_consequences', [])

        # Keep only consequences with the allowed biotypes; if there are none, skip this variant
        consequences = [c for c in consequences if c['biotype'] in acceptable_biotypes]
        if not consequences:
            continue

        # Split by overlapping/distant and by gene
        # TODO
        overlapping_consequences = [c for c in consequences if 'distance' not in c]
        distant_consequences = [c for c in consequences if 'distance' in c]

        # TODO Apply this part of the strategy to only distant consequences
        # Flatten the list of consequence terms and find the most severe one
        all_consequence_terms = [term for c in consequences for term in c['consequence_terms']]
        all_consequence_terms.sort(key=lambda term: consequence_term_severity_rank[term])
        most_severe_consequence_term = all_consequence_terms[0]

        # Keep only consequences which include the most severe consequence term; sort by increasing order of distance.
        # If there is no 'distance' attribute in VEP results, it means that it is not applicable as the variant resides
        # *inside* the gene; hence, in this case the distance is set to 0.
        consequences = [c for c in consequences if most_severe_consequence_term in c['consequence_terms']]
        consequences.sort(key=lambda consequence: abs(consequence.get('distance', 0)))

        # Return a subset of fields (required for output) of filtered consequences
        results_by_variant[variant_identifier].extend([
            (variant_identifier, c['gene_id'], c.get('gene_symbol', ''), most_severe_consequence_term)
            for c in consequences
        ])

    return results_by_variant


def get_variants_without_consequences(results_by_variant):
    """Returns a list of variants for which no consequences were found."""
    return sorted({
        variant_id
        for variant_id, list_of_consequences in results_by_variant.items()
        if len(list_of_consequences) == 0
    })


def process_variants(variants):
    """Given a list of variant IDs, return a list of consequence types (each including Ensembl gene name & ID and a
    functional consequence code) for a given variant.
    """

    # Query VEP with default parameters, looking for variants affecting protein coding and miRNA transcripts
    # up to a standard distance (5000 nucleotides either way, which is default for VEP) from the variant.
    vep_results = query_vep(variants=variants)
    results_by_variant = extract_consequences(vep_results=vep_results, acceptable_biotypes={'protein_coding', 'miRNA'})

    # See if there are variants with no consequences up to the default distance
    variants_without_consequences = get_variants_without_consequences(results_by_variant)
    if variants_without_consequences:
        logger.info('Found {} variant(s) without standard consequences: {}'.format(
            len(variants_without_consequences), '|'.join(variants_without_consequences)))

    # Yield all consequences for all variants. Note they are not grouped by variant, all consequences are yielded in a
    # common sequence.
    for variant_id, variant_consequences in results_by_variant.items():
        for consequence_to_yield in deduplicate_list(variant_consequences):
            yield consequence_to_yield


def main():
    # Parse command line arguments
    args = parser.parse_args()

    # Load variants to query from STDIN
    variants_to_query = [colon_based_id_to_vep_id(v) for v in sys.stdin.read().splitlines()]

    # Query VEP with all variants at once (for the purpose of efficiency), print out the consequences to STDOUT.
    consequences = process_variants(variants_to_query)
    for variant_id, gene_id, gene_symbol, consequence_term in consequences:
        print('\t'.join([vep_id_to_colon_id(variant_id), gene_id, gene_symbol, consequence_term]))

    logger.info('Successfully processed {} variants'.format(len(variants_to_query)))


if __name__ == '__main__':
    main()
