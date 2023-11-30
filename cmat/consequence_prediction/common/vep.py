import itertools
import json
import logging
from collections import defaultdict
from functools import lru_cache

import requests

from retry import retry

logging.basicConfig()
logger = logging.getLogger('consequence_mapping')
logger.setLevel(logging.INFO)

# The "distance to the nearest gene" parameters, used to query VEP.
VEP_SHORT_QUERY_DISTANCE = 5000


def deduplicate_list(lst):
    """Removes duplicates from a list containing arbitrary (possibly unhashable) values."""
    return [element for element, _ in itertools.groupby(sorted(lst))]


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


def most_severe_consequence(consequence_terms_with_transcripts, consequence_term_severity_rank):
    return min(consequence_terms_with_transcripts,
               key=lambda term_and_transcript: consequence_term_severity_rank[term_and_transcript[0]])[0]


def most_severe_consequence_per_gene(variant_identifier, consequences, include_transcripts):
    results = []
    consequence_term_severity_rank = load_consequence_severity_rank()
    consequences_per_gene = defaultdict(list)
    for c in consequences:
        key = (c['gene_id'], c.get('gene_symbol', ''))
        # Keep track of consequence term alongside transcript ID
        consequences_per_gene[key].extend((term, c['transcript_id']) for term in c['consequence_terms'])
    for (gene_id, gene_symbol), terms_with_transcripts in consequences_per_gene.items():
        most_severe_consequence_term = most_severe_consequence(terms_with_transcripts, consequence_term_severity_rank)
        # If we're including transcripts, need to include every transcript associated with the most severe consequence
        if include_transcripts:
            for term, transcript_id in terms_with_transcripts:
                if term == most_severe_consequence_term:
                    results.append((variant_identifier, gene_id, gene_symbol, most_severe_consequence_term, transcript_id))
        else:
            results.append((variant_identifier, gene_id, gene_symbol, most_severe_consequence_term))
    return results


def overall_most_severe_consequence(variant_identifier, consequences, include_transcripts):
    results = []
    consequence_term_severity_rank = load_consequence_severity_rank()
    # Flatten the list of consequence terms and find the most severe one
    all_consequence_terms = [(term, c['transcript_id']) for c in consequences for term in c['consequence_terms']]
    most_severe_consequence_term = most_severe_consequence(all_consequence_terms, consequence_term_severity_rank)

    # Keep only consequences which include the most severe consequence term.
    for c in consequences:
        if most_severe_consequence_term in c['consequence_terms']:
            if include_transcripts:
                results.append((variant_identifier, c['gene_id'], c.get('gene_symbol', ''), most_severe_consequence_term, c['transcript_id']))
            else:
                results.append((variant_identifier, c['gene_id'], c.get('gene_symbol', ''), most_severe_consequence_term))
    return results


def extract_consequences(vep_results, acceptable_biotypes, include_transcripts):
    """Given VEP results, return a list of consequences matching certain criteria.

    Args:
        vep_results: results obtained from VEP for a list of variants, in JSON format.
        acceptable_biotypes: a list of transcript biotypes to consider (as defined in Ensembl documentation, see
            https://www.ensembl.org/info/genome/genebuild/biotypes.html). Consequences for other transcript biotypes
            will be ignored.
        include_transcripts: whether to include transcript IDs alongside consequence terms
    """
    results_by_variant = defaultdict(list)
    for result in vep_results:
        variant_identifier = result['input']
        consequences = result.get('transcript_consequences', [])

        # Keep only consequences with the allowed biotypes; if there are none, skip this variant
        consequences = [c for c in consequences if c['biotype'] in acceptable_biotypes]
        if not consequences:
            continue

        # If there is no 'distance' attribute in VEP results, it means the variant overlaps the gene.
        overlapping_consequences = [c for c in consequences if 'distance' not in c]

        # For genes overlapping the variant, we report the most severe consequence per gene.
        if overlapping_consequences:
            consequences_for_variant = most_severe_consequence_per_gene(variant_identifier, overlapping_consequences, include_transcripts)
        # If there are no consequences on overlapping genes, we take the overall most severe consequence and all genes
        # associated with that consequence
        else:
            consequences_for_variant = overall_most_severe_consequence(variant_identifier, consequences, include_transcripts)
        results_by_variant[variant_identifier].extend(consequences_for_variant)

    return results_by_variant
