import logging
from functools import lru_cache

from requests import RequestException

from cmat.clinvar_xml_io.ontology_uri import OntologyUri
from cmat.trait_mapping.utils import json_request
from cmat.trait_mapping.ols import build_ols_query

logger = logging.getLogger(__package__)


@lru_cache
def fetch_eval_data(*, db_iden=None, uri=None, include_neighbors=False, target_ontology='EFO'):
    """
    Query OLS for this ontology identifier and extract the following:
    - Whether the term is obsolete in target ontology (default EFO)
    - Synonyms (replacement terms or exact matches)
    - Parents and children in target ontology
    """
    if db_iden:
        db, iden = db_iden
        ontology_uri = OntologyUri(iden, db).uri
    elif uri:
        ontology_uri = uri
    else:
        logger.warning("Must provide either DB + identifier or full URI")
        return None
    curie = OntologyUri.uri_to_curie(ontology_uri)

    # Defaults to return if OLS query fails or no term in target ontology
    is_obsolete = False
    synonyms = {}
    parents = {}
    children = {}

    url = build_ols_query(ontology_uri)
    try:
        json_response = json_request(url)
    except RequestException:
        logger.warning(f'OLS4 error for {url}, trying OLS3...')
        json_response = json_request(url.replace('/ols4/', '/ols/'))
    if json_response and '_embedded' in json_response:
        for term in json_response['_embedded']['terms']:
            # Get only target ontology terms (even if imported)
            if term['ontology_name'] == target_ontology:
                synonyms, is_obsolete = extract_synonyms_and_obsolete(term)
                # If requested, fetch the parents and children of this term
                if include_neighbors:
                    parents, children = extract_parents_and_children(term)

    if include_neighbors:
        return curie, is_obsolete, synonyms, parents, children
    return curie, is_obsolete, synonyms


def extract_synonyms_and_obsolete(ontology_term):
    synonyms = {ontology_term['iri']}
    is_obsolete = ontology_term['is_obsolete']

    # Add replacement term if this one is obsolete
    if is_obsolete and ontology_term['term_replaced_by']:
        synonyms.add(ontology_term['term_replaced_by'])
    # Also add exact matches
    if 'exactMatch' in ontology_term['annotation']:
        synonyms.update(ontology_term['annotation']['exactMatch'])
    if 'has exact match' in ontology_term['annotation']:
        synonyms.update(ontology_term['annotation']['has exact match'])

    # Synonyms contains current included URIs, convert to DB:ID style
    synonyms = {OntologyUri.uri_to_curie(s) for s in synonyms}
    # Filter out Nones
    synonyms = {s for s in synonyms if s is not None}
    return synonyms, is_obsolete


def extract_parents_and_children(ontology_term):
    links = ontology_term['_links']
    parents = get_all_term_curies(links['parents']['href']) if 'parents' in links else {}
    children = get_all_term_curies(links['children']['href']) if 'children' in links else {}
    return parents, children


def get_all_term_curies(url):
    curies = set()
    json_response = json_request(url)
    if json_response and '_embedded' in json_response:
        for term in json_response['_embedded']['terms']:
            curies.add(term['obo_id'])
    return curies
