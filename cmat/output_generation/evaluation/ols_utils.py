import logging
from functools import lru_cache

from cmat.clinvar_xml_io.ontology_uri import OntologyUri
from cmat.trait_mapping.utils import json_request
from cmat.trait_mapping.ols import build_ols_query

logger = logging.getLogger(__package__)


@lru_cache
def fetch_eval_data(*, db_iden=None, uri=None, include_neighbors=False):
    """
    Query OLS for this ontology identifier and extract the following:
    - Whether the term is obsolete in EFO
    - Synonyms (replacement terms or exact matches)
    - Parents and children in EFO
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

    # Defaults to return if OLS query fails
    is_obsolete = False
    synonyms = {curie}
    parents = {}
    children = {}

    url = build_ols_query(ontology_uri)
    json_response = json_request(url)
    if json_response and '_embedded' in json_response:
        for term in json_response['_embedded']['terms']:
            # Get only EFO terms (even if imported)
            if term['ontology_name'] == 'efo':
                synonyms, is_obsolete = extract_synonyms_and_obsolete(term)
                # If requested, fetch the parents and children of this term
                if include_neighbors:
                    parents, children = extract_parents_and_children(term)

    if include_neighbors:
        return curie, is_obsolete, synonyms, parents, children
    return curie, is_obsolete, synonyms


def extract_synonyms_and_obsolete(efo_term):
    synonyms = {efo_term['iri']}
    is_obsolete = efo_term['is_obsolete']

    # Add replacement term if this one is obsolete
    if is_obsolete and efo_term['term_replaced_by']:
        synonyms.add(efo_term['term_replaced_by'])
    # Also add exact matches
    if 'exactMatch' in efo_term['annotation']:
        synonyms.update(efo_term['annotation']['exactMatch'])
    if 'has exact match' in efo_term['annotation']:
        synonyms.update(efo_term['annotation']['has exact match'])

    # Synonyms contains current EFO-included URIs, convert to DB:ID style
    synonyms = {OntologyUri.uri_to_curie(s) for s in synonyms}
    # Filter out Nones
    synonyms = {s for s in synonyms if s is not None}
    return synonyms, is_obsolete


def extract_parents_and_children(efo_term):
    links = efo_term['_links']
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
