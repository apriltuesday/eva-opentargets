from functools import lru_cache

from eva_cttv_pipeline.clinvar_xml_io.clinvar_xml_io.ontology_uri import OntologyUri
from eva_cttv_pipeline.trait_mapping.ols import build_ols_query
from eva_cttv_pipeline.trait_mapping.utils import json_request


@lru_cache
def get_canonical_id(db, iden):
    """
    Choose a canonical representative of the set of synonyms (replacement terms or exact matches)
    for this ontology identifier, using the following ranking:
    1) Prefer current terms over obsolete terms
    2) Among current terms, choose the one with DB:ID lexicographically first
    """
    synonyms = set()
    ontology_uri = OntologyUri(iden, db)
    url = build_ols_query(str(ontology_uri))
    json_response = json_request(url)
    if not json_response or '_embedded' not in json_response:
        return None
    for term in json_response['_embedded']['terms']:
        # Get only EFO terms (even if imported)
        if term['ontology_name'] == 'efo':
            # Check whether current term is obsolete
            if term['is_obsolete'] and term['term_replaced_by']:
                synonyms.add(term['term_replaced_by'])
            else:
                synonyms.add(term['iri'])
            # Also add exact matches
            # TODO does not check if exact matches are themselves obsolete!
            if 'exactMatch' in term['annotation']:
                synonyms.update(term['annotation']['exactMatch'])

    # Synonyms contains current EFO-included URIs, convert to DB:ID style
    synonyms = {uri_to_curie(s) for s in synonyms}
    # Filter out Nones and sort lexicographically
    synonyms = sorted([s for s in synonyms if s is not None])
    # If there's nothing left just return the original identifier as is, otherwise return the first
    if synonyms:
        return synonyms[0]
    elif db.lower() == 'orphanet':
        return f'Orphanet:{iden}'
    return iden


def uri_to_curie(uri):
    """Convert an ontology uri to a DB:ID format."""
    URI_DB_TO_DB_DICT = {
        "ordo": "Orphanet",
        "orphanet": "Orphanet",
        "omim": "OMIM",
        "efo": "EFO",
        "hp": "HP",
        "mondo": "MONDO",
    }

    if not any(x in uri.lower() for x in URI_DB_TO_DB_DICT.keys()):
        return None
    uri = uri.rstrip("/")
    uri_list = uri.split("/")
    if "identifiers.org" in uri:
        db = uri_list[-2]
        id_ = uri_list[-1]
    elif "omim.org" in uri:
        db = "OMIM"
        id_ = uri_list[-1]
    else:
        db, id_ = uri_list[-1].split("_")
    db = URI_DB_TO_DB_DICT[db.lower()]
    return "{}:{}".format(db, id_)
