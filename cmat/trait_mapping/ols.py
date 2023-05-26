from functools import lru_cache
import logging
import requests
import urllib

from cmat.trait_mapping.utils import json_request


OLS_EFO_SERVER = 'https://www.ebi.ac.uk/ols'
# The setting for local OLS installation should be uncommented if necessary. Note that the link
# for the local deployment is different from the production link in three regards: (1) it must use
# HTTP instead of HTTPS; (2) it must include the port which you used when deploying the Docker
# container; (3) it does *not* include /ols in its path.
# OLS_EFO_SERVER = 'http://127.0.0.1:8080'

logger = logging.getLogger(__package__)


def build_ols_query(ontology_uri: str) -> str:
    """Build a url to query OLS for a given ontology uri."""
    return "https://www.ebi.ac.uk/ols/api/terms?iri={}".format(ontology_uri)


@lru_cache(maxsize=16384)
def get_ontology_label_from_ols(ontology_uri: str) -> str:
    """
    Using provided ontology URI, build an OLS URL with which to make a request to find the term label for this URI.

    :param ontology_uri: A URI for a term in an ontology.
    :return: Term label for the ontology URI provided in the parameters.
    """
    url = build_ols_query(ontology_uri)
    json_response = json_request(url)

    if not json_response:
        return ''

    # If the '_embedded' section is missing from the response, it means that the term is not found in OLS
    if '_embedded' not in json_response:
        if '/medgen/' not in url and '/omim/' not in url:
            logger.warning('OLS queried OK but did not return any results for URL {}'.format(url))
        return ''

    # Go through all terms found by the requested identifier and try to find the one where the _identifier_ and the
    # _term_ come from the same ontology (marked by a special flag). Example of such a situation would be a MONDO term
    # in the MONDO ontology. Example of a reverse situation is a MONDO term in EFO ontology (being *imported* into it
    # at some point).
    for term in json_response["_embedded"]["terms"]:
        if term["is_defining_ontology"]:
            return term["label"]

    if '/medgen/' not in url and '/omim/' not in url:
        logger.warning('OLS queried OK, but there is no defining ontology in its results for URL {}'.format(url))
    return ''


def double_encode_uri(uri: str) -> str:
    """Double encode a given uri."""
    return urllib.parse.quote(urllib.parse.quote(uri, safe=""), safe="")


def ols_efo_query(uri: str) -> requests.Response:
    """
    Query EFO using OLS for a given ontology uri, returning the response from the request.

    :param uri: Ontology uri to use in querying EFO using OLS
    :return: Response from OLS
    """
    double_encoded_uri = double_encode_uri(uri)
    return requests.get(
        "{}/api/ontologies/efo/terms/{}".format(OLS_EFO_SERVER, double_encoded_uri))


@lru_cache(maxsize=16384)
def is_current_and_in_efo(uri: str) -> bool:
    """
    Checks whether given ontology uri is a valid and non-obsolete term in EFO.

    :param uri: Ontology uri to use in querying EFO using OLS
    :return: Boolean value, true if ontology uri is valid and non-obsolete term in EFO
    """
    response = ols_efo_query(uri)
    if response.status_code != 200:
        return False
    response_json = response.json()
    return not response_json["is_obsolete"]


@lru_cache(maxsize=16384)
def is_in_efo(uri: str) -> bool:
    """
    Checks whether given ontology uri is a valid term in EFO.

    :param uri: Ontology uri to use in querying EFO using OLS
    :return: Boolean value, true if ontology uri is valid term in EFO
    """
    response = ols_efo_query(uri)
    return response.status_code == 200


@lru_cache(maxsize=16384)
def get_replacement_term(uri: str) -> str:
    """
    Finds replacement term in EFO (if present) for the given ontology uri.

    :param uri: Ontology uri to use in querying EFO using OLS
    :return: Replacement term URI or empty string if not obsolete
    """
    response = ols_efo_query(uri)
    if response.status_code != 200:
        return ""
    response_json = response.json()
    if response_json["term_replaced_by"] is not None:
        return response_json["term_replaced_by"]
    return ""
