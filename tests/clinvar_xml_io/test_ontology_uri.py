from cmat.clinvar_xml_io.ontology_uri import OntologyUri


def test_uri_to_curie():
    assert OntologyUri.uri_to_curie('http://www.orpha.net/ORDO/Orphanet_713') == 'Orphanet:713'
    assert OntologyUri.uri_to_curie('http://purl.obolibrary.org/obo/HP_0002930') == 'HP:0002930'
    assert OntologyUri.uri_to_curie('https://omim.org/entry/300377') == 'OMIM:300377'
    # If for some reason we're already passed a CURIE (probably a mistake in OLS response), return it as-is
    assert OntologyUri.uri_to_curie('HP:0002505') == 'HP:0002505'
    # Not a supported db
    assert OntologyUri.uri_to_curie('http://purl.obolibrary.org/obo/DOID_10652') == None
