
class OntologyUri:
    db_to_uri_dict = {
        "orphanet": "http://www.orpha.net/ORDO/Orphanet_{}",
        "omim": "http://identifiers.org/omim/{}",
        "efo": "http://www.ebi.ac.uk/efo/EFO_{}",
        "mesh": "http://identifiers.org/mesh/{}",
        "medgen": "http://identifiers.org/medgen/{}",
        "hp": "http://purl.obolibrary.org/obo/HP_{}",
        "doid": "http://purl.obolibrary.org/obo/DOID_{}",
        "mondo": "http://purl.obolibrary.org/obo/MONDO_{}",
    }

    def __init__(self, id_, db):
        self.id_ = id_
        self.db = db
        self.uri = self.db_to_uri_dict[self.db.lower()].format(self.id_)

    def __str__(self):
        return self.uri
