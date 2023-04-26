import logging

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class OntologyUri:
    # ClinVar stores cross-references in very different formats. This provides their conversion to full IRIs, along with
    # some examples of how this looks like in ClinVar data.
    db_to_uri_conversion = {
        'orphanet': lambda x: f'http://www.orpha.net/ORDO/Orphanet_{x}',  # <XRef ID="1756" DB="Orphanet"/>
        'omim': lambda x: f'https://www.omim.org/entry/{x}',  # <XRef Type="MIM" ID="612773" DB="OMIM"/>
        'efo': lambda x: f'http://www.ebi.ac.uk/efo/{x}',  # <XRef ID="EFO_0005137" DB="EFO"/>
        'mesh': lambda x: f'http://identifiers.org/mesh/{x}',  # <XRef ID="D065630" DB="MeSH"/>
        'medgen': lambda x: f'http://identifiers.org/medgen/{x}',  # <XRef ID="C0235833" DB="MedGen"/>
        # <XRef ID="MONDO:0013353" DB="MONDO"/>
        'mondo': lambda x: 'http://purl.obolibrary.org/obo/{}'.format(x.replace(':', '_')),
        # <XRef ID="HP:0011147" DB="Human Phenotype Ontology"/>
        'hp': lambda x: 'http://purl.obolibrary.org/obo/{}'.format(x.replace(':', '_')),
    }

    def __init__(self, id_, db):
        self.id_ = id_
        self.db = db if db.lower() != 'human phenotype ontology' else 'HP'
        self.uri = self.db_to_uri_conversion[self.db.lower()](self.id_)

    def __str__(self):
        return self.uri

    @property
    def curie(self):
        return self.uri_to_curie(self.uri)

    @staticmethod
    def uri_to_curie(uri):
        """Convert an ontology uri to a DB:ID format."""
        uri_db_to_curie_db = {
            "ordo": "Orphanet",
            "orphanet": "Orphanet",
            "omim": "OMIM",
            "efo": "EFO",
            "hp": "HP",
            "mondo": "MONDO",
        }
        if not any(x in uri.lower() for x in uri_db_to_curie_db.keys()):
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
            last_component = uri_list[-1]
            if ":" in last_component:
                return last_component
            elif "_" in last_component:
                db, id_ = last_component.split("_")
            else:
                logger.warning(f"Could not convert URI to CURIE: {uri}")
                return None
        db = uri_db_to_curie_db[db.lower()]
        return "{}:{}".format(db, id_)
