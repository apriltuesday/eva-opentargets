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
