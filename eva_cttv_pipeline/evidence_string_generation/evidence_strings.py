import jsonschema


class EvidenceString(dict):
    """Evidence string class. Subclass of dict to use indexing."""

    def __init__(self, allele_origin, clinvar_record, clinvar_trait, ontology_id, consequence_type):
        super().__init__({})

        # Section A. ALLELE ORIGIN ATTRIBUTES. Differentiate between germline and somatic variants.
        # FIXME 1: this seems redundant.
        # FIXME 2: what to do with RCVs with nonstandard values?
        # FIXME 3: what to do with RCVs with both values?
        if allele_origin == 'germline':
            self.allele_origins = ['germline']
            self.datasource_id = 'eva'
            self.datatype_id = 'genetic_association'
        elif allele_origin == 'somatic':
            self.allele_origins = ['somatic']
            self.datasource_id = 'eva_somatic'
            self.datatype_id = 'somatic_mutation'
        else:
            raise AssertionError(f'Unknown allele origin: {allele_origin}')

        # ASSOCIATION ATTRIBUTES.
        # List of patterns of inheritance reported for the variant.
        self.allelic_requirements = clinvar_record.mode_of_inheritance
        # Levels of clinical significance reported for the variant.
        self.clinical_significances = clinvar_record.clinical_significance_list
        # Confidence (review status).
        self.confidence = clinvar_record.review_status
        # Literature.
        # FIXME: rewrite for passing all traits and all variants
        self.literature = sorted(set(clinvar_trait.pubmed_refs +             # Trait-specific references
                                     clinvar_record.measure.pubmed_refs +    # Variant-specific references
                                     clinvar_record.observed_pubmed_refs))   # "ObservedIn" references
        # RCV identifier.
        self.study_id = clinvar_record.accession

        # VARIANT ATTRIBUTES.
        self.target_from_source_id = consequence_type.ensembl_gene_id
        self.variant_functional_consequence_id = consequence_type.so_term.name
        self.variant_id = clinvar_record.measure.vcf_full_coords  # CHROM_POS_REF_ALT notation
        self.variant_rs_id = clinvar_record.measure.rs_id

        # PHENOTYPE ATTRIBUTES.
        # The alphabetical list of all disease names from ClinVar
        # FIXME: currently only one trait is being passed.
        self.cohort_phenotypes = sorted([trait.name.lower() for trait in clinvar_record.traits])
        # The first item in the list of all disease names
        self.disease_from_source = self.cohort_phenotypes[0]
        # The internal identifier of that first disease
        # FIXME: not currently populated
        self.disease_from_source_id = 'None'
        # The EFO identifier to which we mapped that first disease
        # FIXME: what to do with multiple?
        self.disease_from_source_mapped_id = ontology_id

    # ALLELIC ORIGIN DEFINITIONS
    # alleleOrigins: whether the variant is observed as germline and/or somatic.
    @property
    def allele_origins(self):
        return self.get('alleleOrigins')

    @allele_origins.setter
    def allele_origins(self, value):
        self['alleleOrigins'] = value

    # datasourceId splits all variants into germline and somatic.
    @property
    def datasource_id(self):
        return self.get('datasourceId')

    @datasource_id.setter
    def datasource_id(self, value):
        self['datasourceId'] = value

    # datatypeId also splits all variants into germline and somatic.
    @property
    def datatype_id(self):
        return self['datatypeId']

    @datatype_id.setter
    def datatype_id(self, value):
        self['datatypeId'] = value

    # PHENOTYPE DEFINITIONS
    # cohortPhenotypes: full list of original disease names from ClinVar.
    @property
    def cohort_phenotypes(self):
        return self.get('cohortPhenotypes')

    @cohort_phenotypes.setter
    def cohort_phenotypes(self, value):
        assert value is not None, 'Missing list of cohort phenotypes'
        self['cohortPhenotypes'] = value

    # diseaseFromSource: the first disease name alphabetically.
    @property
    def disease_from_source(self):
        return self['diseaseFromSource']

    @disease_from_source.setter
    def disease_from_source(self, value):
        self['diseaseFromSource'] = value

    # diseaseFromSourceId: the original ontology identifier for the disease in the diseaseFromSource field.
    @property
    def disease_from_source_id(self):
        return self['diseaseFromSourceId']

    @disease_from_source_id.setter
    def disease_from_source_id(self, value):
        self['diseaseFromSourceId'] = value

    # diseaseFromSourceMappedId: which EFO term we mapped the first disease to.
    @property
    def disease_from_source_mapped_id(self):
        return self['diseaseFromSourceMappedId']

    @disease_from_source_mapped_id.setter
    def disease_from_source_mapped_id(self, value):
        self['diseaseFromSourceMappedId'] = value

    # OTHER VARIANT INFORMATION DEFINITIONS
    # allelicRequirements: list of patterns of inheritance observed for the variant.
    @property
    def allelic_requirements(self):
        return self.get('allelicRequirements')

    @allelic_requirements.setter
    def allelic_requirements(self, value):
        if value:
            self['allelicRequirements'] = value

    # clinicalSignificances: list of clinical significance levels reported for the variant.
    @property
    def clinical_significances(self):
        return self.get('clinicalSignificances')

    @clinical_significances.setter
    def clinical_significances(self, value):
        if value:
            self['clinicalSignificances'] = value

    # confidence: free text description.
    @property
    def confidence(self):
        return self['confidence']

    @confidence.setter
    def confidence(self, value):
        self['confidence'] = value

    # literature: PubMed references.
    @property
    def literature(self):
        return self['literature']

    @literature.setter
    def literature(self, value):
        if value:
            self['literature'] = value

    # studyId: RCV identifier.
    @property
    def study_id(self):
        return self['studyId']

    @study_id.setter
    def study_id(self, value):
        self['studyId'] = value

    # TARGET INFORMATION.
    # targetFromSourceId: Ensembl gene ID
    @property
    def target_from_source_id(self):
        return self['targetFromSourceId']

    @target_from_source_id.setter
    def target_from_source_id(self, value):
        self['targetFromSourceId'] = value

    # variantFunctionalConsequenceId
    @property
    def variant_functional_consequence_id(self):
        return self['variantFunctionalConsequenceId']

    @variant_functional_consequence_id.setter
    def variant_functional_consequence_id(self, value):
        self['variantFunctionalConsequenceId'] = value

    # variantId
    @property
    def variant_id(self):
        return self['variantId']

    @variant_id.setter
    def variant_id(self, value):
        self['variantId'] = value

    # variantRsId
    @property
    def variant_rs_id(self):
        return self['variantRsId']

    @variant_rs_id.setter
    def variant_rs_id(self, value):
        self['variantRsId'] = value

    def validate(self, ot_schema_contents):
        jsonschema.validate(self, ot_schema_contents, format_checker=jsonschema.FormatChecker())
        return True
