from datetime import datetime
from collections import UserDict


class ClinvarRecord(UserDict):
    """
    Class of which instances hold data on individual clinvar records. Subclass of UserDict rather
    than dict in order to use attributes
    """

    # A score for the review status of the assigned clinical significance ranges from 0 to 4 and corresponds to the
    # number of gold stars displayed on ClinVar website. See details here:
    # https://www.ncbi.nlm.nih.gov/clinvar/docs/details/#review_status
    score_map = {
        "CRITERIA_PROVIDED_SINGLE_SUBMITTER": 1,
        "CRITERIA_PROVIDED_CONFLICTING_INTERPRETATIONS": 1,
        "CRITERIA_PROVIDED_MULTIPLE_SUBMITTERS_NO_CONFLICTS": 2,
        "REVIEWED_BY_EXPERT_PANEL": 3,
        "PRACTICE_GUIDELINE": 4,
    }

    def __init__(self, cellbase_dict):
        """Initialise a ClinVar record object from JSON data. See /clinvar-variant-types/README.md for the in-depth
        explanation of ClinVar data model. See also issue #127 for the most recent discussions on changing support of
        different ClinVar record types.
        """
        UserDict.__init__(self, cellbase_dict)
        if 'measureSet' in self.data['referenceClinVarAssertion']:
            # MeasureSet provides information on a variant or a set of variants located on the same chromosomal copy.
            if self.data['referenceClinVarAssertion']['measureSet']['type'] == 'Variant':
                # The measure "list" actually only contains a single variant. This is the only case we are currently
                # supporting. As of July 2020, it accounts for >99.7% of all ClinVar records.
                measure_list = self.data['referenceClinVarAssertion']['measureSet']['measure']
            else:
                # Uncommon record types, such as "Haplotype", "Phase unknown", or "Distinct chromosomes".
                # Not currently supported.
                measure_list = []
        elif 'measureSet' in self.data['referenceClinVarAssertion']['genotypeSet']:
            # The record contains a GenotypeSet, a rare subtype which contains an assertion about a group of variants
            # from several chromosome copies. This could be either a CompoundHeterozygote or a Diplotype, and those
            # types are currently not processed.
            measure_list = []
        else:
            raise KeyError('ClinVar record contains neither a MeasureSet, nor a GenotypeSet')

        self.measures = [ClinvarRecordMeasure(measure_dict, self) for measure_dict in measure_list]

    @property
    def date(self):
        return datetime.utcfromtimestamp(
            self.data['referenceClinVarAssertion']['dateLastUpdated'] / 1000).isoformat()

    @property
    def score(self):
        """Returns a score for the review status of the assigned clinical significance. See score_map above. It should
        be noted that currently this property is not used, but this might change in the future."""
        return self.score_map.get(self.data['referenceClinVarAssertion']['clinicalSignificance']['reviewStatus'], 0)

    @property
    def accession(self):
        return self.data['referenceClinVarAssertion']['clinVarAccession']['acc']

    @property
    def traits(self):
        trait_list = []
        for trait in self.data['referenceClinVarAssertion']['traitSet']['trait']:
            trait_list.append([])
            for name in trait['name']:
                # First trait name in the list will always be the "Preferred" one
                if name['elementValue']['type'] == 'Preferred':
                    trait_list[-1] = [name['elementValue']['value']] + trait_list[-1]
                elif name['elementValue']['type'] in ["EFO URL", "EFO id", "EFO name"]:
                    continue  # if the trait name not originally from clinvar
                else:
                    trait_list[-1].append(name['elementValue']['value'])

        return trait_list

    @property
    def trait_pubmed_refs(self):
        pubmed_refs_list = []
        for trait in self.data['referenceClinVarAssertion']['traitSet']['trait']:
            pubmed_refs_list.append([])
            if 'citation' in trait:
                for citation in trait['citation']:
                    if ('id' in citation) and citation['id'] is not None:
                        for citation_id in citation['id']:
                            if citation_id['source'] == 'PubMed':
                                pubmed_refs_list[-1].append(int(citation_id['value']))

        return pubmed_refs_list

    @property
    def observed_pubmed_refs(self):
        pubmed_refs_list = []
        if 'observedIn' in self.data['referenceClinVarAssertion']:
            for observed_in in self.data['referenceClinVarAssertion']['observedIn']:
                for observed_data in observed_in['observedData']:
                    if 'citation' in observed_data:
                        for citation in observed_data['citation']:
                            if ('id' in citation) and citation['id'] is not None:
                                for citation_id in citation['id']:
                                    if citation_id['source'] == 'PubMed':
                                        pubmed_refs_list.append(int(citation_id['value']))
        return pubmed_refs_list

    @property
    def trait_refs_list(self):
        return [['http://europepmc.org/abstract/MED/' + str(ref) for ref in ref_list]
                for ref_list in self.trait_pubmed_refs]

    @property
    def observed_refs_list(self):
        return ['http://europepmc.org/abstract/MED/' + str(ref)
                for ref in self.observed_pubmed_refs]

    @property
    def clinical_significance(self):
        return \
            self.data['referenceClinVarAssertion']['clinicalSignificance']['description']

    @property
    def allele_origins(self):
        allele_origins = set()
        for clinvar_assertion_document in self.data['clinVarAssertion']:
            for observed_in_document in clinvar_assertion_document['observedIn']:
                allele_origins.add(observed_in_document['sample']['origin'].lower())

        return list(allele_origins)


class ClinvarRecordMeasure(UserDict):

    def __init__(self, clinvar_measure_dict, clinvar_record):
        UserDict.__init__(self, clinvar_measure_dict)
        self.clinvar_record = clinvar_record

    @property
    def rs_id(self):
        if "xref" in self.data:
            for xref in self.data["xref"]:
                if xref["db"].lower() == "dbsnp":
                    return "rs{}".format(xref["id"])
        return None

    @property
    def nsv_id(self):
        if "xref" in self.data:
            for xref in self.data["xref"]:
                if xref["db"].lower() == "dbvar" and xref["id"].lower()[:3] in ("nsv", "esv"):
                    return xref["id"]
        return None

    @property
    def hgvs(self):
        hgvs_list = []
        for attribute_set in self.data['attributeSet']:
            if attribute_set['attribute']['type'].startswith('HGVS'):
                hgvs_list.append(attribute_set['attribute']['value'])

        return hgvs_list

    @property
    def variant_type(self):
        return self.data['type']

    @property
    def pubmed_refs(self):
        pubmed_refs_list = []
        if 'citation' in self.data:
            for citation in self.data['citation']:
                if 'id' in citation and citation['id'] is not None:
                    for citation_id in citation['id']:
                        if citation_id['source'] == 'PubMed':
                            pubmed_refs_list.append(int(citation_id['value']))
        return pubmed_refs_list

    @property
    def refs_list(self):
        return ['http://europepmc.org/abstract/MED/' + str(ref)
                for ref in self.pubmed_refs]

    @property
    def chr(self):
        return self.sequence_location_helper("chr")

    @property
    def vcf_pos(self):
        return self.sequence_location_helper("positionVCF")

    @property
    def vcf_ref(self):
        return self.sequence_location_helper("referenceAlleleVCF")

    @property
    def vcf_alt(self):
        return self.sequence_location_helper("alternateAlleleVCF")

    @property
    def has_complete_coordinates(self):
        return self.chr and self.vcf_pos and self.vcf_ref and self.vcf_alt

    def sequence_location_helper(self, attr):
        if "sequenceLocation" in self.data:
            for sequence_location in self.data["sequenceLocation"]:
                if sequence_location["assembly"].lower() == "grch38":
                    if attr in sequence_location:
                        return sequence_location[attr]
        return None

