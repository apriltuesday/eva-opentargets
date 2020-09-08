import copy
import json
import re

import jsonschema

from eva_cttv_pipeline.evidence_string_generation import config
from eva_cttv_pipeline.evidence_string_generation import utilities


def get_cttv_variant_type(clinvar_record_measure):
    if clinvar_record_measure.vcf_ref is not None and clinvar_record_measure.vcf_alt is not None:
        if len(clinvar_record_measure.vcf_ref) < 2 and len(clinvar_record_measure.vcf_alt) < 2:
            cttv_variant_type = "snp single"
        elif len(clinvar_record_measure.vcf_ref) > 50 or len(clinvar_record_measure.vcf_alt) > 50:
            cttv_variant_type = "structural variant"
        else:
            cttv_variant_type = "snp single"  # Sam asked for this in his email 21/05/2015
            # cttv_variant_type = 'snp multiple'
    else:
        if clinvar_record_measure.rs_id is not None:
            cttv_variant_type = "snp single"
        elif clinvar_record_measure.nsv_id is not None:
            cttv_variant_type = "structural variant"
        else:
            cttv_variant_type = "snp single"

    return cttv_variant_type


def process_clinical_significance(clin_sig):
    """
    Processes ClinVar clinical significance string into a format suitable for OT JSON schema.

    Namely, splits multiple clinical significance levels into an array and normalises names (to lowercase, using only
    spaces for delimiters). Multiple levels of clinical significance are separated using two delimieters: ('/', ', ').
    See /clinvar-variant-types/README.md for further explanation. The output array is sorted alphabetically.

    Example: 'Benign/Likely benign, risk_factor' → ['benign', 'likely benign', 'risk factor'].
    """
    return sorted(re.split('/|, ', clin_sig.lower().replace('_', ' ')))


class CTTVEvidenceString(dict):
    """
    Base evidence string class. Holds variables and methods common between somatic and genetic
    evidence strings.
    Subclass of dict to use indexing.
    """

    def __init__(self, a_dictionary, clinvar_record=None, ref_list=None,
                 ensembl_gene_id=None, report=None, trait=None):
        super().__init__(a_dictionary)

        if ensembl_gene_id:
            self.add_unique_association_field('gene', ensembl_gene_id)
        if clinvar_record:
            self.add_unique_association_field('clinvarAccession', clinvar_record.accession)

        if ensembl_gene_id:
            ensembl_gene_id_uri = get_ensembl_gene_id_uri(ensembl_gene_id)
            self.set_target(ensembl_gene_id_uri)

        if ref_list and len(ref_list) > 0:
            self.top_level_literature = ref_list

        self.disease_id = trait.ontology_id
        self.add_unique_association_field('phenotype', trait.ontology_id)

        if trait.clinvar_name:
            self.disease_source_name = trait.clinvar_name

        if trait.ontology_label:
            self.disease_name = trait.ontology_label

    def add_unique_association_field(self, key, value):
        self['unique_association_fields'][key] = value

    def _clear_target(self):
        self['target']['id'] = []

    def set_target(self, target_id):
        self['target']['id'] = target_id

    @property
    def disease_name(self):
        return self['disease']['name']

    @disease_name.setter
    def disease_name(self, value):
        self['disease']['name'] = value

    @property
    def disease_source_name(self):
        return self['disease']['source_name']

    @disease_source_name.setter
    def disease_source_name(self, value):
        self['disease']['source_name'] = value

    @property
    def disease_id(self):
        return self['disease']['id']

    @disease_id.setter
    def disease_id(self, value):
        self['disease']['id'] = value

    @property
    def evidence_codes(self):
        return self['evidence']['evidence_codes']

    @evidence_codes.setter
    def evidence_codes(self, ev_code_list):
        self['evidence']['evidence_codes'] = ev_code_list

    @property
    def top_level_literature(self):
        return self['literature']['references']

    @top_level_literature.setter
    def top_level_literature(self, reference_list):
        self['literature'] = \
            {'references': [{'lit_id': reference} for reference in reference_list]}

    def validate(self, ot_schema_contents):
        jsonschema.validate(self, ot_schema_contents, format_checker=jsonschema.FormatChecker())
        return True


class CTTVGeneticsEvidenceString(CTTVEvidenceString):
    """
    Class for genetics evidence string specifically.
    Holds information required for Open Target's evidence strings for genetic information.
    """

    with utilities.open_file(utilities.get_resource_file(__package__, config.GEN_EV_STRING_JSON),
                             "rt") as gen_json_file:
        base_json = json.load(gen_json_file)

    def __init__(self, clinvar_record, clinvar_record_measure, report, trait, consequence_type):

        a_dictionary = copy.deepcopy(self.base_json)

        ref_list = list(set(clinvar_record.trait_refs_list[trait.trait_counter] +
                            clinvar_record.observed_refs_list +
                            clinvar_record_measure.refs_list))
        ref_list.sort()

        super().__init__(a_dictionary, clinvar_record, ref_list, consequence_type.ensembl_gene_id, report, trait)

        variant_type = get_cttv_variant_type(clinvar_record_measure)

        self.add_unique_association_field('alleleOrigin', 'germline')
        if clinvar_record_measure.rs_id:
            self.set_variant('http://identifiers.org/dbsnp/' + clinvar_record_measure.rs_id,
                             variant_type)
            self.add_unique_association_field('variant_id', clinvar_record_measure.rs_id)
        elif clinvar_record_measure.nsv_id:
            self.set_variant('http://identifiers.org/dbsnp/' + clinvar_record_measure.nsv_id,
                             variant_type)
            self.add_unique_association_field('variant_id', clinvar_record_measure.nsv_id)
        else:
            self.set_variant('http://www.ncbi.nlm.nih.gov/clinvar/' + clinvar_record.accession,
                             variant_type)
            self.add_unique_association_field('variant_id', clinvar_record.accession)
        self.date = clinvar_record.date
        self.db_xref_url = 'http://identifiers.org/clinvar.record/' + clinvar_record.accession
        self.url = 'http://www.ncbi.nlm.nih.gov/clinvar/' + clinvar_record.accession
        # See https://github.com/opentargets/platform/issues/1139#issuecomment-682592678
        self.association = True
        self.gene_2_var_ev_codes = ['http://identifiers.org/eco/cttv_mapping_pipeline']
        most_severe_so_term = consequence_type.so_term
        if most_severe_so_term.accession is None:
            self.gene_2_var_func_consequence = 'http://targetvalidation.org/sequence/' + \
                                               most_severe_so_term.so_name
        else:
            self.gene_2_var_func_consequence = 'http://purl.obolibrary.org/obo/' + \
                                               most_severe_so_term.accession.replace(':', '_')

        if len(ref_list) > 0:
            self.set_var_2_disease_literature(ref_list)
            # Arbitrarily select only one reference among all
            self.unique_reference = ref_list[0]

        if clinvar_record.clinical_significance:
            self.clinical_significance = process_clinical_significance(clinvar_record.clinical_significance)

        # Populate star rating and review status
        star_rating, review_status = clinvar_record.score
        self.clinvar_rating = (star_rating, review_status)

        # Populate mode of inheritance (if present)
        self.mode_of_inheritance = clinvar_record.mode_of_inheritance

    @property
    def db_xref_url(self):
        if self['evidence']['gene2variant']['provenance_type']['database']['dbxref']['url'] \
                == self['evidence']['variant2disease']['provenance_type']['database']['dbxref']['url']:
            return \
                self['evidence']['variant2disease']['provenance_type']['database']['dbxref']['url']
        else:
            raise Exception("db_xref_url attributes different")

    @db_xref_url.setter
    def db_xref_url(self, url):
        self['evidence']['gene2variant']['provenance_type']['database']['dbxref']['url'] = url
        self['evidence']['variant2disease']['provenance_type']['database']['dbxref']['url'] = url

    @property
    def url(self):
        if self['evidence']['gene2variant']['urls'][0]['url'] \
                == self['evidence']['variant2disease']['urls'][0]['url']:
            return self['evidence']['gene2variant']['urls'][0]['url']
        else:
            raise Exception("url attributes different")

    @url.setter
    def url(self, url):
        self['evidence']['gene2variant']['urls'][0]['url'] = url
        self['evidence']['variant2disease']['urls'][0]['url'] = url

    @property
    def gene_2_var_ev_codes(self):
        return self['evidence']['gene2variant']['evidence_codes']

    @gene_2_var_ev_codes.setter
    def gene_2_var_ev_codes(self, gene_2_var_ev_codes):
        self['evidence']['gene2variant']['evidence_codes'] = gene_2_var_ev_codes

    @property
    def gene_2_var_func_consequence(self):
        return self['evidence']['gene2variant']['functional_consequence']

    @gene_2_var_func_consequence.setter
    def gene_2_var_func_consequence(self, so_term):
        self['evidence']['gene2variant']['functional_consequence'] = so_term

    def set_var_2_disease_literature(self, ref_list):
        self['evidence']['variant2disease']['provenance_type']['literature'] = \
            {'references': [{'lit_id': reference} for reference in ref_list]}

    @property
    def association(self):
        if self['evidence']['gene2variant']['is_associated'] \
                == self['evidence']['variant2disease']['is_associated']:
            return self['evidence']['gene2variant']['is_associated']
        else:
            raise Exception("association attributes different")

    @association.setter
    def association(self, is_associated):
        self['evidence']['gene2variant']['is_associated'] = is_associated
        self['evidence']['variant2disease']['is_associated'] = is_associated

    def _clear_variant(self):
        self['variant']['id'] = []
        self['variant']['type'] = []

    def set_variant(self, var_id, var_type):
        self['variant']['id'] = var_id
        self['variant']['type'] = var_type

    @property
    def unique_reference(self):
        return self['evidence']['variant2disease']['unique_experiment_reference']

    @unique_reference.setter
    def unique_reference(self, reference):
        self['evidence']['variant2disease']['unique_experiment_reference'] = reference

    @property
    def date(self):
        if self['evidence']['gene2variant']['date_asserted'] == \
                self['evidence']['variant2disease']['date_asserted']:
            return self['evidence']['gene2variant']['date_asserted']
        else:
            raise Exception("date attributes have different values")

    @date.setter
    def date(self, date_string):
        self['evidence']['gene2variant']['date_asserted'] = date_string
        self['evidence']['variant2disease']['date_asserted'] = date_string

    @property
    def clinical_significance(self):
        return self['evidence']['variant2disease']['clinical_significance']

    @clinical_significance.setter
    def clinical_significance(self, clinical_significance):
        self['evidence']['variant2disease']['clinical_significance'] = clinical_significance

    @property
    def clinvar_rating(self):
        return self['evidence']['variant2disease']['clinvar_rating']

    @clinvar_rating.setter
    def clinvar_rating(self, clinvar_rating_data):
        star_rating, review_status = clinvar_rating_data
        self['evidence']['variant2disease']['clinvar_rating'] = {
            'star_rating': star_rating,
            'review_status': review_status,
        }

    @property
    def mode_of_inheritance(self):
        return self['evidence']['variant2disease'].get('mode_of_inheritance')

    @mode_of_inheritance.setter
    def mode_of_inheritance(self, mode_of_inheritance):
        if mode_of_inheritance:
            self['evidence']['variant2disease']['mode_of_inheritance'] = mode_of_inheritance


class CTTVSomaticEvidenceString(CTTVEvidenceString):

    """
    Class for somatic evidence string specifically.
    Holds information required for Open Target's evidence strings for somatic information.
    """

    with utilities.open_file(utilities.get_resource_file(__package__, config.SOM_EV_STRING_JSON),
                             "rt") as som_json_file:
        base_json = json.load(som_json_file)

    def __init__(self, clinvar_record, clinvar_record_measure, report, trait, consequence_type):

        a_dictionary = copy.deepcopy(self.base_json)

        ref_list = list(set(clinvar_record.trait_refs_list[trait.trait_counter] +
                            clinvar_record.observed_refs_list +
                            clinvar_record_measure.refs_list))
        ref_list.sort()

        super().__init__(a_dictionary, clinvar_record, ref_list, consequence_type.ensembl_gene_id, report, trait)

        self.add_unique_association_field('alleleOrigin', 'somatic')
        if clinvar_record_measure.rs_id:
            self.add_unique_association_field('variant_id', clinvar_record_measure.rs_id)
        elif clinvar_record_measure.nsv_id:
            self.add_unique_association_field('variant_id', clinvar_record_measure.nsv_id)
        else:
            self.add_unique_association_field('variant_id', clinvar_record.accession)

        self.date = clinvar_record.date
        self.db_xref_url = 'http://identifiers.org/clinvar.record/' + clinvar_record.accession
        self.url = 'http://www.ncbi.nlm.nih.gov/clinvar/' + clinvar_record.accession
        # See https://github.com/opentargets/platform/issues/1139#issuecomment-682592678
        self.association = True

        self.set_known_mutations(consequence_type.so_term)

        if len(ref_list) > 0:
            self.evidence_literature = ref_list

        if clinvar_record.clinical_significance:
            self.clinical_significance = process_clinical_significance(clinvar_record.clinical_significance)

        # Populate star rating and review status
        star_rating, review_status = clinvar_record.score
        self.clinvar_rating = (star_rating, review_status)

    @property
    def db_xref_url(self):
        return self['evidence']['provenance_type']['database']['dbxref']['url']

    @db_xref_url.setter
    def db_xref_url(self, url):
        self['evidence']['provenance_type']['database']['dbxref']['url'] = url

    @property
    def url(self):
        return self['evidence']['urls'][0]['url']

    @url.setter
    def url(self, url):
        self['evidence']['urls'][0]['url'] = url

    @property
    def evidence_literature(self):
        return self['evidence']['provenance_type']['literature']['references']

    @evidence_literature.setter
    def evidence_literature(self, ref_list):
        self['evidence']['provenance_type']['literature'] = \
            {'references': [{'lit_id': reference} for reference in ref_list]}

    @property
    def association(self):
        return self['evidence']['is_associated']

    @association.setter
    def association(self, is_associated):
        self['evidence']['is_associated'] = is_associated

    @property
    def date(self):
        return self['evidence']['date_asserted']

    @date.setter
    def date(self, date_string):
        self['evidence']['date_asserted'] = date_string

    def _clear_known_mutations(self):
        self['evidence']['known_mutations'] = []

    def add_known_mutation(self, new_functional_consequence, so_name):
        new_known_mutation = \
            {'functional_consequence': new_functional_consequence, 'preferred_name': so_name}
        self['evidence']['known_mutations'].append(new_known_mutation)

    def set_known_mutations(self, so_term):
        if so_term.accession:
            new_functional_consequence = \
                "http://purl.obolibrary.org/obo/" + so_term.accession.replace(':', '_')
        else:
            new_functional_consequence = \
                'http://targetvalidation.org/sequence/' + so_term.so_name
        self.add_known_mutation(new_functional_consequence, so_term.so_name)

    @property
    def clinical_significance(self):
        return self['evidence']['clinical_significance']

    @clinical_significance.setter
    def clinical_significance(self, clinical_significance):
        self['evidence']['clinical_significance'] = clinical_significance

    @property
    def clinvar_rating(self):
        return self['evidence']['clinvar_rating']

    @clinvar_rating.setter
    def clinvar_rating(self, clinvar_rating_data):
        star_rating, review_status = clinvar_rating_data
        self['evidence']['clinvar_rating'] = {
            'star_rating': star_rating,
            'review_status': review_status,
        }


def get_ensembl_gene_id_uri(ensembl_gene_id):
    return 'http://identifiers.org/ensembl/' + ensembl_gene_id
