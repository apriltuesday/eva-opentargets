import logging
import itertools
import copy
import json
import re
import sys
import os
from collections import defaultdict

import jsonschema

from eva_cttv_pipeline import clinvar_xml_utils
from eva_cttv_pipeline.evidence_string_generation import consequence_type as CT

logger = logging.getLogger(__package__)

# Illegal allele regular expression
ILLEGAL_ALLELE_SEQUENCE = re.compile(r'[^ACGT]')

# Output settings
EVIDENCE_STRINGS_FILE_NAME = 'evidence_strings.json'
EVIDENCE_RECORDS_FILE_NAME = 'evidence_records.tsv'
UNMAPPED_TRAITS_FILE_NAME = 'unmappedTraits.tsv'
UNAVAILABLE_EFO_FILE_NAME = 'unavailableefo.tsv'
NSV_LIST_FILE = 'nsvlist.txt'


class Report:
    """Holds counters and other records of a pipeline run. Includes method to write to output files, and __str__ shows
    the summary of the report. One instance of this class is instantiated in the running of the pipeline."""

    def __init__(self, trait_mappings=None):
        if trait_mappings is None:
            self.trait_mappings = {}
        else:
            self.trait_mappings = copy.deepcopy(trait_mappings)

        self.unrecognised_clin_sigs = set()
        self.ensembl_gene_id_uris = set()
        self.traits = set()
        self.n_unrecognised_allele_origin = defaultdict(int)
        self.nsv_list = []
        self.unmapped_traits = defaultdict(int)
        self.evidence_string_count = 0
        self.counters = self.__get_counters()

    def __str__(self):

        report_strings = [
            str(self.counters["record_counter"]) + ' ClinVar records in total',
            str(self.evidence_string_count) + ' evidence string jsons generated',
            str(self.counters["n_processed_clinvar_records"]) +
            ' ClinVar records generated at least one evidence string',
            str(len(self.unrecognised_clin_sigs)) +
            " Clinical significance string(s) not found " +
            "among those described in ClinVar documentation:",
            str(self.unrecognised_clin_sigs),
            str(self.counters["n_same_ref_alt"]) +
            ' ClinVar records with allowed clinical significance ' +
            'did present the same reference and alternate and were skipped',
            'Activities of those ClinVar records with ' +
            'unrecognized clinical significances were set to "unknown".',
            str(len(self.ensembl_gene_id_uris)) +
            ' distinct ensembl gene ids appear in generated evidence string json objects',
            str(len(self.traits)) +
            ' distinct trait names found to include in generated evidence string json objects',
            str(self.counters["n_pathogenic_no_rs"]) +
            ' ClinVar records with allowed clinical significance DO NOT have an rs id',
            str(self.counters["n_multiple_evidence_strings"]) +
            ' ClinVar records generated more than one evidence_string',
            str(self.counters["n_germline_somatic"]) +
            ' ClinVar records with germline and somatic origins',
            str(self.counters["n_multiple_allele_origin"]) +
            ' ClinVar records with more than one allele origin',
            'Number valid ClinVar records with unprocessed allele origins:'
        ]

        report_strings.extend(
            [' ' + alleleOrigin + ': ' +
             str(self.n_unrecognised_allele_origin[alleleOrigin])
             for alleleOrigin in self.n_unrecognised_allele_origin])

        report_strings.extend([
            str(self.counters["n_no_variant_to_ensg_mapping"]) +
            ' ClinVar records with allowed clinical significance and valid rs id ' +
            'were skipped due to a lack of Variant->ENSG mapping.',
            str(self.counters["n_missed_strings_unmapped_traits"]) +
            ' ClinVar records with allowed clinical significance, valid rs id and ' +
            'Variant->ENSG mapping were skipped due to a lack of EFO mapping (see ' +
            UNMAPPED_TRAITS_FILE_NAME + ').',
            str(self.counters["n_records_no_recognised_allele_origin"]) +
            ' ClinVar records with allowed clinical significance, ' +
            'valid rs id, valid Variant->ENSG' +
            ' mapping and valid EFO mapping were skipped due to a lack of a valid alleleOrigin.',
            str(self.counters["n_more_than_one_efo_term"]) +
            ' evidence strings with more than one trait mapped to EFO terms',
            str(self.counters["n_valid_rs_and_nsv"]) +
            ' evidence strings were generated from ClinVar records with rs and nsv ids',
            str(self.counters["n_nsvs"]) + ' total nsvs found',
            str(self.counters["n_nsv_skipped_clin_sig"]) +
            ' ClinVar nsvs were skipped because of a different clinical significance',
            str(self.counters["n_nsv_skipped_wrong_ref_alt"]) +
            ' ClinVar nsvs were skipped because of same ref and alt'
        ])

        return '\n'.join(report_strings)

    def write_output(self, dir_out):
        write_string_list_to_file(self.nsv_list, dir_out + '/' + NSV_LIST_FILE)

        # Contains traits without an ontology mapping
        with open(dir_out + '/' + UNMAPPED_TRAITS_FILE_NAME, 'wt') as fdw:
            fdw.write('Trait\tCount\n')
            for trait_list in self.unmapped_traits:
                fdw.write(str(trait_list) + '\t' +
                          str(self.unmapped_traits[trait_list]) + '\n')

    def remove_trait_mapping(self, trait_name):
        if trait_name in self.trait_mappings:
            del self.trait_mappings[trait_name]

    @staticmethod
    def __get_counters():
        return {"n_processed_clinvar_records": 0,
                "n_pathogenic_no_rs": 0,
                "n_multiple_evidence_strings": 0,
                "n_multiple_allele_origin": 0,
                "n_germline_somatic": 0,
                "n_records_no_recognised_allele_origin": 0,
                "n_no_variant_to_ensg_mapping": 0,
                "n_more_than_one_efo_term": 0,
                "n_same_ref_alt": 0,
                "n_missed_strings_unmapped_traits": 0,
                "n_nsvs": 0,
                "n_valid_rs_and_nsv": 0,
                "n_nsv_skipped_clin_sig": 0,
                "n_nsv_skipped_wrong_ref_alt": 0,
                "record_counter": 0,
                "n_total_clinvar_records": 0}


def validate_evidence_string(ev_string, ot_schema_contents):
    try:
        jsonschema.validate(ev_string, ot_schema_contents, format_checker=jsonschema.FormatChecker())
    except jsonschema.exceptions.ValidationError as err:
        logger.error('Error: evidence string does not validate against schema.')
        logger.error(f'Error message: {err}')
        logger.error(f'Complete evidence string: {json.dumps(ev_string)}')
        sys.exit(1)
    except jsonschema.exceptions.SchemaError:
        logger.error('Error: OpenTargets schema file is invalid')
        sys.exit(1)


def launch_pipeline(clinvar_xml_file, efo_mapping_file, gene_mapping_file, ot_schema_file, dir_out):
    os.makedirs(dir_out, exist_ok=True)
    string_to_efo_mappings = load_efo_mapping(efo_mapping_file)
    variant_to_gene_mappings = CT.process_consequence_type_file(gene_mapping_file)
    report = clinvar_to_evidence_strings(
        string_to_efo_mappings, variant_to_gene_mappings, clinvar_xml_file, ot_schema_file,
        output_evidence_strings=os.path.join(dir_out, EVIDENCE_STRINGS_FILE_NAME))
    report.write_output(dir_out)
    print(report)


def clinvar_to_evidence_strings(string_to_efo_mappings, variant_to_gene_mappings, clinvar_xml, ot_schema,
                                output_evidence_strings):
    report = Report(trait_mappings=string_to_efo_mappings)
    ot_schema_contents = json.loads(open(ot_schema).read())
    output_evidence_strings_file = open(output_evidence_strings, 'wt')

    logger.info('Processing ClinVar records')
    for clinvar_record in clinvar_xml_utils.ClinVarDataset(clinvar_xml):
        report.counters['record_counter'] += 1
        if report.counters['record_counter'] % 1000 == 0:
            logger.info('{} records processed'.format(report.counters['record_counter']))
        n_ev_strings_per_record = 0

        if clinvar_record.measure is None:  # No valid variants to process in this record
            continue
        report.counters['n_nsvs'] += (clinvar_record.measure.nsv_id is not None)
        append_nsv(report.nsv_list, clinvar_record.measure)
        report.counters['n_multiple_allele_origin'] += (len(clinvar_record.allele_origins) > 1)

        # Within each ClinVar record, an evidence string is generated for all possible permutations of:
        # 1. Allele origins
        grouped_allele_origins = convert_allele_origins(clinvar_record.allele_origins)
        # 2. EFO mappings
        grouped_diseases = group_diseases_by_efo_mapping(clinvar_record.traits, string_to_efo_mappings)
        # 3. Genes where the variant has effect
        consequence_types = get_consequence_types(clinvar_record.measure, variant_to_gene_mappings)

        # If at least one of the resulting groups is empty, add this to the report.
        if not grouped_allele_origins:
            report.counters['n_records_no_recognised_allele_origin'] += 1
        if not consequence_types:
            report.counters['n_no_variant_to_ensg_mapping'] += 1
        if not grouped_diseases:
            report.counters['n_missed_strings_unmapped_traits'] += 1
            if not clinvar_record.traits:
                unmapped_trait = 'No traits present in ClinVar record'
            elif not clinvar_record.traits[0].preferred_or_other_name:
                unmapped_trait = 'Trait with no name in ClinVar record'
            else:
                unmapped_trait = clinvar_record.traits[0].preferred_or_other_name
            report.unmapped_traits[unmapped_trait] += 1

        for allele_origins, disease_attributes, consequence_attributes in itertools.product(
                grouped_allele_origins, grouped_diseases, consequence_types):
            disease_name, disease_source_id, disease_mapped_efo_id = disease_attributes
            evidence_string = generate_evidence_string(clinvar_record, allele_origins, disease_name, disease_source_id,
                                                       disease_mapped_efo_id, consequence_attributes)

            # Validate and immediately output the evidence string (not keeping everything in memory)
            validate_evidence_string(evidence_string, ot_schema_contents)
            output_evidence_strings_file.write(json.dumps(evidence_string) + '\n')

            report.evidence_string_count += 1
            report.counters['n_valid_rs_and_nsv'] += (clinvar_record.measure.nsv_id is not None)
            report.traits.add(disease_mapped_efo_id)
            report.remove_trait_mapping(disease_name)
            report.ensembl_gene_id_uris.add(consequence_attributes.ensembl_gene_id)
            n_ev_strings_per_record += 1

        if n_ev_strings_per_record > 0:
            report.counters['n_processed_clinvar_records'] += 1
            if n_ev_strings_per_record > 1:
                report.counters['n_multiple_evidence_strings'] += 1

    output_evidence_strings_file.close()
    return report


def generate_evidence_string(clinvar_record, allele_origins, disease_name, disease_source_id, disease_mapped_efo_id,
                             consequence_attributes):
    """Generates an evidence string based on ClinVar record and some additional attributes."""
    is_somatic = allele_origins == ['somatic']
    evidence_string = {
        # ALLELE ORIGIN ATTRIBUTES. There are three attributes which are currently completely redundant (their
        # values are completely correlated between each other). This is intended, see answer to question 1 here:
        # https://github.com/EBIvariation/eva-opentargets/issues/189#issuecomment-782136128.
        'alleleOrigins': allele_origins,
        'datasourceId': 'eva_somatic' if is_somatic else 'eva',
        'datatypeId': 'somatic_mutation' if is_somatic else 'genetic_association',

        # ASSOCIATION ATTRIBUTES.
        # List of patterns of inheritance reported for the variant.
        'allelicRequirements': clinvar_record.mode_of_inheritance,

        # Levels of clinical significance reported for the variant.
        'clinicalSignificances': clinvar_record.clinical_significance_list,

        # Confidence (review status).
        'confidence': clinvar_record.review_status,

        # Literature. ClinVar records provide three types of references: trait-specific; variant-specific; and
        # "observed in" references. Open Targets are interested only in that last category.
        'literature': sorted(set([str(r) for r in clinvar_record.evidence_support_pubmed_refs])),

        # RCV identifier.
        'studyId': clinvar_record.accession,

        # VARIANT ATTRIBUTES.
        'targetFromSourceId': consequence_attributes.ensembl_gene_id,
        'variantFunctionalConsequenceId': consequence_attributes.so_term.accession,
        'variantId': clinvar_record.measure.vcf_full_coords,  # CHROM_POS_REF_ALT notation
        'variantRsId': clinvar_record.measure.rs_id,

        # PHENOTYPE ATTRIBUTES.
        # The alphabetical list of *all* disease names from that ClinVar record
        'cohortPhenotypes': sorted([trait.preferred_or_other_name for trait in clinvar_record.traits
                                    if trait.preferred_or_other_name is not None]),

        # One disease name for this evidence string (see group_diseases_by_efo_mapping)
        'diseaseFromSource': disease_name,

        # The internal identifier of that disease
        'diseaseFromSourceId': disease_source_id,

        # The EFO identifier to which we mapped that first disease. Converting the URI to a compact representation as
        # required by the Open Targets JSON schema.
        'diseaseFromSourceMappedId': disease_mapped_efo_id.split('/')[-1],
    }
    # Remove the attributes with empty values (either None or empty lists)
    evidence_string = {key: value for key, value in evidence_string.items() if value}
    return evidence_string


def get_consequence_types(clinvar_record_measure, consequence_type_dict):
    """Returns the list of functional consequences for a given ClinVar record measure.

    This is the place where ClinVar records are paired with the information about gene and functional consequences.
    This information is produced by two pipelines in the `vep-mapping-pipeline` subdirectory:
    1. The main one, `vep_mapping_pipeline`, runs all records in the ClinVar VCF dump through Variant Effect Predictor
       and outputs the results using a VCF-compatible "CHROM:POS:REF:ALT" identifier.
    2. The auxiliary one, `repeat_expansion_variants`, uses a different approach to extract information about repeat
       expansion variants from the ClinVar TSV dump. The repeat expansion variants have several important features:
       a. Most (but not all) of them are not present in ClinVar VCF dump, only in the TSV dump.
       b. Even the variants which have the coordinates cannot be adequately processed by VEP: it will output the wrong
          functional consequence type.
       c. The "CHROM:POS:REF:ALT" notation is not useful for these variants, because some of them have an indeterminate
          number of repeats, hence there is no single ALT allele.
       This second pipeline outputs the results using the RCV identifier from ClinVar."""

    # As the first option, try to link a variant and its functional consequence by RCV accession. This will only happen
    # for repeat expansion variants, since the main pipeline uses CHROM:POS:REF:ALT identifiers.
    #
    # The reason we don't use RCV accessions more often is because they are, in general, not variant-specific: the same
    # RCV may link to multiple variants with different functional consequences. However, for repeat expansion variants
    # RCV accessions *are* specific and contain only one variant.
    #
    # The reason we pair first by RCV accession and *then* by CHROM:POS:REF:ALT identifiers is that some repeat
    # expansion variants actually *do* appear in the ClinVar VCF dump, will be fed to VEP, and will receive incorrect
    # consequence annotations. By using RCV pairing first, we prioritise results of the variant expansion pipeline over
    # the general VEP pipeline.
    if clinvar_record_measure.clinvar_record.accession in consequence_type_dict:
        return consequence_type_dict[clinvar_record_measure.clinvar_record.accession]

    # If RCV is not present in the consequences file, pair using full variant description (CHROM:POS:REF:ALT)
    if clinvar_record_measure.has_complete_coordinates:
        # This VCF-flavoured identifier is used to pair ClinVar records with functional consequence predictions.
        # Example of such an identifier: 14:23423715:G:A
        coord_id = ':'.join([clinvar_record_measure.chr, str(clinvar_record_measure.vcf_pos),
                             clinvar_record_measure.vcf_ref, clinvar_record_measure.vcf_alt])
        # Log unusual variants with IUPAC ambiguity symbols. Example: 12_32625716_G_H (from RCV000032000)
        if ILLEGAL_ALLELE_SEQUENCE.search(clinvar_record_measure.vcf_ref + clinvar_record_measure.vcf_alt):
            logger.warning(f'Observed variant with non-ACGT allele sequences: {coord_id}')
        if coord_id in consequence_type_dict:
            return consequence_type_dict[coord_id]

    # Previously, the pairing was also attempted based on rsID and nsvID. This is not reliable because of lack of allele
    # specificity, and has been removed.
    return []


def write_string_list_to_file(string_list, filename):
    with open(filename, 'wt') as out_file:
        out_file.write('\n'.join(string_list))


def append_nsv(nsv_list, clinvar_record_measure):
    nsv = clinvar_record_measure.nsv_id
    if nsv is not None:
        nsv_list.append(nsv)
    return nsv_list


def load_efo_mapping(efo_mapping_file):
    trait_2_efo = defaultdict(list)
    n_efo_mappings = 0

    with open(efo_mapping_file, 'rt') as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('#') or not line:
                continue
            line_list = line.split('\t')
            assert len(line_list) == 3, f'Incorrect string to EFO mapping format for line {line}'
            clinvar_name, ontology_id, ontology_label = line_list
            trait_2_efo[clinvar_name.lower()].append((ontology_id, ontology_label))
            n_efo_mappings += 1
    logger.info('{} EFO mappings loaded'.format(n_efo_mappings))
    return trait_2_efo


def get_terms_from_file(terms_file_path):
    if terms_file_path is not None:
        print('Loading list of terms...')
        with open(terms_file_path, 'rt') as terms_file:
            terms_list = [line.rstrip() for line in terms_file]
        print(str(len(terms_file_path)) + ' terms found at ' + terms_file_path)
    else:
        terms_list = []

    return terms_list


def convert_allele_origins(orig_allele_origins):
    """Splits the original list of allele origins from ClinVar into up to two groups: one for 'somatic' (if present),
    and another one for all other types, which are considered 'germline' (if present)."""
    orig_allele_origins = {item.lower() for item in orig_allele_origins}
    converted_allele_origins = []
    if 'somatic' in orig_allele_origins:
        converted_allele_origins.append(['somatic'])
        orig_allele_origins.remove('somatic')
    if orig_allele_origins:  # Is there something remaining for the second (non-somatic) group?
        converted_allele_origins.append(sorted(orig_allele_origins))
    return converted_allele_origins


def group_diseases_by_efo_mapping(clinvar_record_traits, string_to_efo_mappings):
    """Processes and groups diseases from a ClinVar record in two ways. First, diseases mapping to the same EFO term are
    grouped together. Second, if a group of diseases has multiple EFO mappings, they are split into separate groups.
    For each group, the tuple of the following values is returned: (1) the original name of the disease which comes
    first lexicographically; (2) the original MedGen ID of that disease; (3) the mapped EFO term of that disease.
    Diseases without mappings are skipped.

    For example, for a ClinVar record with the following mappings:
        * Diseases A, B, C -> EFO_1
        * Disease D -> EFO_2 & EFO_3
        * Diseases E, F -> EFO_4 & EFO_5
        * Disease G -> No mapping

    The following tuples will be generated:
        * (A, MedGen_A, EFO_1)
        * (D, MedGen_D, EFO_2)
        * (D, MedGen_D, EFO_3)
        * (E, MedGen_E, EFO_4)
        * (E, MedGen_E, EFO_5)"""

    # Group traits by their EFO mappings and explode multiple mappings
    efo_to_traits = defaultdict(list)  # Key: EFO ID, value: list of traits mapped to that ID
    for trait in clinvar_record_traits:
        # Try to match using all trait names.
        for trait_name in trait.all_names:
            for efo_id, efo_label in string_to_efo_mappings.get(trait_name.lower(), []):
                efo_to_traits[efo_id].append(trait)

    # Generate tuples by keeping only one disease from each group
    grouped_tuples = []
    for efo_id, traits in efo_to_traits.items():
        traits = sorted(traits, key=lambda t: t.preferred_or_other_name)
        selected_trait = traits[0]
        grouped_tuples.append((selected_trait.preferred_or_other_name, selected_trait.medgen_id, efo_id))
    return grouped_tuples
