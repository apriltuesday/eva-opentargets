import xml.etree.ElementTree as ET

from eva_cttv_pipeline.clinvar_xml_io.clinvar_xml_io import ClinVarTrait, ClinVarRecordMeasure, ClinVarDataset, \
    ClinVarRecord
from eva_cttv_pipeline.clinvar_xml_io.clinvar_xml_io.xml_parsing import iterate_rcv_from_xml
from eva_cttv_pipeline.evidence_string_generation.clinvar_to_evidence_strings import load_efo_mapping, \
    get_consequence_types
from eva_cttv_pipeline.evidence_string_generation import consequence_type as CT


PROCESSOR = 'CMAT'


class AnnotatingClinVarDataset(ClinVarDataset):
    """This class provides the ability to parse ClinVar records (RCVs) and annotate them with EFO mappings and
    consequence mappings on the fly."""

    def __init__(self, clinvar_xml, string_to_efo_mappings, variant_to_gene_mappings):
        super().__init__(clinvar_xml)
        self.header_attr['ProcessedBy'] = PROCESSOR
        self.string_to_efo_mappings = string_to_efo_mappings
        self.variant_to_gene_mappings = variant_to_gene_mappings

        self.overall_counts = {}
        self.cmat_counts = {}
        self.clinvar_counts = {}
        self.match_counts = {}
        self.f1_scores = {}

    def __iter__(self):
        # Initialise counts
        self.overall_counts = {
            'total': 0,
            'has_supported_measure': 0,
            'has_supported_trait': 0,
        }
        self.cmat_counts = {
            'has_gene': 0,
            'has_consequences': 0,
            'has_efo_mappings': 0
        }
        self.clinvar_counts = {
            'has_gene': 0,
            'has_consequences': 0,
            'has_any_mappings': 0,
            'has_efo_mappings': 0
        }
        self.match_counts = {
            'genes': 0,
            'conseqs': 0,
            'traits': 0
        }
        self.f1_scores = {
            'genes': 0,
            'conseqs': 0,
            'traits': 0
        }
        for rcv in iterate_rcv_from_xml(self.clinvar_xml):
            record = AnnotatedClinVarRecord(rcv)
            self.annotate(record)
            yield record
        # Compute averages for scores
        for key in self.f1_scores:
            self.f1_scores[key] /= self.match_counts[key]

    def annotate(self, record):
        self.overall_counts['total'] += 1
        # Functional consequences for measure
        if record.measure:
            self.overall_counts['has_supported_measure'] += 1
            self.annotate_and_count_measure(record)
        # EFO terms for traits
        if record.traits_with_valid_names:
            self.overall_counts['has_supported_trait'] += 1
            self.annotate_and_count_traits(record)

    def annotate_and_count_measure(self, record):
        consequence_types = get_consequence_types(record.measure, self.variant_to_gene_mappings)
        record.measure.add_ensembl_annotations(consequence_types)

        annotated_genes = {ct.ensembl_gene_id for ct in consequence_types if ct.ensembl_gene_id}
        annotated_conseqs = {EnsemblAnnotatedClinVarMeasure.format_so_term(ct.so_term)
                             for ct in consequence_types if ct.so_term}
        if annotated_genes:
            self.cmat_counts['has_gene'] += 1
        if annotated_conseqs:
            self.cmat_counts['has_consequences'] += 1

        if record.measure.hgnc_ids:
            self.clinvar_counts['has_gene'] += 1
        if record.measure.so_terms:
            self.clinvar_counts['has_consequences'] += 1

        if record.measure.hgnc_ids and annotated_genes:
            # TODO need to map genes - try using gene symbol / gene name
            self.f1_scores['genes'] += self.f1_score(set(record.measure.hgnc_ids), annotated_genes)
            self.match_counts['genes'] += 1
        if record.measure.so_terms and annotated_conseqs:
            self.f1_scores['conseqs'] += self.f1_score(record.measure.so_terms, annotated_conseqs)
            self.match_counts['conseqs'] += 1

    def annotate_and_count_traits(self, record):
        has_existing_id = False
        existing_efo_ids = set()
        annotated_efo_ids = set()

        for trait in record.traits_with_valid_names:
            efo_ids = []
            for trait_name in trait.all_names:
                efo_ids.extend(
                    EfoMappedClinVarTrait.format_efo_id(efo_id)
                    for efo_id, efo_label in self.string_to_efo_mappings.get(trait_name.lower(), []))
            if efo_ids:
                trait.add_efo_mappings(efo_ids)
                annotated_efo_ids.update(efo_ids)

            if trait.xrefs:
                has_existing_id = True
            if trait.efo_aligned_ids:
                existing_efo_ids.update(trait.efo_aligned_ids)

        if has_existing_id:
            self.clinvar_counts['has_any_mappings'] += 1
        if existing_efo_ids:
            self.clinvar_counts['has_efo_mappings'] += 1
        if annotated_efo_ids:
            self.cmat_counts['has_efo_mappings'] += 1

        if annotated_efo_ids and existing_efo_ids:
            # TODO potentially will need to match these, e.g. with OLS
            self.f1_scores['traits'] += self.f1_score(existing_efo_ids, annotated_efo_ids)
            self.match_counts['traits'] += 1

    def report(self):
        print('\nOverall counts:')
        self.print_counter(self.overall_counts)
        print('\nClinVar counts:')
        self.print_counter(self.clinvar_counts)
        print('\nCMAT counts:')
        self.print_counter(self.cmat_counts)
        print('\nRecords with both ClinVar and CMAT annotations present:')
        self.print_counter(self.match_counts)
        print('\nF1 scores (averaged over records with both annotations present):')
        self.print_counter(self.f1_scores)
        print()

    @staticmethod
    def print_counter(counter):
        max_len = len(max(counter.keys(), key=lambda x: len(x)))
        for k, v in counter.items():
            print(f'{k: <{max_len}} {v}')

    @staticmethod
    def f1_score(ground_truth_set, predicted_set):
        if len(ground_truth_set) == 0 and len(predicted_set) == 0:
            return 0
        tp = len(predicted_set & ground_truth_set)
        fp = len(predicted_set - ground_truth_set)
        fn = len(ground_truth_set - predicted_set)
        return 2*tp / (2*tp + fp + fn)


class AnnotatedClinVarRecord(ClinVarRecord):

    def __init__(self, rcv):
        super().__init__(rcv, trait_class=EfoMappedClinVarTrait, measure_class=EnsemblAnnotatedClinVarMeasure)


class EfoMappedClinVarTrait(ClinVarTrait):

    def add_efo_mappings(self, efo_ids):
        efo_elts = []
        for efo_id in efo_ids:
            efo_id = self.format_efo_id(efo_id)
            efo_elts.append(ET.Element('XRef', attrib={'ID': efo_id, 'DB': 'EFO', 'providedBy': PROCESSOR}))
        self.trait_xml.extend(efo_elts)

    @staticmethod
    def format_efo_id(efo_id):
        if efo_id.startswith('http'):
            return efo_id.split('/')[-1].replace('_', ':')
        return efo_id


class EnsemblAnnotatedClinVarMeasure(ClinVarRecordMeasure):

    def add_ensembl_annotations(self, consequences):
        consequence_elts = []
        for consequence_attributes in consequences:
            attr_set_elt = ET.Element('AttributeSet')
            attribute_elt = ET.Element('Attribute', attrib={'Type': 'MolecularConsequence', 'providedBy': PROCESSOR})
            attribute_elt.text = consequence_attributes.so_term.so_name.replace('_', ' ')
            so_elt = ET.Element('XRef', attrib={'ID': self.format_so_term(consequence_attributes.so_term),
                                                'DB': 'Sequence Ontology'})
            ensembl_elt = ET.Element('XRef', attrib={'ID': consequence_attributes.ensembl_gene_id, 'DB': 'Ensembl'})
            attr_set_elt.extend((attribute_elt, so_elt, ensembl_elt))
            consequence_elts.append(attr_set_elt)
        self.measure_xml.extend(consequence_elts)

    @staticmethod
    def format_so_term(so_term):
        return so_term.accession.replace('_', ':')


def generate_annotated_clinvar_xml(clinvar_xml_file, efo_mapping_file, gene_mapping_file, output_xml_file):
    """Generate an annotated XML file of ClinVar RCVs based on EFO mappings file and gene mapping file (as documented in
    clinvar_to_evidence_strings)."""
    string_to_efo_mappings = load_efo_mapping(efo_mapping_file)
    variant_to_gene_mappings = CT.process_consequence_type_file(gene_mapping_file)

    dataset = AnnotatingClinVarDataset(clinvar_xml_file, string_to_efo_mappings, variant_to_gene_mappings)
    dataset.write(output_xml_file)
    dataset.report()
