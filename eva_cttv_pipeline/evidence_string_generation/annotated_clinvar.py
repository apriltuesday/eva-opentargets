import re
import xml.etree.ElementTree as ET
from collections import defaultdict

from eva_cttv_pipeline.clinvar_xml_io.clinvar_xml_io.ontology_uri import OntologyUri

from eva_cttv_pipeline.clinvar_xml_io.clinvar_xml_io import ClinVarTrait, ClinVarRecordMeasure, ClinVarDataset, \
    ClinVarRecord
from eva_cttv_pipeline.clinvar_xml_io.clinvar_xml_io.xml_parsing import iterate_rcv_from_xml
from eva_cttv_pipeline.evidence_string_generation.clinvar_to_evidence_strings import load_efo_mapping, \
    get_consequence_types
from eva_cttv_pipeline.evidence_string_generation import consequence_type as CT
from eva_cttv_pipeline.evidence_string_generation.set_metrics import SetComparisonMetrics

PROCESSOR = 'CMAT'


class AnnotatingClinVarDataset(ClinVarDataset):
    """This class provides the ability to parse ClinVar records (RCVs) and annotate them with EFO mappings and
    consequence mappings on the fly."""

    def __init__(self, clinvar_xml, string_to_efo_mappings, variant_to_gene_mappings,
                 eval_gene_mappings=None, eval_xref_mappings=None):
        super().__init__(clinvar_xml)
        self.header_attr['ProcessedBy'] = PROCESSOR
        self.string_to_efo_mappings = string_to_efo_mappings
        self.variant_to_gene_mappings = variant_to_gene_mappings
        self.eval_gene_mappings = eval_gene_mappings
        self.eval_xref_mappings = eval_xref_mappings
        self.overall_counts = {}
        self.gene_metrics = None
        self.conseq_metrics = None
        self.trait_metrics = None

    def __iter__(self):
        # Initialise counts
        self.overall_counts = {
            'total': 0,
            'has_supported_measure': 0,
            'has_supported_trait': 0,
        }
        self.gene_metrics = SetComparisonMetrics()
        self.conseq_metrics = SetComparisonMetrics()
        self.trait_metrics = SetComparisonMetrics()

        for rcv in iterate_rcv_from_xml(self.clinvar_xml):
            record = AnnotatedClinVarRecord(rcv)
            self.annotate(record)
            yield record

        # Finalise - computes averages, etc.
        self.gene_metrics.finalise()
        self.conseq_metrics.finalise()
        self.trait_metrics.finalise()

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

        if self.eval_gene_mappings:
            existing_ensembl_ids = self.eval_gene_mappings.get(record.accession, [])
            self.gene_metrics.count_and_score(cv_set=existing_ensembl_ids, cmat_set=annotated_genes)
            self.conseq_metrics.count_and_score(cv_set=record.measure.existing_so_terms, cmat_set=annotated_conseqs)

    def annotate_and_count_traits(self, record):
        for trait in record.traits_with_valid_names:
            existing_efo_ids = set()
            annotated_efo_ids = set()
            efo_ids = []
            for trait_name in trait.all_names:
                efo_ids.extend(
                    EfoMappedClinVarTrait.format_efo_id(efo_id)
                    for efo_id, efo_label in self.string_to_efo_mappings.get(trait_name.lower(), []))

            for db, iden, _ in trait.current_efo_aligned_xrefs:
                curie = OntologyUri(iden, db).curie
                if curie:
                    existing_efo_ids.add(curie)
            trait.add_efo_mappings(efo_ids)
            if self.eval_xref_mappings:
                for efo_id in efo_ids:
                    # Attempt to match to an id in ClinVar based on synonyms
                    for cv_id in existing_efo_ids:
                        if efo_id in self.eval_xref_mappings.get(cv_id, {}):
                            annotated_efo_ids.add(cv_id)
                    # If didn't find anything, just use our id
                    if not annotated_efo_ids:
                        annotated_efo_ids.add(efo_id)

                self.trait_metrics.count_and_score(cv_set=existing_efo_ids, cmat_set=annotated_efo_ids)

    def report(self):
        print('\nOverall counts (RCVs):')
        self.print_counter(self.overall_counts)
        if self.eval_gene_mappings:
            print('\nGene annotations:')
            self.gene_metrics.report()
            print('\nFunctional consequences:')
            self.conseq_metrics.report()
        if self.eval_xref_mappings:
            print('\nTrait mappings:')
            self.trait_metrics.report()
        print()

    @staticmethod
    def print_counter(counter):
        max_len = len(max(counter.keys(), key=lambda x: len(x)))
        for k, v in counter.items():
            print(f'{k: <{max_len}} {v}')


class AnnotatedClinVarRecord(ClinVarRecord):

    def __init__(self, rcv):
        super().__init__(rcv, trait_class=EfoMappedClinVarTrait, measure_class=EnsemblAnnotatedClinVarMeasure)


class EfoMappedClinVarTrait(ClinVarTrait):

    def add_efo_mappings(self, efo_ids):
        efo_elts = []
        for efo_id in efo_ids:
            efo_id = self.format_efo_id(efo_id)
            # Include Status attribute so this isn't included among current xrefs
            efo_elts.append(ET.Element('XRef', attrib={
                'ID': efo_id, 'DB': 'EFO', 'Status': 'annotated', 'providedBy': PROCESSOR}))
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


def load_evaluation_gene_mappings(input_path):
    mapping = defaultdict(list)
    with open(input_path) as input_file:
        for line in input_file:
            cols = line.strip().split('\t')
            if len(cols) != 2:
                continue
            mapping[cols[0]].append(cols[1])
    return mapping


def load_evaluation_xref_mappings(input_path):
    mapping = {}
    with open(input_path) as input_file:
        for line in input_file:
            cols = line.strip().split('\t')
            if len(cols) != 2:
                continue
            mapping[cols[0]] = re.sub(r"{|}|'", '', cols[1]).split(', ')
    return mapping


def generate_annotated_clinvar_xml(clinvar_xml_file, efo_mapping_file, gene_mapping_file, output_xml_file,
                                   eval_gene_file=None, eval_xref_file=None):
    """Generate an annotated XML file of ClinVar RCVs based on EFO mappings file and gene mapping file (as documented in
    clinvar_to_evidence_strings)."""
    string_to_efo_mappings = load_efo_mapping(efo_mapping_file)
    variant_to_gene_mappings = CT.process_consequence_type_file(gene_mapping_file)
    # Need both files to do an evaluation
    if eval_gene_file and eval_xref_file:
        eval_gene_mappings = load_evaluation_gene_mappings(eval_gene_file)
        eval_xref_mappings = load_evaluation_xref_mappings(eval_xref_file)
        dataset = AnnotatingClinVarDataset(clinvar_xml_file, string_to_efo_mappings, variant_to_gene_mappings,
                                           eval_gene_mappings=eval_gene_mappings,
                                           eval_xref_mappings=eval_xref_mappings)
    else:
        dataset = AnnotatingClinVarDataset(clinvar_xml_file, string_to_efo_mappings, variant_to_gene_mappings)
    dataset.write(output_xml_file)
    dataset.report()
