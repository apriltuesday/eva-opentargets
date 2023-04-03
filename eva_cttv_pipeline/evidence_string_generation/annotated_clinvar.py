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
        self.record_counts = {}

    def __iter__(self):
        # Initialise counts
        self.record_counts = {
            'total': 0,
            'has_supported_measure': 0,
            'has_consequences': 0,
            'has_supported_trait': 0,
            'has_efo_mappings': 0
        }
        for rcv in iterate_rcv_from_xml(self.clinvar_xml):
            record = AnnotatedClinVarRecord(rcv)
            self.annotate(record)
            yield record

    def annotate(self, record):
        self.record_counts['total'] += 1

        # Functional consequences for measure
        if record.measure:
            self.record_counts['has_supported_measure'] += 1
            consequence_types = get_consequence_types(record.measure, self.variant_to_gene_mappings)
            if consequence_types:
                self.record_counts['has_consequences'] += 1
                record.measure.add_ensembl_annotations(consequence_types)

        # EFO terms for trait
        if record.traits_with_valid_names:
            self.record_counts['has_supported_trait'] += 1
        record_has_efo = False
        for trait in record.traits_with_valid_names:
            efo_ids = []
            for trait_name in trait.all_names:
                efo_ids.extend(efo_id for efo_id, efo_label in self.string_to_efo_mappings.get(trait_name.lower(), []))
            if efo_ids:
                record_has_efo = True
                trait.add_efo_mappings(efo_ids)
        if record_has_efo:
            self.record_counts['has_efo_mappings'] += 1

    def report(self):
        print('\nRecord counts:')
        for key, val in self.record_counts.items():
            print(f'{key: <21} {val}')
        print()


class AnnotatedClinVarRecord(ClinVarRecord):

    def __init__(self, rcv):
        super().__init__(rcv, trait_class=EfoMappedClinVarTrait, measure_class=EnsemblAnnotatedClinVarMeasure)


class EfoMappedClinVarTrait(ClinVarTrait):

    def add_efo_mappings(self, efo_ids):
        efo_elts = []
        for efo_id in efo_ids:
            if efo_id.startswith('http'):
                efo_id = efo_id.split('/')[-1].replace('_', ':')
            efo_elts.append(ET.Element('XRef', attrib={'ID': efo_id, 'DB': 'EFO', 'providedBy': PROCESSOR}))
        self.trait_xml.extend(efo_elts)


class EnsemblAnnotatedClinVarMeasure(ClinVarRecordMeasure):

    def add_ensembl_annotations(self, consequences):
        consequence_elts = []
        for consequence_attributes in consequences:
            attr_set_elt = ET.Element('AttributeSet')
            attribute_elt = ET.Element('Attribute', attrib={'Type': 'MolecularConsequence', 'providedBy': PROCESSOR})
            attribute_elt.text = consequence_attributes.so_term.so_name.replace('_', ' ')
            so_elt = ET.Element('XRef', attrib={'ID': consequence_attributes.so_term.accession.replace('_', ':'),
                                                'DB': 'Sequence Ontology'})
            ensembl_elt = ET.Element('XRef', attrib={'ID': consequence_attributes.ensembl_gene_id, 'DB': 'Ensembl'})
            attr_set_elt.extend((attribute_elt, so_elt, ensembl_elt))
            consequence_elts.append(attr_set_elt)
        self.measure_xml.extend(consequence_elts)


def generate_annotated_clinvar_xml(clinvar_xml_file, efo_mapping_file, gene_mapping_file, output_xml_file):
    """Generate an annotated XML file of ClinVar RCVs based on EFO mappings file and gene mapping file (as documented in
    clinvar_to_evidence_strings)."""
    string_to_efo_mappings = load_efo_mapping(efo_mapping_file)
    variant_to_gene_mappings = CT.process_consequence_type_file(gene_mapping_file)

    dataset = AnnotatingClinVarDataset(clinvar_xml_file, string_to_efo_mappings, variant_to_gene_mappings)
    dataset.write(output_xml_file)
    dataset.report()
