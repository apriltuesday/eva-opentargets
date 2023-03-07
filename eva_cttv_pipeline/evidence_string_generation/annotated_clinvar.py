import xml.etree.ElementTree as ET

from eva_cttv_pipeline.clinvar_xml_io.clinvar_xml_io import ClinVarTrait, ClinVarRecordMeasure, ClinVarDataset, \
    ClinVarRecord
from eva_cttv_pipeline.clinvar_xml_io.clinvar_xml_io.xml_parsing import iterate_rcv_from_xml
from eva_cttv_pipeline.evidence_string_generation.clinvar_to_evidence_strings import load_efo_mapping, \
    get_consequence_types
from eva_cttv_pipeline.evidence_string_generation import consequence_type as CT


PROVIDER = 'ClinVarXmlPipeline'


class AnnotatingClinVarDataset(ClinVarDataset):
    """This class provides the ability to parse ClinVar records (RCVs) and annotate them with EFO mappings and
    consequence mappings on the fly."""

    def __init__(self, clinvar_xml, string_to_efo_mappings, variant_to_gene_mappings):
        super().__init__(clinvar_xml)
        self.header = f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<ReleaseSet Dated="{self.today()}" ProcessedBy="{PROVIDER}" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" Type="full" xsi:noNamespaceSchemaLocation="http://ftp.ncbi.nlm.nih.gov/pub/clinvar/xsd_public/clinvar_public_1.60.xsd">
'''
        self.header = self.header.encode('utf-8')
        self.string_to_efo_mappings = string_to_efo_mappings
        self.variant_to_gene_mappings = variant_to_gene_mappings

    def __iter__(self):
        for rcv in iterate_rcv_from_xml(self.clinvar_xml):
            record = AnnotatedClinVarRecord(rcv)
            self.annotate(record)
            yield record

    def annotate(self, record):
        # Functional consequences for measure
        if record.measure:
            consequence_types = get_consequence_types(record.measure, self.variant_to_gene_mappings)
            record.measure.add_ensembl_annotations(consequence_types)
        # EFO terms for trait
        for trait in record.traits_with_valid_names:
            efo_ids = []
            for trait_name in trait.all_names:
                efo_ids.extend(efo_id for efo_id, efo_label in self.string_to_efo_mappings.get(trait_name.lower(), []))
            trait.add_efo_mappings(efo_ids)


class AnnotatedClinVarRecord(ClinVarRecord):

    def __init__(self, rcv):
        super().__init__(rcv, trait_class=EfoMappedClinVarTrait, measure_class=EnsemblAnnotatedClinVarMeasure)


class EfoMappedClinVarTrait(ClinVarTrait):

    def add_efo_mappings(self, efo_ids):
        efo_elts = []
        for efo_id in efo_ids:
            if efo_id.startswith('http'):
                efo_id = efo_id.split('/')[-1].replace('_', ':')
            efo_elts.append(ET.Element('XRef', attrib={'ID': efo_id, 'DB': 'EFO', 'providedBy': PROVIDER}))
        self.trait_xml.extend(efo_elts)


class EnsemblAnnotatedClinVarMeasure(ClinVarRecordMeasure):

    def add_ensembl_annotations(self, consequences):
        consequence_elts = []
        for consequence_attributes in consequences:
            attr_set_elt = ET.Element('AttributeSet')
            attribute_elt = ET.Element('Attribute', attrib={'Type': 'MolecularConsequence', 'providedBy': PROVIDER})
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
