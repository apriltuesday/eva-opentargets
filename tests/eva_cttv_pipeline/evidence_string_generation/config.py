import os

from eva_cttv_pipeline.clinvar_xml_io import clinvar_xml_io

OT_SCHEMA_VERSION = "2.2.6"

test_dir = os.path.dirname(__file__)
efo_mapping_file = os.path.join(test_dir, 'resources', 'string_to_ontology_mappings.tsv')
snp_2_gene_file = os.path.join(test_dir, 'resources/snp2gene_extract.tsv')


def get_expected_evidence_string(filename):
    return open(os.path.join(test_dir, 'resources', f'expected_{filename}_evidence_string.json')).read()


def get_test_clinvar_record(filename='test_clinvar_record.xml.gz'):
    """The default test file contains an extract of ClinVar XML for the record RCV000002127."""
    test_clinvar_record_file = os.path.join(test_dir, 'resources', filename)
    return [r for r in clinvar_xml_io.ClinVarDataset(test_clinvar_record_file)][0]
