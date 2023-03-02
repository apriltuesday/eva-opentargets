import gzip
import os

from eva_cttv_pipeline.evidence_string_generation.annotated_clinvar import generate_annotated_clinvar_xml

resources_dir = os.path.join(os.path.dirname(__file__), 'resources')


def test_generate_annotated_xml():
    input_file = os.path.join(resources_dir, 'test_annotation_input.xml.gz')
    efo_mapping_file = os.path.join(resources_dir, 'string_to_ontology_mappings.tsv')
    gene_mapping_file = os.path.join(resources_dir, 'snp2gene_extract.tsv')
    output_file = os.path.join(resources_dir, 'test_output.xml.gz')
    expected_output_file = os.path.join(resources_dir, 'expected_annotation_output.xml.gz')

    generate_annotated_clinvar_xml(input_file, efo_mapping_file, gene_mapping_file, output_file)

    expected_output = gzip.open(expected_output_file).read()
    actual_output = gzip.open(output_file).read()
    assert expected_output == actual_output

    if os.path.exists(output_file):
        os.remove(output_file)
