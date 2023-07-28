import gzip
import os
import re

from cmat.output_generation.annotated_clinvar import generate_annotated_clinvar_xml, string_to_set

resources_dir = os.path.join(os.path.dirname(__file__), 'resources')


def test_string_to_set():
    assert string_to_set("{'HP:0002269'}") == {'HP:0002269'}
    assert string_to_set("{'Orphanet:49382', 'MONDO:0018852'}") == {'Orphanet:49382', 'MONDO:0018852'}
    assert string_to_set('{}') == set()


def test_generate_annotated_xml():
    input_file = os.path.join(resources_dir, 'test_annotation_input.xml.gz')
    efo_mapping_file = os.path.join(resources_dir, 'string_to_ontology_mappings.tsv')
    gene_mapping_file = os.path.join(resources_dir, 'snp2gene_extract.tsv')
    output_file = os.path.join(resources_dir, 'test_output.xml.gz')
    expected_output_file = os.path.join(resources_dir, 'expected_annotation_output.xml.gz')

    generate_annotated_clinvar_xml(input_file, efo_mapping_file, gene_mapping_file, output_file)

    expected_output = gzip.open(expected_output_file, 'rt').read()
    actual_output = gzip.open(output_file, 'rt').read()

    # Ignore last processed date
    last_processed_regex = re.compile('(?<=LastProcessed=")[-0-9]+(?=")')
    expected_output = last_processed_regex.sub('', expected_output)
    actual_output = last_processed_regex.sub('', actual_output)
    assert expected_output == actual_output

    if os.path.exists(output_file):
        os.remove(output_file)
