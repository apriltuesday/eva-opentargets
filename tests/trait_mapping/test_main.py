"""Tests for the trait mapping pipeline. Test resources are compressed XML files which contain one or a few records
manually extracted from the main ClinVar XML to check specific cases."""
import csv
import os
import tempfile

import pytest

from cmat.trait_mapping.main import parse_traits, process_traits


def get_test_resource(resource_name):
    """Gets full path to the test resource located in the same directory as the test module."""

    # Full path to this module.
    this_module = os.path.abspath(__file__)

    # Full path to the directory where it is contained.
    module_directory = os.path.dirname(this_module)

    # Full path to the requested resource.
    return os.path.join(module_directory, 'resources', resource_name)


def run_pipeline(resource_name):
    """Runs the pipeline on a given test resource and returns the output traits, automated mappings, and curation terms
    as lists of lists."""
    input_filename = get_test_resource(resource_name)
    traits_file, mappings_file, curation_file = [tempfile.NamedTemporaryFile(delete=False) for _ in range(3)]
    filters = {
        'ontologies': 'efo,ordo,hp,mondo',
        'required': 'cttv,eva-clinvar,clinvar-xrefs,gwas',
        'preferred': 'eva-clinvar,cttv,gwas,clinvar-xrefs',
    }
    parse_traits(
        input_filepath=input_filename,
        output_traits_filepath=traits_file.name,
    )
    process_traits(
        traits_filepath=traits_file.name,
        output_mappings_filepath=mappings_file.name,
        output_curation_filepath=curation_file.name,
        filters=filters,
        zooma_host='https://www.ebi.ac.uk',
        oxo_target_list=['Orphanet', 'efo', 'hp', 'mondo'],
        oxo_distance=3,
        ontology='EFO'
    )
    output_traits = [row for row in csv.reader(open(traits_file.name), delimiter=',')]
    output_mappings = [line.rstrip().split('\t') for line in open(mappings_file.name).read().splitlines()]
    output_curation = [line.rstrip().split('\t') for line in open(curation_file.name).read().splitlines()]
    for temp_file in (traits_file, mappings_file, curation_file):
        os.remove(temp_file.name)
    return output_traits, output_mappings, output_curation


@pytest.mark.integration
@pytest.mark.skip('Redundant with new end-to-end tests, to be replaced with unit tests of trait_mapping.main')
def test_main():
    """Basic sanity test of output files, using a random sample of records."""
    output_traits, output_mappings, output_curation = run_pipeline('sample.xml.gz')
    all_terms = {x[0] for x in output_traits}
    mapped_terms = {x[0] for x in output_mappings}
    curation_terms = {x[0] for x in output_curation}
    assert len(mapped_terms) + len(curation_terms) == len(all_terms)
